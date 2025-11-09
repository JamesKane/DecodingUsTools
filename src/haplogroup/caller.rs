use crate::haplogroup::types::Locus;
use indicatif::ProgressBar;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::error::Error;
use std::collections::HashSet;

pub type BaseCounts = HashMap<char, u32>;

/// Tally base counts at the specified 1-based positions on a given chromosome/contig.
/// - Positions are provided as an iterator of u32 (1-based genomic positions)
/// - Returns a map of position -> base counts (A/C/G/T/N after reverse-complement correction)
pub fn tally_bases_at_positions<R: Read, I: IntoIterator<Item = u32>>(
    bam: &mut R,
    header: &bam::HeaderView,
    chromosome: &str,
    positions: I,
    min_mapq: u8,
) -> Result<HashMap<u32, BaseCounts>, Box<dyn Error>> {
    let target_tid = header
        .tid(chromosome.as_bytes())
        .ok_or_else(|| format!("Chromosome {} not found in BAM header", chromosome))? as i32;

    let mut wanted: HashSet<u32> = positions.into_iter().collect();
    if wanted.is_empty() {
        return Ok(HashMap::new());
    }
    let mut coverage: HashMap<u32, BaseCounts> = HashMap::new();

    for r in bam.records() {
        let record = r?;
        if record.tid() != target_tid { continue; }
        if record.mapq() < min_mapq { continue; }
        let start_pos = record.pos() as u32; // 0-based
        let read_len = record.seq_len() as usize;
        let cigar = record.cigar();
        let mut ref_pos = start_pos; // 0-based
        let mut read_pos: usize = 0;

        for op in cigar.iter() {
            match op.char() {
                'M' | '=' | 'X' => {
                    let len = op.len() as usize;
                    for i in 0..len {
                        let vcf_pos = ref_pos + 1; // 1-based
                        if wanted.contains(&vcf_pos) {
                            let rp = read_pos + i;
                            if rp < read_len {
                                // Per-base quality filter and N masking
                                let base_qual = record.qual()[rp] as u8;
                                if base_qual < min_mapq { continue; }
                                let enc = record.seq().encoded_base(rp);
                                if enc == 15 { continue; } // skip 'N'
                                let base = match enc { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', 15 => 'N', _ => 'N' };
                                // Do NOT reverse-complement: base calls are in reference orientation per CIGAR/alignment
                                *coverage.entry(vcf_pos).or_default().entry(base).or_insert(0) += 1;
                            }
                        }
                        ref_pos += 1;
                    }
                    read_pos += len;
                }
                'D' | 'N' => { ref_pos += op.len(); }
                'I' | 'S' => { read_pos += op.len() as usize; }
                _ => {}
            }
        }
    }

    Ok(coverage)
}

pub fn collect_snp_calls<R: Read>(
    min_depth: u32,
    min_quality: u8,
    bam: &mut R,
    header: &bam::HeaderView,
    _build_id: String,
    chromosome: String,
    positions: &mut HashMap<u32, Vec<(&str, &Locus)>>,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
) -> Result<(), Box<dyn Error>> {
    // Env-gated debug logging
    let debug = std::env::var("DECODINGUS_DEBUG_CALLER").ok().filter(|v| !v.is_empty() && v != "0").is_some();

    // Allow env overrides while preserving existing API behavior by default
    let env_min_mapq = std::env::var("DECODINGUS_MIN_MAPQ").ok().and_then(|s| s.parse::<u8>().ok());
    let env_min_baseq = std::env::var("DECODINGUS_MIN_BASEQ").ok().and_then(|s| s.parse::<u8>().ok());
    let env_call_freq = std::env::var("DECODINGUS_CALL_FREQ").ok().and_then(|s| s.parse::<f64>().ok());

    let min_mapq: u8 = env_min_mapq.unwrap_or(min_quality);
    let min_baseq: u8 = env_min_baseq.unwrap_or(min_quality);
    let call_freq: f64 = env_call_freq.unwrap_or(0.70);

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing position {pos}")
            .unwrap(),
    );

    let target_tid = match header.tid(chromosome.as_bytes()) {
        Some(tid) => tid as i32,
        None => {
            if debug {
                // List available reference names to aid alias debugging
                let names: Vec<String> = (0..header.target_count())
                    .filter_map(|i| std::str::from_utf8(header.tid2name(i)).ok().map(|s| s.to_string()))
                    .collect();
                eprintln!("[caller.debug] Chromosome '{}' not found. Available targets: {:?}", chromosome, names);
            }
            return Err(format!("Chromosome {} not found in BAM header", chromosome).into());
        }
    };

    if debug {
        let positions_in = positions.len();
        eprintln!(
            "[caller.debug] target='{}' tid={} positions_in={} min_depth={} min_mapq={} min_baseq={} call_freq={}",
            chromosome, target_tid, positions_in, min_depth, min_mapq, min_baseq, call_freq
        );
    }

    process_region(
        bam,
        target_tid,
        positions,
        min_depth,
        min_mapq,
        min_baseq,
        call_freq,
        snp_calls,
        &progress,
        debug,
    )?;

    progress.finish_and_clear();
    Ok(())
}


fn process_region<R: Read>(
    bam: &mut R,
    target_tid: i32,
    positions: &HashMap<u32, Vec<(&str, &Locus)>>,
    min_depth: u32,
    min_mapq: u8,
    min_baseq: u8,
    call_freq: f64,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
    progress: &ProgressBar,
    debug: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    // Debug specific position for P312 troubleshooting
    let debug_p312 = std::env::var("DECODINGUS_DEBUG_P312").ok().filter(|v| !v.is_empty() && v != "0").is_some();
    let p312_pos = 20_901_962u32; // CHM13v2.0 position for P312

    if debug_p312 {
        if positions.contains_key(&p312_pos) {
            eprintln!("[caller.P312] positions HashMap CONTAINS {}", p312_pos);
        } else {
            eprintln!("[caller.P312] positions HashMap DOES NOT CONTAIN {} !!", p312_pos);
        }
    }
    
    let mut coverage: HashMap<u32, Vec<char>> = HashMap::new();
    let mut seen_any_cov: HashSet<u32> = HashSet::new();

    // Track raw observations at P312 position if debugging
    let mut p312_raw_bases: Vec<(char, u8, bool)> = Vec::new(); // (base, qual, is_reverse_strand)

    for r in bam.records() {
        let record = r?;
        if record.tid() != target_tid { continue; }
        let start_pos = record.pos() as u32;
        progress.set_position(start_pos as u64);

        if record.mapq() >= min_mapq {
            let read_len = record.seq_len() as usize;
            let cigar = record.cigar();
            let is_reverse = record.is_reverse();
            let mut ref_pos = start_pos; // 0-based
            let mut read_pos: usize = 0;

            for op in cigar.iter() {
                match op.char() {
                    'M' | '=' | 'X' => {
                        let len = op.len() as usize;
                        for i in 0..len {
                            let vcf_pos = ref_pos + i as u32 + 1; // 1-based
                            if let Some(_loci) = positions.get(&vcf_pos) {
                                let rp = read_pos + i;
                                if rp < read_len {
                                    // Filter by per-base quality and skip 'N'
                                    let base_qual = record.qual()[rp] as u8;
                                    if base_qual < min_baseq { continue; }
                                    let enc = record.seq().encoded_base(rp);
                                    if enc == 15 { continue; }
                                    let base = match enc { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', 15 => 'N', _ => 'N' };

                                    // Debug P312 position specifically
                                    if debug_p312 && vcf_pos == p312_pos {
                                        let read_name = String::from_utf8_lossy(record.qname());
                                        eprintln!("[caller.P312] READ NAME: {}", read_name);

                                        // Get a few bases around this position for context
                                        let mut context = String::new();
                                        for offset in -2..=2i32 {
                                            let ctx_pos = (rp as i32 + offset) as usize;
                                            if ctx_pos < read_len {
                                                let ctx_enc = record.seq().encoded_base(ctx_pos);
                                                let ctx_base = match ctx_enc { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', 15 => 'N', _ => 'N' };
                                                context.push(ctx_base);
                                            }
                                        }

                                        p312_raw_bases.push((base, base_qual, is_reverse));
                                        eprintln!("[caller.P312] RAW: ref_pos={} (1-based={}), read_pos={}, base='{}', qual={}, reverse={}, enc={}, context={}",
                                                  ref_pos + i as u32, vcf_pos, rp, base, base_qual, is_reverse, enc, context);
                                    }

                                    seen_any_cov.insert(vcf_pos);
                                    coverage.entry(vcf_pos).or_default().push(base);
                                }
                            }
                        }
                        ref_pos += len as u32;
                        read_pos += len;
                    }
                    'D' | 'N' => { ref_pos += op.len(); }
                    'I' | 'S' => { read_pos += op.len() as usize; }
                    _ => {}
                }
            }
        }
    }

    // P312 debugging output before calling
    if debug_p312 && !p312_raw_bases.is_empty() {
        eprintln!("[caller.P312] Raw bases observed at pos {}:", p312_pos);
        let mut base_counts: HashMap<char, u32> = HashMap::new();
        let mut strand_info: HashMap<char, (u32, u32)> = HashMap::new(); // (forward_count, reverse_count)
        for (base, qual, is_rev) in &p312_raw_bases {
            eprintln!("[caller.P312]   base={} qual={} reverse_strand={}", base, qual, is_rev);
            *base_counts.entry(*base).or_insert(0) += 1;
            let entry = strand_info.entry(*base).or_insert((0, 0));
            if *is_rev {
                entry.1 += 1;
            } else {
                entry.0 += 1;
            }
        }
        eprintln!("[caller.P312] Base counts:");
        for (base, count) in base_counts.iter() {
            let (fwd, rev) = strand_info.get(base).unwrap_or(&(0, 0));
            eprintln!("[caller.P312]   {}: {} (forward={}, reverse={})", base, count, fwd, rev);
        }
        if let Some(bases) = coverage.get(&p312_pos) {
            eprintln!("[caller.P312] Coverage vector at {}: {:?}", p312_pos, bases);
        }
    }

    let mut passing_min_depth: Vec<u32> = Vec::new();
    let mut emitting_calls: Vec<u32> = Vec::new();

    for (pos, bases) in coverage.iter() {
        if bases.len() >= min_depth as usize {
            passing_min_depth.push(*pos);
            let mut base_counts: HashMap<char, u32> = HashMap::new();
            for base in bases.iter() { *base_counts.entry(*base).or_insert(0) += 1; }
            if let Some((_, &count)) = base_counts.iter().max_by_key(|&(_, c)| c) {
                let total = bases.len() as u32;
                let freq = count as f64 / total as f64;
                if freq >= call_freq { emitting_calls.push(*pos); }
            }
        }
    }

    // Emit calls after diagnostics collection
    for (pos, bases) in coverage.into_iter() {
        if bases.len() >= min_depth as usize {
            let mut base_counts: HashMap<char, u32> = HashMap::new();
            for base in &bases { *base_counts.entry(*base).or_insert(0) += 1; }
            if let Some((&base, &count)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                let total = bases.len() as u32;
                let freq = count as f64 / total as f64;
                if freq >= call_freq {
                    snp_calls.insert(pos, (base, total, freq));
                }
            }
        }
    }

    if debug {
        let positions_in = positions.len() as u32;
        let seen_cov = seen_any_cov.len() as u32;
        let pass_depth = passing_min_depth.len() as u32;
        let emit = emitting_calls.len() as u32;
        eprintln!(
            "[caller.debug] positions_in={} seen_any_coverage={} passing_min_depth={} emitting_calls={}",
            positions_in, seen_cov, pass_depth, emit
        );
        // Print first N example tallies
        let mut printed = 0;
        let maxn = 5;
        for pos in passing_min_depth.into_iter().take(maxn) {
            // reconstruct counts for printing
            // Note: since coverage moved, re-tally from snp_calls or skip if absent
            if let Some((base, total, freq)) = snp_calls.get(&pos) {
                eprintln!(
                    "[caller.debug] example pos={} call={} total={} freq={:.3}",
                    pos, base, total, freq
                );
            }
            printed += 1;
            if printed >= maxn { break; }
        }
    }

    Ok(())
}
