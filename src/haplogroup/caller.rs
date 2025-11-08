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
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        indicatif::ProgressStyle::default_spinner()
            .template("{spinner:.green} Processing position {pos}")
            .unwrap(),
    );

    let target_tid = header
        .tid(chromosome.as_bytes())
        .ok_or_else(|| format!("Chromosome {} not found in BAM header", chromosome))? as i32;

    process_region(
        bam,
        target_tid,
        positions,
        min_depth,
        min_quality,
        snp_calls,
        &progress,
    )?;

    progress.finish_and_clear();
    Ok(())
}


fn process_region<R: Read>(
    bam: &mut R,
    target_tid: i32,
    positions: &HashMap<u32, Vec<(&str, &Locus)>>,
    min_depth: u32,
    min_quality: u8,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
    progress: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut coverage: HashMap<u32, Vec<char>> = HashMap::new();

    for r in bam.records() {
        let record = r?;
        if record.tid() != target_tid { continue; }
        let start_pos = record.pos() as u32;
        progress.set_position(start_pos as u64);

        if record.mapq() >= min_quality {
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
                            if let Some(_loci) = positions.get(&vcf_pos) {
                                let rp = read_pos + i;
                                if rp < read_len {
                                    // Filter by per-base quality and skip 'N'
                                    let base_qual = record.qual()[rp] as u8;
                                    if base_qual < min_quality { continue; }
                                    let enc = record.seq().encoded_base(rp);
                                    if enc == 15 { continue; }
                                    let base = match enc { 1 => 'A', 2 => 'C', 4 => 'G', 8 => 'T', 15 => 'N', _ => 'N' };
                                    // Do NOT reverse-complement per-read bases; alignment already orients them
                                    coverage.entry(vcf_pos).or_default().push(base);
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
    }

    for (pos, bases) in coverage {
        if bases.len() >= min_depth as usize {
            let mut base_counts: HashMap<char, u32> = HashMap::new();

            for base in &bases {
                *base_counts.entry(*base).or_insert(0) += 1;
            }

            if let Some((&base, &count)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                let total = bases.len() as u32;
                let freq = count as f64 / total as f64;

                if freq >= 0.7 {
                    snp_calls.insert(pos, (base, total, freq));
                }
            }
        }
    }

    Ok(())
}
