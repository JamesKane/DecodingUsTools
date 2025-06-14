use bio::io::fasta;

use crate::haplogroup::types::Locus;
use indicatif::ProgressBar;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;

pub fn collect_snp_calls(
    min_depth: u32,
    min_quality: u8,
    fasta_reader: &mut fasta::IndexedReader<File>,
    mut bam: bam::IndexedReader,
    build_id: String,
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
    let header = bam.header().clone();

    // Create fetch definition for the entire chromosome
    let fetch_def = format!("{}:1-", chromosome);
    bam.fetch(&fetch_def)?;

    process_region(
        &mut bam,
        &header,
        fasta_reader,
        positions,
        &build_id,
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
    header: &bam::HeaderView,
    fasta: &mut fasta::IndexedReader<File>,
    positions: &HashMap<u32, Vec<(&str, &Locus)>>,
    build_id: &str,
    min_depth: u32,
    min_quality: u8,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
    progress: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut coverage: HashMap<u32, Vec<u8>> = HashMap::new();

    for r in bam.records() {
        let record = r?;
        let start_pos = record.pos() as u32;
        progress.set_position(start_pos as u64);

        if record.mapq() >= min_quality {
            let sequence = record.seq().as_bytes();
            let cigar = record.cigar();
            let mut ref_pos = start_pos;
            let mut read_pos = 0;

            let tid = record.tid();
            let ref_name = String::from_utf8_lossy(header.tid2name(tid as u32));

            for op in cigar.iter() {
                match op.char() {
                    'M' | '=' | 'X' => {
                        let len = op.len() as usize;
                        for i in 0..len {
                            let vcf_pos = ref_pos + 1;
                            // Check if this position has any loci and if they match the current chromosome
                            if let Some(loci) = positions.get(&vcf_pos) {
                                let relevant_loci = loci.iter().any(|(_, locus)| {
                                    if let Some(coord) = locus.coordinates.get(build_id) {
                                        coord.chromosome == ref_name
                                    } else {
                                        false
                                    }
                                });

                                if relevant_loci && read_pos + i < sequence.len() {
                                    let base_index = read_pos + i;
                                    let base = sequence[base_index].to_ascii_uppercase();

                                    if let Ok(_) =
                                        fasta.fetch(&ref_name, ref_pos as u64, (ref_pos + 1) as u64)
                                    {
                                        let mut ref_seq = Vec::new();
                                        fasta.read(&mut ref_seq)?;
                                        coverage.entry(vcf_pos).or_default().push(base);
                                    }
                                }
                            }
                            ref_pos += 1;
                        }
                        read_pos += len;
                    }
                    'D' | 'N' => {
                        ref_pos += op.len();
                    }
                    'I' | 'S' => {
                        read_pos += op.len() as usize;
                    }
                    _ => {}
                }
            }
        }
    }

    for (pos, bases) in coverage {
        if bases.len() >= min_depth as usize {
            let mut base_counts: HashMap<char, u32> = HashMap::new();

            for &base in &bases {
                *base_counts.entry(base as char).or_insert(0) += 1;
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
