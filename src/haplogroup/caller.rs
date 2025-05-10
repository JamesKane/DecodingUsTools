use crate::haplogroup::types::Snp;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;

pub(crate) fn collect_snp_calls(
    min_depth: u32,
    min_quality: u8,
    mut fasta_reader: &mut bio::io::fasta::IndexedReader<File>,
    mut bam: IndexedReader,
    positions: &mut HashMap<u32, Vec<(&str, &Snp)>>,
    mut snp_calls: &mut HashMap<u32, (char, u32, f64)>,
) -> Result<(), Box<dyn Error>> {
    // Create a progress bar for BAM processing
    let progress_style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
        .unwrap()
        .progress_chars("#>-");

    let header = bam.header().clone();

    // Determine which chromosome we need to process
    let need_y = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "chrY"));
    let need_mt = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "chrM"));

    // Process one chromosome based on the tree type
    if need_y {
        let tid = bam.header().tid(b"chrY").ok_or("chrY not found in BAM")?;
        let total_len = bam.header().target_len(tid).unwrap_or(0);

        let progress_bar = ProgressBar::new(total_len);
        progress_bar.set_style(progress_style.clone());
        progress_bar.set_message("Processing chromosome Y...");

        // Create the region string for fetch
        let region = format!("chrY:1-{}", total_len);
        bam.fetch(&region)?;
        process_region(
            &mut bam,
            &header,
            &mut fasta_reader,
            &positions,
            min_depth,
            min_quality,
            &mut snp_calls,
            &progress_bar,
        )?;

        progress_bar.finish_with_message("Chromosome Y processing complete");
    } else if need_mt {
        let tid = bam.header().tid(b"chrM").ok_or("chrM not found in BAM")?;
        let total_len = bam.header().target_len(tid).unwrap_or(0);

        let progress_bar = ProgressBar::new(total_len);
        progress_bar.set_style(progress_style);
        progress_bar.set_message("Processing mitochondrial DNA...");

        // Create the region string for fetch
        let region = format!("chrM:1-{}", total_len);
        bam.fetch(&region)?;
        process_region(
            &mut bam,
            &header,
            &mut fasta_reader,
            &positions,
            min_depth,
            min_quality,
            &mut snp_calls,
            &progress_bar,
        )?;

        progress_bar.finish_with_message("Mitochondrial DNA processing complete");
    }
    Ok(())
}

fn process_region<R: Read>(
    bam: &mut R,
    header: &rust_htslib::bam::HeaderView,
    fasta: &mut bio::io::fasta::IndexedReader<File>,
    positions: &HashMap<u32, Vec<(&str, &Snp)>>,
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
                            if positions.contains_key(&vcf_pos) && read_pos + i < sequence.len() {
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
