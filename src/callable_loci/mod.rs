mod options;
pub(crate) mod profilers;
pub(crate) mod report;
pub(crate) mod types;
mod utils;

use crate::callable_loci::types::CalledState;
pub use options::CallableOptions;
use profilers::{
    bam_stats::BamStats, callable_profiler::CallableProfiler, contig_profiler::ContigProfiler,
};
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::Read;
use rust_htslib::{bam, faidx};
use std::collections::HashMap;
use std::error::Error;

fn process_position(pileup: &bam::pileup::Pileup, options: &CallableOptions) -> (u32, u32, u32) {
    let mut raw_depth = 0;
    let mut qc_depth = 0;
    let mut low_mapq_count = 0;

    for align in pileup.alignments() {
        raw_depth += 1;
        let record = align.record();

        if record.mapq() <= options.max_low_mapq {
            low_mapq_count += 1;
        }

        if record.mapq() >= options.min_mapping_quality {
            if let Some(qpos) = align.qpos() {
                if let Some(&qual) = record.qual().get(qpos) {
                    if qual >= options.min_base_quality || align.is_del() {
                        qc_depth += 1;
                    }
                }
            }
        }
    }

    (raw_depth, qc_depth, low_mapq_count)
}

pub fn process_single_contig(
    bam: &mut bam::IndexedReader,
    fasta: &mut faidx::Reader,
    header: &bam::HeaderView,
    counter: &mut CallableProfiler,
    contig_stats: &mut HashMap<usize, ContigProfiler>,
    options: &CallableOptions,
    tid: usize,
) -> Result<(), Box<dyn Error>> {
    let contig_len = header.target_len(tid as u32).unwrap_or(0);
    bam.fetch((tid as u32, 0, contig_len))?;
    let mut pileup = bam.pileup();
    pileup.set_max_depth(if options.max_depth > 0 {
        options.max_depth
    } else {
        500
    });

    let contig = std::str::from_utf8(header.tid2name(tid as u32))?;
    let mut current_pos = 0;

    for p in pileup {
        let pileup = p?;
        if pileup.tid() as usize != tid {
            break;
        }

        let pos = pileup.pos();

        // Process any gaps (NO_COVERAGE positions) before this pileup position
        while current_pos < pos {
            if let Some(stats) = contig_stats.get_mut(&tid) {
                stats.progress_bar.set_position(current_pos as u64);
            }

            let seq = fasta.fetch_seq(contig, current_pos as usize, current_pos as usize)?;
            let ref_base = seq.first().copied().unwrap_or(b'N');

            counter.process_position(
                contig,
                current_pos,
                ref_base,
                0, // raw_depth
                0, // qc_depth
                0, // low_mapq_count
                options,
            )?;

            current_pos += 1;
        }

        // Process the current pileup position
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.progress_bar.set_position(pos as u64);
        }

        let seq = fasta.fetch_seq(contig, pos as usize, pos as usize)?;
        let ref_base = seq.first().copied().unwrap_or(b'N');

        let (raw_depth, qc_depth, low_mapq_count) = process_position(&pileup, options);

        counter.process_position(
            contig,
            pos,
            ref_base,
            raw_depth,
            qc_depth,
            low_mapq_count,
            options,
        )?;

        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.process_position(raw_depth, pileup.alignments());
        }

        current_pos = pos + 1;
    }

    // Process any remaining positions after the last pileup
    while current_pos < contig_len as u32 {
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.progress_bar.set_position(current_pos as u64);
        }

        let seq = fasta.fetch_seq(contig, current_pos as usize, current_pos as usize)?;
        let ref_base = seq.first().copied().unwrap_or(b'N');

        counter.process_position(
            contig,
            current_pos,
            ref_base,
            0, // raw_depth
            0, // qc_depth
            0, // low_mapq_count
            options,
        )?;

        current_pos += 1;
    }

    let stats = &contig_stats[&tid];
    counter.finish_contig(&stats.name, stats.length as u32)?;
    Ok(())
}

fn detect_aligner(header: &HeaderView) -> String {
    let header_text = String::from_utf8_lossy(header.as_bytes()).to_string();
    let header_lower = header_text.to_lowercase();

    // Check for common aligners in header
    if header_lower.contains("@pg\tid:bwa") {
        "BWA".to_string()
    } else if header_lower.contains("@pg\tid:minimap2") {
        "minimap2".to_string()
    } else if header_lower.contains("@pg\tid:bowtie2") {
        "Bowtie2".to_string()
    } else if header_lower.contains("@pg\tid:star") {
        "STAR".to_string()
    } else if header_lower.contains("bwa") {
        "BWA".to_string()
    } else if header_lower.contains("minimap2") {
        "minimap2".to_string()
    } else if header_lower.contains("bowtie2") {
        "Bowtie2".to_string()
    } else if header_lower.contains("star") {
        "STAR".to_string()
    } else {
        "Unknown".to_string()
    }
}
