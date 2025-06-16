mod options;
mod profilers;
mod report;
mod types;
mod utils;

use crate::callable_loci::types::CalledState;
use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
pub use options::CallableOptions;
use profilers::{
    bam_stats::BamStats, callable_profiler::CallableProfiler, contig_profiler::ContigProfiler,
};
use rust_htslib::bam;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::Write;

pub fn run(
    bam_file: String,
    reference_file: String,
    output_bed: String,
    output_summary: String,
    mut options: CallableOptions,
    contigs: Option<Vec<String>>,
) -> Result<(), Box<dyn Error>> {
    // Create and collect BAM statistics first
    let mut bam_stats = BamStats::new(10000);
    bam_stats.collect_stats(&bam_file)?;

    // Print BAM statistics
    let stats = bam_stats.get_stats();
    println!("\nBAM Statistics:");
    if let Some(avg_read_length) = stats.get("average_read_length") {
        println!("Average read length: {:.1} bp", avg_read_length);
    }
    if let Some(paired_percentage) = stats.get("paired_percentage") {
        println!("Paired reads: {:.1}%", paired_percentage);
    }
    if let Some(avg_insert_size) = stats.get("average_insert_size") {
        println!("Average insert size: {:.1} bp", avg_insert_size);
    }
    println!(); // Add a blank line for better formatting

    options = options.with_contigs(contigs);

    let mut bam = bam::IndexedReader::from_path(&bam_file)?;
    let header = bam.header().clone();
    let mut fasta = IndexedReader::from_file(&reference_file)?;

    let (multi_progress, main_progress) = setup_progress_bars();
    let mut contig_stats = initialize_contig_stats(&header, &options, &multi_progress)?;
    validate_contig_selection(&contig_stats, &options)?;

    let mut counter = initialize_counter(&contig_stats, &output_bed)?;
    process_contigs(
        &mut bam,
        &mut fasta,
        &header,
        &mut counter,
        &mut contig_stats,
        &options,
        &main_progress,
    )?;
    finish_progress_bars(&contig_stats, &main_progress);

    // Create the report model
    let report = report::collect_analysis_report(&bam_file, &bam_stats, &contig_stats, &counter)?;

    // Generate HTML report
    report::write_html_report(&report, &output_summary)?;

    Ok(())
}

fn setup_progress_bars() -> (MultiProgress, ProgressBar) {
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    main_progress.set_message("Analyzing coverage...");
    (multi_progress, main_progress)
}

fn initialize_contig_stats(
    header: &bam::HeaderView,
    options: &CallableOptions,
    multi_progress: &MultiProgress,
) -> Result<HashMap<usize, ContigProfiler>, Box<dyn Error>> {
    let mut contig_stats = HashMap::new();
    let selected_contigs: Option<std::collections::HashSet<String>> = options
        .selected_contigs
        .as_ref()
        .map(|contigs| contigs.iter().cloned().collect());

    for tid in 0..header.target_names().len() {
        let contig_name = std::str::from_utf8(header.tid2name(tid as u32))?;
        if let Some(ref selected) = selected_contigs {
            if !selected.contains(contig_name) {
                continue;
            }
        }
        let length = header.target_len(tid as u32).unwrap_or(0) as usize;
        contig_stats.insert(
            tid,
            ContigProfiler::new(
                contig_name.to_string(),
                length,
                multi_progress,
                options.clone(),
            ),
        );
    }
    Ok(contig_stats)
}

fn validate_contig_selection(
    contig_stats: &HashMap<usize, ContigProfiler>,
    options: &CallableOptions,
) -> Result<(), Box<dyn Error>> {
    if let Some(ref selected) = options.selected_contigs {
        if contig_stats.is_empty() {
            return Err(format!(
                "None of the specified contigs ({}) were found in the BAM file",
                selected
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )
            .into());
        }
    }
    Ok(())
}

fn initialize_counter(
    contig_stats: &HashMap<usize, ContigProfiler>,
    output_bed: &str,
) -> std::io::Result<CallableProfiler> {
    let largest_contig_length = contig_stats
        .values()
        .filter(|stats| stats.name != "chrM")
        .map(|stats| stats.length)
        .max()
        .unwrap_or(0) as u32;
    CallableProfiler::new(output_bed, largest_contig_length)
}

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

fn process_contigs(
    bam: &mut bam::IndexedReader,
    fasta: &mut IndexedReader<File>,
    header: &bam::HeaderView,
    counter: &mut CallableProfiler,
    contig_stats: &mut HashMap<usize, ContigProfiler>,
    options: &CallableOptions,
    main_progress: &ProgressBar,
) -> Result<(), Box<dyn Error>> {
    let mut contig_tids: Vec<_> = contig_stats.keys().cloned().collect();
    contig_tids.sort();

    for &tid in contig_tids.iter() {
        let contig_name = contig_stats[&tid].name.clone();
        main_progress.set_message(format!("Processing {}...", contig_name));

        process_single_contig(bam, fasta, header, counter, contig_stats, options, tid)?;
    }
    Ok(())
}

fn process_single_contig(
    bam: &mut bam::IndexedReader,
    fasta: &mut IndexedReader<File>,
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

            let mut seq = Vec::new();
            fasta.fetch(contig, current_pos as u64, (current_pos + 1) as u64)?;
            fasta.read(&mut seq)?;
            let ref_base = seq.first().map(|&b| b).unwrap_or(b'N');

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

        let mut seq = Vec::new();
        fasta.fetch(contig, pos as u64, (pos + 1) as u64)?;
        fasta.read(&mut seq)?;
        let ref_base = seq.first().map(|&b| b).unwrap_or(b'N');

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

        let mut seq = Vec::new();
        fasta.fetch(contig, current_pos as u64, (current_pos + 1) as u64)?;
        fasta.read(&mut seq)?;
        let ref_base = seq.first().map(|&b| b).unwrap_or(b'N');

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

fn finish_progress_bars(
    contig_stats: &HashMap<usize, ContigProfiler>,
    main_progress: &ProgressBar,
) {
    for (_, stats) in contig_stats.iter() {
        stats.progress_bar.finish_with_message("Complete");
    }
    main_progress.finish_with_message("Coverage analysis complete!");
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
