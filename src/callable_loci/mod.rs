mod options;
mod profilers;
mod types;
mod utils;

use crate::callable_loci::types::CalledState;
use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
pub use options::CallableOptions;
use profilers::{callable_profiler::CallableProfiler, contig_profiler::ContigProfiler};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufWriter;

pub fn run(
    bam_file: String,
    reference_file: String,
    output_bed: String,
    output_summary: String,
    mut options: CallableOptions,
    contigs: Option<Vec<String>>,
) -> Result<(), Box<dyn Error>> {
    println!("Requested contigs: {:?}", contigs);
    options = options.with_contigs(contigs);

    let mut bam = bam::IndexedReader::from_path(&bam_file)?;
    let header = bam.header().clone();
    let mut fasta = IndexedReader::from_file(&reference_file)?;

    // Setup progress tracking
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    main_progress.set_message("Analyzing coverage...");

    // Initialize contig tracking
    let mut contig_stats = HashMap::new();
    let selected_contigs: Option<std::collections::HashSet<String>> = options
        .selected_contigs
        .as_ref()
        .map(|contigs| contigs.iter().cloned().collect());

    // Initialize stats for each contig we'll process
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
                &multi_progress,
                options.clone(),
            ),
        );
    }

    // Validate contig selection
    if let Some(ref selected) = selected_contigs {
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

    // Initialize counter
    let largest_contig_length = contig_stats
        .values()
        .filter(|stats| stats.name != "chrM")
        .map(|stats| stats.length)
        .max()
        .unwrap_or(0) as u32;
    let mut counter = CallableProfiler::new(&output_bed, largest_contig_length)?;

    // Process each selected contig
    let mut contig_tids: Vec<_> = contig_stats.keys().cloned().collect();
    contig_tids.sort();
    for &tid in contig_tids.iter() {
        let contig_name = contig_stats[&tid].name.clone();
        main_progress.set_message(format!("Processing {}...", contig_name));

        // Seek to contig and set up pileup
        bam.fetch((tid as u32, 0, header.target_len(tid as u32).unwrap_or(0)))?;
        let mut pileup = bam.pileup();
        pileup.set_max_depth(if options.max_depth > 0 {
            options.max_depth
        } else {
            500
        });

        // Process pileup for this contig
        for p in pileup {
            let pileup = p?;
            if pileup.tid() as usize != tid {
                break; // Move to next contig
            }

            let pos = pileup.pos();
            let contig = std::str::from_utf8(header.tid2name(pileup.tid()))?;

            // Update progress for current contig
            if let Some(stats) = contig_stats.get_mut(&tid) {
                stats.progress_bar.set_position(pos as u64);
            }

            // Get reference base
            let mut seq = Vec::new();
            fasta.fetch(contig, pos as u64, (pos + 1) as u64)?;
            fasta.read(&mut seq)?;
            let ref_base = seq.first().map(|&b| b).unwrap_or(b'N');

            let mut raw_depth = 0;
            let mut qc_depth = 0;
            let mut low_mapq_count = 0;

            // Process alignments
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

            counter.process_position(
                contig,
                pos,
                ref_base,
                raw_depth,
                qc_depth,
                low_mapq_count,
                &options,
            )?;

            // Update contig stats
            if let Some(stats) = contig_stats.get_mut(&tid) {
                stats.process_position(raw_depth, pileup.alignments());
            }
        }

        // Finish this contig
        let stats = &contig_stats[&tid];
        counter.finish_contig(&stats.name, stats.length as u32)?;
    }

    // Finish all progress bars
    for (_, stats) in contig_stats.iter() {
        stats.progress_bar.finish_with_message("Complete");
    }
    main_progress.finish_with_message("Coverage analysis complete!");

    // Write summary
    use std::io::Write;
    let mut summary_writer = BufWriter::new(File::create(output_summary)?);

    // Write header
    writeln!(
        summary_writer,
        "Contig|Start|Stop|UniqueReads|RefN|NoCoverage|LowCoverage|ExcessiveCoverage|PoorMappingQuality|Callable|CoveragePercent|AvgDepth|AvgMapQ|AvgBaseQ"
    )?;

    // Write stats for each contig
    let mut sorted_stats: Vec<_> = contig_stats.iter().collect();
    sorted_stats.sort_by(|a, b| utils::natural_sort::natural_cmp(&a.1.name, &b.1.name));
    for (_, stats) in sorted_stats {
        let counts = counter.get_contig_counts(&stats.name);
        writeln!(
            summary_writer,
            "{}|1|{}|{}|{}|{}|{}|{}|{}|{}|{:.2}|{:.2}|{:.1}|{:.1}",
            stats.name,
            stats.length,
            stats.n_selected_reads,
            counts[CalledState::REF_N as usize],
            counts[CalledState::NO_COVERAGE as usize],
            counts[CalledState::LOW_COVERAGE as usize],
            counts[CalledState::EXCESSIVE_COVERAGE as usize],
            counts[CalledState::POOR_MAPPING_QUALITY as usize],
            counts[CalledState::CALLABLE as usize],
            (stats.n_covered_bases as f64 / stats.length as f64) * 100.0,
            stats.summed_coverage as f64 / stats.length as f64,
            if stats.n_selected_reads > 0 {
                stats.summed_mapq as f64 / stats.n_selected_reads as f64
            } else {
                0.0
            },
            if stats.quality_bases > 0 {
                stats.summed_baseq as f64 / stats.quality_bases as f64
            } else {
                0.0
            }
        )?;
    }
    summary_writer.flush()?;

    Ok(())
}
