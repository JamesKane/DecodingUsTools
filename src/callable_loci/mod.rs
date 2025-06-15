mod options;
mod profilers;
mod utils;

use crate::callable_loci::profilers::callable_profiler::CalledState;
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
    options: CallableOptions,
) -> Result<(), Box<dyn Error>> {
    let mut bam = bam::IndexedReader::from_path(&bam_file)?;
    let header = bam.header().clone();
    let mut fasta = IndexedReader::from_file(&reference_file)?;

    // Create progress bars
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    main_progress.set_message("Analyzing coverage...");

    // Create contig progress bars map using usize as key
    let mut contig_stats = HashMap::new();

    // First pass to get contig lengths and set up progress bars
    for tid in 0..header.target_names().len().min(25) {
        let contig_name = std::str::from_utf8(header.tid2name(tid as u32))?;
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

    // Find largest non-chrM contig length
    let largest_contig_length = contig_stats
        .values()
        .filter(|stats| stats.name != "chrM")
        .map(|stats| stats.length)
        .max()
        .unwrap_or(0) as u32;

    let mut counter = CallableProfiler::new(&output_bed, largest_contig_length)?;

    let mut pileup = bam.pileup();
    pileup.set_max_depth(if options.max_depth > 0 {
        options.max_depth
    } else {
        500
    });

    let mut current_contig = String::new();
    // Process pileup
    for p in pileup {
        let pileup = p?;
        let tid = pileup.tid() as usize; // Convert tid to usize here
        let pos = pileup.pos();
        let contig = std::str::from_utf8(header.tid2name(pileup.tid()))?;

        // Detect contig transition
        if current_contig != contig && !current_contig.is_empty() {
            // Use the length from contig_stats for the current contig
            let length = contig_stats
                .values()
                .find(|stats| stats.name == current_contig)
                .map(|stats| stats.length)
                .unwrap_or(0) as u32;
            counter.finish_contig(&current_contig, length)?;
        }
        current_contig = contig.to_string();

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
        // First determine the state based on the parameters
        let state = if ref_base == b'N' || ref_base == b'n' {
            CalledState::REF_N
        } else if raw_depth == 0 {
            CalledState::NO_COVERAGE
        } else if raw_depth >= options.min_depth_for_low_mapq
            && (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction
        {
            CalledState::POOR_MAPPING_QUALITY
        } else if qc_depth < options.min_depth {
            CalledState::LOW_COVERAGE
        } else if options.max_depth > 0 && qc_depth > options.max_depth {
            CalledState::EXCESSIVE_COVERAGE
        } else {
            CalledState::CALLABLE
        };

        let is_low_mapq = raw_depth >= options.min_depth_for_low_mapq
            && (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction;

        // Now call process_position with the correct parameters
        counter.process_position(contig, pos, qc_depth, is_low_mapq, state)?;

        // Update contig stats
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.process_position(raw_depth, pileup.alignments());
        }
    }

    // Finish all progress bars
    for (_, stats) in contig_stats.iter() {
        stats.progress_bar.finish_with_message("Complete");
    }
    main_progress.finish_with_message("Coverage analysis complete!");

    // Finish final contig
    let final_length = contig_stats
        .values()
        .find(|stats| stats.name == current_contig)
        .map(|stats| stats.length)
        .unwrap_or(0) as u32;
    counter.finish_contig(&current_contig, final_length)?;

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
