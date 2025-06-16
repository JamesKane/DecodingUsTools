mod options;
mod profilers;
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
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

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

    write_summary(&contig_stats, &counter, &output_summary, &bam_stats)?;
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
                0,  // raw_depth
                0,  // qc_depth
                0,  // low_mapq_count
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
            0,  // raw_depth
            0,  // qc_depth
            0,  // low_mapq_count
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

use std::path::Path;
fn write_summary(
    contig_stats: &HashMap<usize, ContigProfiler>,
    counter: &CallableProfiler,
    output_summary: &str,
    bam_stats: &BamStats,
) -> Result<(), Box<dyn Error>> {
    let mut html = String::new();

    html.push_str(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BAM Analysis Report</title>
    <style>
        :root {
            --primary-color: #007bff;
            --border-color: #ddd;
            --bg-light: #f5f5f5;
        }
        body { 
            font-family: system-ui, -apple-system, sans-serif;
            margin: 0;
            padding: 2rem;
        }
        main {
            max-width: 1200px;
            margin: 0 auto;
        }
        .stats-box { 
            background: var(--bg-light);
            padding: 1.25rem;
            border-radius: 0.5rem;
            margin-bottom: 2rem;
        }
        .sample-note {
            color: #666;
            font-style: italic;
            font-size: 0.9em;
        }
        .tabs {
            margin-top: 1.25rem;
        }
        .tab-list {
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
            margin-bottom: 1.25rem;
            list-style: none;
            padding: 0;
        }
        .tab-button {
            padding: 0.625rem 1.25rem;
            border: none;
            background: var(--bg-light);
            cursor: pointer;
            border-radius: 0.25rem;
            font-size: 1rem;
        }
        .tab-button:hover {
            background: #e0e0e0;
        }
        .tab-button[aria-selected="true"] {
            background: var(--primary-color);
            color: white;
        }
        .tab-panel {
            display: none;
            padding: 1.25rem;
            border: 1px solid var(--border-color);
            border-radius: 0.5rem;
        }
        .tab-panel[aria-hidden="false"] {
            display: block;
        }
        table {
            border-collapse: collapse;
            width: 100%;
            margin-bottom: 1.25rem;
        }
        th, td {
            border: 1px solid var(--border-color);
            padding: 0.5rem;
            text-align: left;
        }
        th { background-color: var(--bg-light); }
        .coverage-plot {
            width: 100%;
            margin-top: 1.25rem;
        }
        .coverage-plot img {
            width: 100%;
            height: auto;
            display: block;
            max-width: 100%;
            object-fit: contain;
        }
        figure {
            margin: 0;
        }
        figcaption {
            text-align: center;
            margin-top: 0.5rem;
            color: #666;
        }
                .contig-selector {
            margin: 1em 0;
        }
        .contig-selector select {
            width: 300px;
            padding: 8px;
            border: 1px solid var(--border-color);
            border-radius: 0.5rem;
            font-size: 1rem;
            background: var(--bg-light);
        }
        .contig-selector select:focus {
            outline: none;
            border-color: var(--primary-color);
            box-shadow: 0 0 0 2px rgba(0, 123, 255, 0.2);
        }
    </style>
</head>
<body>
<main>
    <h1>BAM Analysis Report</h1>
"#,
    );

    // Add BAM Statistics section
    html.push_str(r#"<section class='stats-box'>"#);
    html.push_str(&format!(
        r#"<h2>BAM Statistics <span class='sample-note'>(based on first {} reads)</span></h2>"#,
        bam_stats.max_samples
    ));

    html.push_str("<dl>");
    let stats = bam_stats.get_stats();
    html.push_str(&format!(
        r#"<dt>Average read length</dt><dd>{:.1} bp</dd>
        <dt>Paired reads</dt><dd>{:.1}%</dd>
        <dt>Average insert size</dt><dd>{:.1} bp</dd>"#,
        stats.get("average_read_length").unwrap_or(&0.0),
        stats.get("paired_percentage").unwrap_or(&0.0),
        stats.get("average_insert_size").unwrap_or(&0.0)
    ));
    html.push_str("</dl></section>");
    
    let mut sorted_stats: Vec<_> = contig_stats.iter().collect();
    sorted_stats.sort_by(|a, b| utils::natural_sort::natural_cmp(&a.1.name, &b.1.name));

    html.push_str(r#"<section class='tabs' role='tablist'>"#);

    // Add the contig selector
    html.push_str(r#"<div class="contig-selector">
    <select id="contig-select" onchange="switchToContig(this.value)" aria-label="Select contig">"#);

    // Add options for each contig
    for (i, (_, stats)) in sorted_stats.iter().enumerate() {
        html.push_str(&format!(
            r#"<option value="{}" {}>{}</option>"#,
            i,
            if i == 0 { "selected" } else { "" },
            stats.name
        ));
    }

    html.push_str("</select></div>");

    // Add tab panels
    for (i, (_, stats)) in sorted_stats.iter().enumerate() {
        html.push_str(&format!(
            r#"<div class="tab-panel" role="tabpanel" id="panel-{}" aria-labelledby="tab-{}" aria-hidden="{}">"#,
            i,
            i,
            i != 0
        ));

        // Add contig statistics table
        html.push_str("<table>");
        html.push_str("<thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>");

        let counts = counter.get_contig_counts(&stats.name);
        add_table_row(&mut html, "Unique Reads", stats.n_selected_reads);
        add_table_row(&mut html, "Reference N", counts[CalledState::REF_N as usize]);
        add_table_row(&mut html, "No Coverage", counts[CalledState::NO_COVERAGE as usize]);
        add_table_row(&mut html, "Low Coverage", counts[CalledState::LOW_COVERAGE as usize]);
        add_table_row(&mut html, "Excessive Coverage", counts[CalledState::EXCESSIVE_COVERAGE as usize]);
        add_table_row(&mut html, "Poor Mapping Quality", counts[CalledState::POOR_MAPPING_QUALITY as usize]);
        add_table_row(&mut html, "Callable", counts[CalledState::CALLABLE as usize]);

        let coverage_percent = (stats.n_covered_bases as f64 / stats.length as f64) * 100.0;
        let avg_depth = stats.summed_coverage as f64 / stats.length as f64;
        let avg_mapq = if stats.n_selected_reads > 0 {
            stats.summed_mapq as f64 / stats.n_selected_reads as f64
        } else {
            0.0
        };
        let avg_baseq = if stats.quality_bases > 0 {
            stats.summed_baseq as f64 / stats.quality_bases as f64
        } else {
            0.0
        };

        add_table_row(&mut html, "Coverage Percent", format!("{:.2}%", coverage_percent));
        add_table_row(&mut html, "Average Depth", format!("{:.2}×", avg_depth));
        add_table_row(&mut html, "Average MapQ", format!("{:.1}", avg_mapq));
        add_table_row(&mut html, "Average BaseQ", format!("{:.1}", avg_baseq));

        html.push_str("</tbody></table>");


        // Add coverage plot if it exists
        let plot_path = format!("{}_coverage.svg", stats.name);
        if Path::new(&plot_path).exists() {
            html.push_str(&format!(
                r#"<figure class='coverage-plot'>
                    <img src="{}" alt="Coverage distribution for {}" loading="lazy">
                    <figcaption>Coverage distribution for {}</figcaption>
                </figure>"#,
                plot_path, stats.name, stats.name
            ));
        }

        html.push_str("</div>");
    }

    html.push_str(
        r#"</section></main>
<script>
function switchToContig(index) {
    // Hide all panels
    document.querySelectorAll('[role="tabpanel"]').forEach(panel => {
        panel.setAttribute('aria-hidden', 'true');
    });
    
    // Show selected panel
    const selectedPanel = document.getElementById('panel-' + index);
    selectedPanel.setAttribute('aria-hidden', 'false');
}

// Initial setup
document.addEventListener('DOMContentLoaded', () => {
    const select = document.getElementById('contig-select');
    if (select) {
        switchToContig(select.value);
    }
});
</script>
</body>
</html>"#,
    );

    let mut writer = BufWriter::new(File::create(output_summary)?);
    writer.write_all(html.as_bytes())?;
    writer.flush()?;

    Ok(())
}

fn add_table_row(html: &mut String, label: &str, value: impl std::fmt::Display) {
    html.push_str(&format!("<tr><td>{}</td><td>{}</td></tr>", label, value));
}
