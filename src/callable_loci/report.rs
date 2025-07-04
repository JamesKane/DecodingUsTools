use super::types::{
    BamAnalysisReport, BamMetadata, BamStatistics, ContigAnalysis, ContigCoverageStats,
    ContigQualityStats, ContigStateCounts,
};
use super::{detect_aligner, utils, BamStats, CallableProfiler, CalledState, ContigProfiler};
use crate::haplogroup::types::ReferenceGenome;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn collect_analysis_report(
    bam_file: &str,
    bam_stats: &BamStats,
    contig_stats: &HashMap<usize, ContigProfiler>,
    counter: &CallableProfiler,
) -> Result<BamAnalysisReport, Box<dyn Error>> {
    let bam = bam::IndexedReader::from_path(bam_file)?;
    let header = bam.header();

    // Collect metadata
    let metadata = BamMetadata {
        reference_build: ReferenceGenome::from_header(header)
            .map(|g| g.name().to_string())
            .unwrap_or_else(|| "Unknown".to_string()),
        aligner: detect_aligner(header),
        sample_count: bam_stats.max_samples,
    };

    // Collect BAM statistics
    let stats = BamStatistics {
        average_read_length: *bam_stats
            .get_stats()
            .get("average_read_length")
            .unwrap_or(&0.0),
        paired_percentage: *bam_stats
            .get_stats()
            .get("paired_percentage")
            .unwrap_or(&0.0),
        average_insert_size: *bam_stats
            .get_stats()
            .get("average_insert_size")
            .unwrap_or(&0.0),
    };

    // Collect contig analyses
    let mut contig_analyses = Vec::new();
    for (_, stats) in contig_stats.iter() {
        let counts = counter.get_contig_counts(&stats.name);
        let coverage_percent = (stats.n_covered_bases as f64 / stats.length as f64) * 100.0;
        let average_depth = stats.summed_coverage as f64 / stats.length as f64;

        let analysis = ContigAnalysis {
            name: stats.name.clone(),
            length: stats.length,
            coverage_stats: ContigCoverageStats {
                unique_reads: stats.n_selected_reads as u64,
                state_counts: ContigStateCounts {
                    reference_n: counts[CalledState::REF_N as usize],
                    no_coverage: counts[CalledState::NO_COVERAGE as usize],
                    low_coverage: counts[CalledState::LOW_COVERAGE as usize],
                    excessive_coverage: counts[CalledState::EXCESSIVE_COVERAGE as usize],
                    poor_mapping_quality: counts[CalledState::POOR_MAPPING_QUALITY as usize],
                    callable: counts[CalledState::CALLABLE as usize],
                },
                coverage_percent,
                average_depth,
            },
            quality_stats: ContigQualityStats {
                average_mapq: if stats.n_selected_reads > 0 {
                    stats.summed_mapq as f64 / stats.n_selected_reads as f64
                } else {
                    0.0
                },
                average_baseq: if stats.quality_bases > 0 {
                    stats.summed_baseq as f64 / stats.quality_bases as f64
                } else {
                    0.0
                },
            },
            coverage_plot: {
                let plot_path = format!("{}_coverage.svg", stats.name);
                if Path::new(&plot_path).exists() {
                    Some(plot_path)
                } else {
                    None
                }
            },
        };
        contig_analyses.push(analysis);
    }

    // Sort contig analyses by name
    contig_analyses.sort_by(|a, b| utils::natural_sort::natural_cmp(&a.name, &b.name));

    Ok(BamAnalysisReport {
        metadata,
        bam_stats: stats,
        contig_analyses,
    })
}

pub fn write_html_report(
    report: &BamAnalysisReport,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    let mut html = String::new();

    // Add header template
    html.push_str(include_str!("templates/report_header.html"));

    // Add BAM Statistics section
    write_bam_stats_section(&mut html, report)?;

    // Add Contig Analysis section
    write_contig_analysis_section(&mut html, &report.contig_analyses)?;

    // Add footer template and scripts
    html.push_str(include_str!("templates/report_footer.html"));

    let mut writer = BufWriter::new(File::create(output_path)?);
    writer.write_all(html.as_bytes())?;
    writer.flush()?;

    Ok(())
}

fn write_bam_stats_section(
    html: &mut String,
    report: &BamAnalysisReport,
) -> Result<(), Box<dyn Error>> {
    html.push_str(r#"<section class='stats-box'>"#);
    html.push_str(&format!(
        r#"<h2>BAM Statistics <span class='sample-note'>(based on first {} reads)</span></h2>"#,
        report.metadata.sample_count
    ));

    html.push_str("<dl>");
    html.push_str(&format!(
        r#"<dt>Reference Build</dt><dd>{}</dd>
        <dt>Aligner</dt><dd>{}</dd>
        <dt>Average read length</dt><dd>{:.1} bp</dd>
        <dt>Paired reads</dt><dd>{:.1}%</dd>
        <dt>Average insert size</dt><dd>{:.1} bp</dd>"#,
        report.metadata.reference_build,
        report.metadata.aligner,
        report.bam_stats.average_read_length,
        report.bam_stats.paired_percentage,
        report.bam_stats.average_insert_size
    ));
    html.push_str("</dl></section>");

    Ok(())
}

fn write_contig_analysis_section(html: &mut String, contigs: &[ContigAnalysis]) -> Result<(), Box<dyn Error>> {
    html.push_str(r#"<div class="contig-analysis">"#);

    // Add the contig selector
    html.push_str(
        r#"<div class="contig-selector">
        <select id="contig-select" onchange="switchToContig(this.value)" aria-label="Select contig">"#,
    );

    // Add options for each contig
    for (i, contig) in contigs.iter().enumerate() {
        html.push_str(&format!(
            r#"<option value="panel-{}" {}>{}</option>"#,
            i,
            if i == 0 { "selected" } else { "" },
            contig.name
        ));
    }

    html.push_str("</select></div>");

    // Container for panels
    html.push_str(r#"<div class="contig-panels">"#);

    // Add tab panels for each contig
    for (i, contig) in contigs.iter().enumerate() {
        write_contig_panel(html, contig, i)?;
    }

    html.push_str("</div></div>");

    Ok(())
}

fn write_contig_panel(
    html: &mut String,
    contig: &ContigAnalysis,
    index: usize,
) -> Result<(), Box<dyn Error>> {
    html.push_str(&format!(
        r#"<div class="tab-panel {}" id="panel-{}">"#,
        if index == 0 { "active" } else { "" },
        index,
    ));

    // Add contig statistics table
    html.push_str("<table>");
    html.push_str("<thead><tr><th>Metric</th><th>Value</th></tr></thead><tbody>");

    add_table_row(html, "Unique Reads", contig.coverage_stats.unique_reads);
    add_table_row(
        html,
        "Reference N",
        contig.coverage_stats.state_counts.reference_n,
    );
    add_table_row(
        html,
        "No Coverage",
        contig.coverage_stats.state_counts.no_coverage,
    );
    add_table_row(
        html,
        "Low Coverage",
        contig.coverage_stats.state_counts.low_coverage,
    );
    add_table_row(
        html,
        "Excessive Coverage",
        contig.coverage_stats.state_counts.excessive_coverage,
    );
    add_table_row(
        html,
        "Poor Mapping Quality",
        contig.coverage_stats.state_counts.poor_mapping_quality,
    );
    add_table_row(
        html,
        "Callable",
        contig.coverage_stats.state_counts.callable,
    );
    add_table_row(
        html,
        "Coverage Percent",
        format!("{:.2}%", contig.coverage_stats.coverage_percent),
    );
    add_table_row(
        html,
        "Average Depth",
        format!("{:.2}×", contig.coverage_stats.average_depth),
    );
    add_table_row(
        html,
        "Average MapQ",
        format!("{:.1}", contig.quality_stats.average_mapq),
    );
    add_table_row(
        html,
        "Average BaseQ",
        format!("{:.1}", contig.quality_stats.average_baseq),
    );

    html.push_str("</tbody></table>");

    // Add coverage plot if available
    if let Some(plot_path) = &contig.coverage_plot {
        html.push_str(&format!(
            r#"<figure class='coverage-plot'>
                <img src="{}" alt="Coverage distribution for {}" loading="lazy">
                <figcaption>Coverage distribution for {}</figcaption>
            </figure>"#,
            plot_path, contig.name, contig.name
        ));
    }

    html.push_str("</div>");
    Ok(())
}

fn add_table_row(html: &mut String, label: &str, value: impl std::fmt::Display) {
    html.push_str(&format!(r#"<tr><td>{}</td><td>{}</td></tr>"#, label, value));
}
