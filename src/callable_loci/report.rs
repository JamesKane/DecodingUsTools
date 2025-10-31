use super::BamStats;
use crate::callable_loci::profilers::callable_profiler::CallableProfiler;
use crate::callable_loci::profilers::contig_profiler::ContigProfiler;
use crate::export::formats::coverage::{
    ContigExport, CoverageExport,
};

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Builds a structured CoverageExport from analysis results
pub fn build_coverage_export(
    contig_stats: &HashMap<usize, ContigProfiler>,
    counter: &CallableProfiler,
    bam_stats: &BamStats,
) -> Result<CoverageExport, Box<dyn Error>> {
    use crate::callable_loci::types::CalledState;
    use crate::export::formats::coverage::{
        ContigExport, CoverageExport, CoverageSummary, QualityMetrics, StateDistribution,
    };
    use std::sync::Arc;

    let mut total_bases = 0u64;
    let mut callable_bases = 0u64;
    let mut total_depth = 0f64;
    let mut total_mapq = 0f64;
    let mut total_baseq = 0f64;
    let mut q30_bases = 0u64;
    let mut total_quality_positions = 0u64;
    let mut total_unique_reads = 0u64;

    let mut contig_exports = Vec::new();

    let mut sorted_contigs: Vec<_> = contig_stats.values().collect();
    sorted_contigs.sort_by(|a, b| compare_contig_names(&a.name, &b.name));

    for stats in sorted_contigs {
        let counts = counter.get_contig_counts(&stats.name);
        let quality_stats = stats.get_quality_stats();

        let coverage_percent = if stats.length > 0 {
            (stats.n_covered_bases as f64 / stats.length as f64) * 100.0
        } else {
            0.0
        };

        let average_depth = if stats.n_covered_bases > 0 {
            stats.summed_coverage as f64 / stats.n_covered_bases as f64
        } else {
            0.0
        };

        total_bases += stats.length as u64;
        callable_bases += counts[CalledState::CALLABLE as usize];
        total_depth += average_depth * stats.length as f64;
        total_mapq += quality_stats.average_mapq * stats.length as f64;
        total_baseq += quality_stats.average_baseq * stats.length as f64;
        q30_bases += (quality_stats.q30_percentage / 100.0 * stats.length as f64) as u64;
        total_quality_positions += stats.length as u64;
        total_unique_reads += stats.n_reads as u64;

        let state_distribution = StateDistribution {
            ref_n: counts[CalledState::REF_N as usize],
            callable: counts[CalledState::CALLABLE as usize],
            no_coverage: counts[CalledState::NO_COVERAGE as usize],
            low_coverage: counts[CalledState::LOW_COVERAGE as usize],
            excessive_coverage: counts[CalledState::EXCESSIVE_COVERAGE as usize],
            poor_mapping_quality: counts[CalledState::POOR_MAPPING_QUALITY as usize],
        };

        contig_exports.push(ContigExport {
            name: stats.name.clone(),
            length: stats.length,
            unique_reads: stats.n_reads as u64,
            coverage_percent,
            average_depth,
            covered_bases: stats.n_covered_bases,
            total_bases: stats.length as u64,
            quality_stats,
            state_distribution,
            coverage_histogram: None,
        });
    }

    let average_depth = if total_bases > 0 {
        total_depth / total_bases as f64
    } else {
        0.0
    };

    let summary = CoverageSummary {
        aligner: bam_stats.aligner().to_string(),
        reference_build: bam_stats.reference_build().to_string(),
        sequencing_platform: bam_stats.infer_platform(),
        read_length: bam_stats.average_read_length(),
        total_bases,
        callable_bases,
        callable_percentage: if total_bases > 0 {
            (callable_bases as f64 / total_bases as f64) * 100.0
        } else {
            0.0
        },
        average_depth,
        contigs_analyzed: contig_stats.len(),
    };

    let quality_metrics = QualityMetrics {
        average_mapq: if total_quality_positions > 0 {
            total_mapq / total_quality_positions as f64
        } else {
            0.0
        },
        average_baseq: if total_quality_positions > 0 {
            total_baseq / total_quality_positions as f64
        } else {
            0.0
        },
        q30_percentage: if total_quality_positions > 0 {
            (q30_bases as f64 / total_quality_positions as f64) * 100.0
        } else {
            0.0
        },
    };

    Ok(CoverageExport {
        summary,
        contigs: Arc::new(contig_exports),
        quality_metrics,
        total_unique_reads
    })
}

pub fn write_html_report(
    export: &CoverageExport,
    bam_stats: &BamStats,
    output_path: &str,
) -> Result<(), Box<dyn Error>> {
    let mut html = String::new();

    // Add header template
    html.push_str(include_str!("templates/report_header.html"));

    // Add BAM Statistics section
    write_bam_stats_section(&mut html, export, bam_stats)?;

    // Add Contig Analysis section
    write_contig_analysis_section(&mut html, &export.contigs)?;

    // Add footer template and scripts
    html.push_str(include_str!("templates/report_footer.html"));

    let mut writer = BufWriter::new(File::create(output_path)?);
    writer.write_all(html.as_bytes())?;
    writer.flush()?;

    Ok(())
}

fn write_bam_stats_section(
    html: &mut String,
    export: &CoverageExport,
    bam_stats: &BamStats,
) -> Result<(), Box<dyn Error>> {
    html.push_str(r#"<section class='stats-box'>"#);
    html.push_str(&format!(
        r#"<h2>BAM Statistics <span class='sample-note'>(based on first {} reads)</span></h2>"#,
        bam_stats.max_samples
    ));

    html.push_str("<div class='stats-columns'>");
    html.push_str("<dl>");
    html.push_str(&format!(
        r#"<dt>Reference Build</dt><dd>{}</dd>
            <dt>Aligner</dt><dd>{}</dd>
            <dt>Sequencing Platform</dt><dd>{}</dd>
            <dt>Average read length</dt><dd>{} bp</dd>
            <dt>Total Unique Reads</dt><dd>{}</dd>
            <dt>Total Bases</dt><dd>{}</dd>
            "#,
        export.summary.reference_build,
        export.summary.aligner,
        export.summary.sequencing_platform,
        export.summary.read_length,
        export.total_unique_reads,
        export.summary.total_bases
    ));
    html.push_str("</dl>");

    html.push_str("<dl>");
    html.push_str(&format!(
        r#"<dt>Callable Bases</dt><dd>{}</dd>
            <dt>Callable Percentage</dt><dd>{:.2}%</dd>
            <dt>Average Depth</dt><dd>{:.2}×</dd>
            <dt>Contigs Analyzed</dt><dd>{}</dd>
            <dt>Average MapQ</dt><dd>{:.1}</dd>
            <dt>Average BaseQ</dt><dd>{:.1}</dd>"#,
        export.summary.callable_bases,
        export.summary.callable_percentage,
        export.summary.average_depth,
        export.summary.contigs_analyzed,
        export.quality_metrics.average_mapq,
        export.quality_metrics.average_baseq
    ));
    html.push_str("</dl>");
    html.push_str("</div>");
    html.push_str("</section>");

    Ok(())
}

fn write_contig_analysis_section(
    html: &mut String,
    contigs: &[ContigExport],
) -> Result<(), Box<dyn Error>> {
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
    contig: &ContigExport,
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

    add_table_row(html, "Length", format!("{} bp", contig.length));
    add_table_row(html, "Unique Reads", contig.unique_reads);
    add_table_row(html, "Covered Bases", contig.covered_bases);
    add_table_row(
        html,
        "Coverage Percent",
        format!("{:.2}%", contig.coverage_percent),
    );
    add_table_row(
        html,
        "Average Depth",
        format!("{:.2}×", contig.average_depth),
    );

    // Quality stats
    html.push_str(r#"<tr><td colspan="2" style="font-weight: bold; background-color: #f5f5f5;">Quality Metrics</td></tr>"#);
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
    add_table_row(
        html,
        "Q30 Percentage",
        format!("{:.2}%", contig.quality_stats.q30_percentage),
    );

    // State distribution
    html.push_str(r#"<tr><td colspan="2" style="font-weight: bold; background-color: #f5f5f5;">State Distribution</td></tr>"#);
    add_table_row(html, "Reference N", contig.state_distribution.ref_n);
    add_table_row(html, "Callable", contig.state_distribution.callable);
    add_table_row(html, "No Coverage", contig.state_distribution.no_coverage);
    add_table_row(html, "Low Coverage", contig.state_distribution.low_coverage);
    add_table_row(
        html,
        "Excessive Coverage",
        contig.state_distribution.excessive_coverage,
    );
    add_table_row(
        html,
        "Poor Mapping Quality",
        contig.state_distribution.poor_mapping_quality,
    );

    html.push_str("</tbody></table>");

    // Add coverage plot if available
    let plot_path = format!("{}_coverage.svg", contig.name);
    if Path::new(&plot_path).exists() {
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

/// Compare contig names for natural sorting (chr1, chr2, ..., chr10, chrX, chrY, chrM)
/// Supports both "chr" prefixed (b38/hs1) and unprefixed (b37) chromosome names
fn compare_contig_names(a: &str, b: &str) -> std::cmp::Ordering {
    use std::cmp::Ordering;

    // Extract prefix and suffix
    let (a_prefix, a_suffix) = split_contig_name(a);
    let (b_prefix, b_suffix) = split_contig_name(b);

    // First compare prefixes (e.g., "chr" vs "chr")
    match a_prefix.cmp(b_prefix) {
        Ordering::Equal => {}
        other => return other,
    }

    // Define chromosome order priority
    let order = |s: &str| -> (usize, u32) {
        // Try to parse as number
        if let Ok(num) = s.parse::<u32>() {
            return (0, num); // Numeric chromosomes come first
        }

        // Special chromosomes
        match s {
            "X" => (1, 0),
            "Y" => (2, 0),
            "M" | "MT" => (3, 0),
            _ => (4, 0), // Everything else comes last, alphabetically
        }
    };

    let (a_cat, a_num) = order(a_suffix);
    let (b_cat, b_num) = order(b_suffix);

    match a_cat.cmp(&b_cat) {
        Ordering::Equal => {
            if a_cat == 0 {
                // Both are numeric, compare numbers
                a_num.cmp(&b_num)
            } else {
                // Same category, use lexicographic comparison
                a_suffix.cmp(b_suffix)
            }
        }
        other => other,
    }
}

/// Split contig name into prefix and suffix (e.g., "chr10" -> ("chr", "10"))
fn split_contig_name(name: &str) -> (&str, &str) {
    // Find where the numeric/identifier part starts
    let split_pos = name
        .find(|c: char| c.is_ascii_digit() || c == 'X' || c == 'Y' || c == 'M')
        .unwrap_or(name.len());

    name.split_at(split_pos)
}
