use crate::callable_loci;
use crate::export::formats::coverage::CoverageExport;
use crate::api::{ApiResult, ProgressCallback, ProgressEvent};
use crate::callable_loci::profilers::{
    bam_stats::BamStats,
    callable_profiler::CallableProfiler,
    contig_profiler::ContigProfiler,
};
use rust_htslib::{bam, faidx};
use std::collections::HashMap;
use std::path::PathBuf;
use rust_htslib::bam::Read;
use crate::utils::bam_reader::BamReaderFactory;
use crate::utils::progress_manager::ProgressManager;

pub struct CoverageAnalyzer {
    progress_callback: Option<ProgressCallback>,
}

impl CoverageAnalyzer {
    pub fn new() -> Self {
        Self { progress_callback: None }
    }

    pub fn with_progress(mut self, callback: ProgressCallback) -> Self {
        self.progress_callback = Some(callback);
        self
    }

    pub async fn analyze(
        &self,
        input: CoverageInput,
    ) -> ApiResult<CoverageOutput> {
        self.emit_progress(ProgressEvent::Started {
            task: "Coverage Analysis".to_string(),
        });

        let result = self.run_analysis(input).await?;

        self.emit_progress(ProgressEvent::Completed {
            task: "Coverage Analysis".to_string(),
        });

        Ok(result)
    }

    async fn run_analysis(&self, input: CoverageInput) -> ApiResult<CoverageOutput> {

        // Collect BAM statistics
        let mut bam_stats = BamStats::new(10000);
        let progress_mgr = ProgressManager::new();
        bam_stats.collect_stats(&input.bam_file, Some(&input.reference_file), &progress_mgr)
            .map_err(|e| crate::api::ApiError::from(format!("Failed to collect BAM stats: {}", e)))?;

        // Prepare output paths
        let output_bed = input.output_bed.clone()
            .unwrap_or_else(|| "callable_regions.bed".to_string());
        let output_summary = input.output_summary.clone()
            .unwrap_or_else(|| "summary.html".to_string());

        let options = input.options.clone().with_contigs(input.contigs.clone());

        let mut bam = BamReaderFactory::open_indexed(&input.bam_file, Some(&input.reference_file))
            .map_err(|e| crate::api::ApiError::from(format!("Failed to open BAM file: {}", e)))?;
        let header = bam.header().clone();

        let mut fasta = faidx::Reader::from_path(&input.reference_file)
            .map_err(|e| crate::api::ApiError::from(format!("Failed to open reference: {}", e)))?;

        // Initialize contig stats (without progress bars for API usage)
        let mut contig_stats = initialize_contig_stats(&header, &options)
            .map_err(|e| crate::api::ApiError::from(e))?;
        validate_contig_selection(&contig_stats, &options)
            .map_err(|e| crate::api::ApiError::from(e))?;

        let mut counter = initialize_counter(&contig_stats, &output_bed)
            .map_err(|e| crate::api::ApiError::from(e))?;

        // Process all contigs
        process_contigs_api(
            &mut bam,
            &mut fasta,
            &header,
            &mut counter,
            &mut contig_stats,
            &options,
        ).map_err(|e| crate::api::ApiError::from(e))?;

        // Build the export structure
        let export = build_coverage_export(&contig_stats, &counter, &bam_stats)
            .map_err(|e| crate::api::ApiError::from(e))?;

        // Collect coverage plot files
        let coverage_plots = collect_coverage_plots(&contig_stats)
            .map_err(|e| crate::api::ApiError::from(e))?;

        // Generate HTML report
        callable_loci::report::write_html_report(&export, &bam_stats, &output_summary)
            .map_err(|e| crate::api::ApiError::from(e.to_string()))?;

        Ok(CoverageOutput {
            export,
            files: OutputFiles {
                bed_file: output_bed,
                summary_html: output_summary,
                coverage_plots,
            },
        })
    }

    fn emit_progress(&self, event: ProgressEvent) {
        if let Some(callback) = &self.progress_callback {
            callback(event);
        }
    }
}

#[derive(Debug, Clone)]
pub struct CoverageInput {
    pub bam_file: String,
    pub reference_file: String,
    pub contigs: Option<Vec<String>>,
    pub options: callable_loci::CallableOptions,
    pub output_bed: Option<String>,
    pub output_summary: Option<String>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct CoverageOutput {
    pub export: CoverageExport,
    pub files: OutputFiles,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct OutputFiles {
    pub bed_file: String,
    pub summary_html: String,
    pub coverage_plots: Vec<String>,
}

// Helper functions adapted from callable_loci::mod.rs

fn initialize_contig_stats(
    header: &bam::HeaderView,
    options: &callable_loci::CallableOptions,
) -> Result<HashMap<usize, ContigProfiler>, String> {
    let mut contig_stats = HashMap::new();
    let selected_contigs: Option<std::collections::HashSet<String>> = options
        .selected_contigs
        .as_ref()
        .map(|contigs| contigs.iter().cloned().collect());

    // Create a ProgressManager even for API usage (progress bars won't be visible in non-terminal contexts)
    let progress_mgr = ProgressManager::new();

    for tid in 0..header.target_names().len() {
        let contig_name = std::str::from_utf8(header.tid2name(tid as u32))
            .map_err(|e| format!("Invalid contig name: {}", e))?;

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
                &progress_mgr,
                options.clone(),
            ),
        );
    }
    Ok(contig_stats)
}

fn validate_contig_selection(
    contig_stats: &HashMap<usize, ContigProfiler>,
    options: &callable_loci::CallableOptions,
) -> Result<(), String> {
    if let Some(ref selected) = options.selected_contigs {
        if contig_stats.is_empty() {
            return Err(format!(
                "None of the specified contigs ({}) were found in the BAM file",
                selected
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
    }
    Ok(())
}

fn initialize_counter(
    contig_stats: &HashMap<usize, ContigProfiler>,
    output_bed: &str,
) -> Result<CallableProfiler, String> {
    let largest_contig_length = contig_stats
        .values()
        .filter(|stats| stats.name != "chrM")
        .map(|stats| stats.length)
        .max()
        .unwrap_or(0) as u32;

    CallableProfiler::new(output_bed, largest_contig_length)
        .map_err(|e| format!("Failed to create CallableProfiler: {}", e))
}

fn process_contigs_api(
    bam: &mut bam::IndexedReader,
    fasta: &mut faidx::Reader,
    header: &bam::HeaderView,
    counter: &mut CallableProfiler,
    contig_stats: &mut HashMap<usize, ContigProfiler>,
    options: &callable_loci::CallableOptions,
) -> Result<(), String> {
    let mut contig_tids: Vec<_> = contig_stats.keys().cloned().collect();
    contig_tids.sort();

    for &tid in contig_tids.iter() {
        process_single_contig_api(bam, fasta, header, counter, contig_stats, options, tid)?;
    }
    Ok(())
}

fn process_single_contig_api(
    bam: &mut bam::IndexedReader,
    fasta: &mut faidx::Reader,
    header: &bam::HeaderView,
    counter: &mut CallableProfiler,
    contig_stats: &mut HashMap<usize, ContigProfiler>,
    options: &callable_loci::CallableOptions,
    tid: usize,
) -> Result<(), String> {
    // This is adapted from callable_loci::mod.rs::process_single_contig
    // Refer to that implementation for the full logic
    callable_loci::process_single_contig(
        bam, fasta, header, counter, contig_stats, options, tid
    ).map_err(|e| format!("Error processing contig: {}", e))
}

fn build_coverage_export(
    contig_stats: &HashMap<usize, ContigProfiler>,
    counter: &CallableProfiler,
    bam_stats: &BamStats,
) -> Result<CoverageExport, String> {
    callable_loci::report::build_coverage_export(contig_stats, counter, bam_stats)
        .map_err(|e| format!("Failed to build coverage export: {}", e))
}

fn collect_coverage_plots(
    contig_stats: &HashMap<usize, ContigProfiler>,
) -> Result<Vec<String>, String> {
    let mut plots = Vec::new();
    for stats in contig_stats.values() {
        let plot_path = format!("{}_coverage.svg", stats.name);
        if PathBuf::from(&plot_path).exists() {
            plots.push(plot_path);
        }
    }
    Ok(plots)
}