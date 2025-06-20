use crate::graph_coverage::coverage::GraphCoverageProcessor;
use crate::sequence_processor::core::SequenceReader;
use crate::sequence_processor::readers::GamReader;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

pub fn run(
    input_file: String,
    min_mapping_quality: u8,
    min_coverage: u32,
    output_file: Option<String>,
) -> Result<()> {
    let input_path = PathBuf::from(&input_file);

    if input_path.extension().and_then(|ext| ext.to_str()) != Some("gam") {
        anyhow::bail!("Input file must be a .gam file");
    }

    let progress = ProgressBarBuilder::new("Processing alignments")
        .with_template("{spinner:.green} [{elapsed_precise}] {pos} alignments processed ({per_sec})")
        .build()?;

    let mut processor = GraphCoverageProcessor::new(min_mapping_quality, min_coverage)
        .with_progress(progress.clone());

    let num_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);

    let mut reader = GamReader::new(&input_path)?;
    let stats = reader.read_sequences_with_threads(&mut processor, &progress, num_threads)?;

    // Print summary statistics
    let coverage_stats = processor.get_stats();
    println!("\nGraph Coverage Summary:");
    println!("Total reads processed: {}", coverage_stats.total_reads);
    println!("Mapped reads: {}", coverage_stats.mapped_reads);
    println!("Unmapped reads: {}", coverage_stats.unmapped_reads);
    println!("Low quality reads: {}", coverage_stats.low_quality_reads);
    println!("Number of nodes covered: {}", coverage_stats.nodes.len());
    println!("Number of edges covered: {}", coverage_stats.edges.len());

    // Output detailed results if requested
    if let Some(output_path) = output_file {
        let mut file = File::create(Path::new(&output_path))?;
        let json = serde_json::to_string_pretty(&coverage_stats)?;
        file.write_all(json.as_bytes())?;
        println!("\nDetailed coverage statistics written to: {}", output_path);
    }

    println!("\nProcessing completed:");
    println!("Processed sequences: {}", stats.processed);
    println!("Errors encountered: {}", stats.errors);
    println!("Sequences too short: {}", stats.too_short);

    Ok(())
}

