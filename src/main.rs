mod bam_fixer;
mod callable_loci;
mod cli;
mod commands;
mod config;
mod haplogroup;
mod utils;
mod vendor;
mod vg;
mod generated;
pub mod sequence_processor;
mod export;
mod api;
mod types;

use clap::Parser;

use crate::bam_fixer::BamFixer;
use anyhow::Context;
use std::error::Error;
use crate::api::coverage::{CoverageAnalyzer, CoverageInput};
use crate::utils::cache::TreeType;

fn main() -> Result<(), Box<dyn Error>> {
    let args = cli::Args::parse();

    match args.command {
        cli::Commands::Coverage {
            bam_file,
            reference_file,
            output_file,
            summary_file,
            contigs,
            min_depth,
            max_depth,
            min_mapping_quality,
            min_base_quality,
            min_depth_for_low_mapq,
            max_low_mapq,
            max_low_mapq_fraction,
        } => {
            let options = callable_loci::CallableOptions::new(
                min_depth,
                max_depth,
                min_mapping_quality,
                min_base_quality,
                min_depth_for_low_mapq,
                max_low_mapq,
                max_low_mapq_fraction,
            );
            let input = CoverageInput {
                bam_file: bam_file.clone(),
                reference_file: reference_file.clone(),
                contigs,
                options,
                output_bed: Some(output_file),
                output_summary: Some(summary_file),
            };

            let analyzer = CoverageAnalyzer::new();
            let rt = tokio::runtime::Builder::new_current_thread()
                .enable_all()
                .build()
                .unwrap();
            let result = rt.block_on(analyzer.analyze(input))?;
        }
        cli::Commands::FindYBranch {
            bam_file,
            reference_file,
            output_file,
            min_depth,
            min_quality,
            provider,
            show_snps,
        } => {
            commands::find_branch::run(
                bam_file,
                reference_file,
                output_file,
                min_depth,
                min_quality,
                provider,
                show_snps,
                TreeType::YDNA
            )?;
        }
        cli::Commands::FindMtBranch {
            bam_file,
            reference_file,
            output_file,
            min_depth,
            min_quality,
            provider,
            show_snps,
        } => {
            commands::find_branch::run(
                bam_file,
                reference_file,
                output_file,
                min_depth,
                min_quality,
                provider,
                show_snps,
                TreeType::MTDNA
            )?;
        }
        cli::Commands::FixSurjectedBam {
            input_bam,
            reference_file,
            output_bam,
            keep_temp,
        } => {
            let fixer = BamFixer::new(reference_file, input_bam, output_bam, keep_temp)
                .context("Failed to initialize BAM fixer")?;
            fixer.run()?;
        }
        cli::Commands::Fingerprint {
            input_file,
            reference_file,
            ksize,
            scaled,
            max_frequency,
            output_file,
            region,
        } => {
            commands::fingerprint::run(
                input_file,
                reference_file,
                ksize,
                scaled,
                max_frequency,
                output_file,
                region,
            )?;
        }
    }

    Ok(())
}
