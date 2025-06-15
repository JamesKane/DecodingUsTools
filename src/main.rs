mod bam_fixer;
mod cli;
mod commands;
mod haplogroup;
mod utils;
mod vendor;
mod config;

use clap::Parser;

use crate::bam_fixer::BamFixer;
use std::error::Error;
use anyhow::Context;

fn main() -> Result<(), Box<dyn Error>> {
    let args = cli::Args::parse();

    match args.command {
        cli::Commands::Coverage {
            bam_file,
            reference_file,
            output_file,
            min_depth,
            max_depth,
            min_mapping_quality,
            min_base_quality,
            min_depth_for_low_mapq,
            max_low_mapq,
            max_low_mapq_fraction,
        } => {
            let options = commands::coverage::CallableOptions::new(
                min_depth,
                max_depth,
                min_mapping_quality,
                min_base_quality,
                min_depth_for_low_mapq,
                max_low_mapq,
                max_low_mapq_fraction,
            );
            let summary_file = format!("{}.summary", output_file);
            commands::coverage::run(bam_file, reference_file, output_file, summary_file, options)?;
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
            commands::find_y_branch::run(
                bam_file,
                reference_file,
                output_file,
                min_depth,
                min_quality,
                provider,
                show_snps,
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
            commands::find_mt_branch::run(
                bam_file,
                reference_file,
                output_file,
                min_depth,
                min_quality,
                provider,
                show_snps,
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
    }

    Ok(())
}
