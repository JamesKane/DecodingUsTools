mod cli;
mod commands;
mod haplogroup;
mod utils;

use clap::Parser;

use std::error::Error;

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
            commands::coverage::run(bam_file, reference_file, output_file, options)?;
        }
        cli::Commands::FindYBranch {
            bam_file,
            reference_file,
            output_file,
            min_depth,
            min_quality,
        } => {
            commands::find_y_branch::run(bam_file, reference_file, output_file, min_depth, min_quality)?;
        }
        cli::Commands::FindMtBranch {
            bam_file,
            reference_file,
            output_file,
            min_depth,
            min_quality,
        } => {
            commands::find_mt_branch::run(bam_file, reference_file, output_file, min_depth, min_quality)?;
        }

    }

    Ok(())
}
