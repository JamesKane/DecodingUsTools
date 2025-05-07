mod cli;
mod commands;
mod utils;
mod haplogroup;

use clap::Parser;

fn main() {
    let args = cli::Args::parse();

    let result = match args.command {
        cli::Commands::Coverage {
            bam_file,
            output_file,
        } => commands::coverage::run(bam_file, output_file),
        cli::Commands::FindYBranch {
            bam_file,
            output_file,
            min_depth,
            min_quality,
        } => commands::find_y_branch::run(bam_file, output_file, min_depth, min_quality),
        cli::Commands::FindMtBranch {
            bam_file,
            output_file,
            min_depth,
            min_quality,
        } => commands::find_mt_branch::run(bam_file, output_file, min_depth, min_quality),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
