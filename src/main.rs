mod cli;
mod commands;
mod utils;

use clap::Parser;

fn main() {
    let args = cli::Args::parse();

    let result = match args.command {
        cli::Commands::Coverage { bam_file, output_file } => {
            commands::coverage::run(bam_file, output_file)
        }
        cli::Commands::FindYBranch { bam_file, output_file, min_depth, min_quality } => {
            commands::y_branch::run(bam_file, output_file, min_depth, min_quality)
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
