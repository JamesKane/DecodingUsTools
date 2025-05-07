mod cli;
mod commands;

use clap::Parser;

fn main() {
    let args = cli::Args::parse();

    let result = match args.command {
        cli::Commands::Coverage { bam_file, output_file } => {
            commands::coverage::run(bam_file, output_file)
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
