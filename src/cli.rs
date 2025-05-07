use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Analyze BAM file coverage and callability
    Coverage {
        /// Path to the BAM file
        bam_file: String,

        /// Output file for the coverage report
        #[arg(short = 'o', long = "output", default_value = "cov_report.txt")]
        output_file: String,
    },
}