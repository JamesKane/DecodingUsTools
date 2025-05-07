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

    /// Find closest Y-DNA branch for a sample (Requires hg38 aligned BAM.)
    FindYBranch {
        /// Input BAM file
        bam_file: String,
        /// Output file for haplogroup results
        output_file: String,
        /// Minimum read depth for SNP calling (default: 10)
        #[arg(long, default_value = "10")]
        min_depth: u32,
        /// Minimum mapping quality (default: 20)
        #[arg(long, default_value = "20")]
        min_quality: u8,
    },

    /// Find closest MT-DNA branch for a sample (Requires hg38 aligned BAM.)
    FindMtBranch {
        /// Input BAM file
        bam_file: String,
        /// Output file for haplogroup results
        output_file: String,
        /// Minimum read depth for SNP calling (default: 10)
        #[arg(long, default_value = "10")]
        min_depth: u32,
        /// Minimum mapping quality (default: 20)
        #[arg(long, default_value = "20")]
        min_quality: u8,
    },

}