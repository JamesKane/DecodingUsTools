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

        /// Path to the reference FASTA file
        #[arg(short = 'r', long = "reference")]
        reference_file: String,

        /// Output file for the coverage report
        #[arg(short = 'o', long = "output", default_value = "cov_report.txt")]
        output_file: String,

        /// Minimum depth for callable regions
        #[arg(long, default_value = "4")]
        min_depth: u32,

        /// Maximum depth before considering excessive coverage
        #[arg(long, default_value_t = u32::MAX)]
        max_depth: u32,

        /// Minimum mapping quality for callable regions
        #[arg(long, default_value = "10")]
        min_mapping_quality: u8,

        /// Minimum base quality for callable regions
        #[arg(long, default_value = "20")]
        min_base_quality: u8,

        /// Minimum depth before considering poor mapping quality
        #[arg(long, default_value = "10")]
        min_depth_for_low_mapq: u32,

        /// Maximum MAPQ to be considered low quality
        #[arg(long, default_value = "1")]
        max_low_mapq: u8,

        /// Maximum fraction of low MAPQ reads allowed
        #[arg(long, default_value = "0.1")]
        max_low_mapq_fraction: f64,
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
