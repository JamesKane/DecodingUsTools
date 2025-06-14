use crate::cli;
use crate::haplogroup;
use crate::utils::cache::TreeType;

pub fn run(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    provider: cli::TreeProvider,
    show_snps: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    haplogroup::analyze_haplogroup(
        bam_file,
        reference_file,
        output_file,
        min_depth,
        min_quality,
        TreeType::YDNA,
        provider,
        show_snps,
    )
}
