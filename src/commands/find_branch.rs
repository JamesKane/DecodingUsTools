use crate::cli::TreeProvider;
use crate::haplogroup;
use crate::utils::cache::TreeType;
use std::error::Error;

pub fn run(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    provider: TreeProvider,
    show_snps: bool,
    tree_type: TreeType,
) -> Result<(), Box<dyn Error>> {
    haplogroup::analyze_haplogroup(
        bam_file,
        reference_file,
        output_file,
        min_depth,
        min_quality,
        tree_type,
        provider,
        show_snps,
    )
}