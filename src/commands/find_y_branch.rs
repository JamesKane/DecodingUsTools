use crate::haplogroup;
use crate::utils::cache::TreeType;

pub fn run(
    bam_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    haplogroup::analyze_haplogroup(
        bam_file,
        output_file,
        min_depth,
        min_quality,
        TreeType::YDNA,
    )
}
