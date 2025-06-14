use rust_htslib::bam::Read;
use crate::haplogroup::types::{ReferenceGenome};
use crate::utils::cache::TreeType;

pub fn validate_reference<R: Read>(
    bam: &R,
    tree_type: TreeType,
) -> Result<ReferenceGenome, Box<dyn std::error::Error>> {
    let header = bam.header();

    // Detect reference genome from BAM header
    let genome = ReferenceGenome::from_header(header)
        .ok_or("Could not determine reference genome from BAM header")?;

    // Get required chromosome based on tree type
    let required_chr = match tree_type {
        TreeType::YDNA => "chrY",
        TreeType::MTDNA => "chrM",
    };

    // Get accession info for the required chromosome
    let accession = genome.get_accession(required_chr)
        .ok_or_else(|| format!("No accession found for {} in {}", required_chr, genome.name()))?;

    // Check if either the common name or accession is present in the BAM header
    let has_required_seq = header.tid(accession.common_name.as_bytes()).is_some() ||
        header.tid(accession.accession.as_bytes()).is_some();

    if !has_required_seq {
        return Err(format!(
            "BAM file missing required sequence: {} ({})",
            accession.common_name,
            accession.accession
        ).into());
    }

    Ok(genome)
}
