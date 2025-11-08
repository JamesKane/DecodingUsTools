use crate::types::ReferenceGenome;
use crate::utils::cache::TreeType;
use rust_htslib::bam::Read;

pub fn validate_reference<R: Read>(
    bam: &R,
    tree_type: TreeType,
) -> Result<(ReferenceGenome, String), Box<dyn std::error::Error>> {
    let header = bam.header();
    
    // Detect reference genome from BAM header
    let genome = ReferenceGenome::from_header(header)
        .ok_or("Could not determine reference genome from BAM header")?;

    // Define possible chromosome names based on build and type
    let possible_names = match (tree_type, genome.name()) {
        (TreeType::YDNA, "GRCh38") => vec!["chrY", "Y", "NC_000024.10", "CM000686.2"],
        (TreeType::YDNA, "GRCh37") => vec!["Y", "chrY"], // hg19 == GRCh37; allow both
        (TreeType::YDNA, "T2T-CHM13v2.0") => vec!["Y", "chrY", "CP086569.2", "NC_060948.1"],
        (TreeType::MTDNA, _) => vec!["chrM", "MT", "M", "NC_001807"],
        _ => vec![]
    };

    // Find the first matching sequence name
    let sequence_name = possible_names.iter()
        .find(|&name| header.tid(name.as_bytes()).is_some())
        .ok_or_else(|| format!(
            "No valid sequence found in BAM. Tried: {}",
            possible_names.join(", ")
        ))?;

    println!("Found sequence as: {}", sequence_name);

    Ok((genome, sequence_name.to_string()))
}
