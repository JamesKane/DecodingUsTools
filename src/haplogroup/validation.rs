use rust_htslib::bam::Read;
use crate::haplogroup::types::Snp;

pub fn validate_hg38_reference<R: Read>(bam: &R) -> Result<(), Box<dyn std::error::Error>> {
    let header = bam.header();
    let sq_count = header.target_count();

    if sq_count == 0 {
        return Err("BAM file has no reference sequences".into());
    }

    // Check if chrY or chrM exists based on need
    let has_chry = header.tid(b"chrY").is_some();
    let has_chrm = header.tid(b"chrM").is_some();

    if !has_chry && !has_chrm {
        return Err("BAM file missing required chromosome (chrY or chrM)".into());
    }

    Ok(())
}

pub fn is_valid_snp(snp: &Snp) -> bool {
    match snp.chromosome.as_str() {
        "chrM" => true,
        "chrY" => snp.build == "hg38",
        _ => false,
    }
}