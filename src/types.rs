#[derive(clap::ValueEnum, Clone, Debug)]
pub enum Region {
    #[value(name = "full")]
    Full,
    #[value(name = "chrY")]
    Ychr,
    #[value(name = "chrM")]
    Mtchr,
}

impl Region {
    pub fn to_chromosome_names(&self) -> Vec<String> {
        match self {
            Region::Full => vec![],  // Empty vec means process all chromosomes
            Region::Ychr => vec!["chrY".to_string(), "Y".to_string()],  // Try both names
            Region::Mtchr => vec!["chrM".to_string(), "MT".to_string(), "M".to_string()],
        }
    }

    pub fn to_output_name(&self) -> &'static str {
        match self {
            Region::Full => "full",
            Region::Ychr => "chrY",  // Changed from "Ychr" to match BAM conventions
            Region::Mtchr => "chrM",
        }
    }
}

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum TreeProvider {
    #[value(name = "ftdna")]
    FTDNA,
    #[value(name = "decodingus")]
    DecodingUs,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ReferenceGenome {
    GRCh38,
    CHM13v2,
    GRCh37,
}

#[derive(Debug, Clone)]
pub struct GenomeAccession {
    pub genome: ReferenceGenome,
    pub chromosome: String,
    pub accession: String,
    pub common_name: String,
}

impl ReferenceGenome {
    pub fn accessions(&self) -> Vec<GenomeAccession> {
        match self {
            ReferenceGenome::GRCh38 => vec![
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrY".to_string(),
                    accession: "CM000686.2".to_string(), // GRCh38 Y chromosome
                    common_name: "chrY".to_string(),
                },
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrM".to_string(),
                    accession: "J01415.2".to_string(), // rCRS mitochondrial
                    common_name: "chrM".to_string(),
                },
            ],
            ReferenceGenome::CHM13v2 => vec![
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrY".to_string(),
                    accession: "CP086569.2".to_string(), // CHM13 Y chromosome
                    common_name: "chrY".to_string(),
                },
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrM".to_string(),
                    accession: "J01415.2".to_string(), // rCRS mitochondrial
                    common_name: "chrM".to_string(),
                },
            ],
            ReferenceGenome::GRCh37 => vec![
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrY".to_string(),
                    accession: "CM000686.1".to_string(), // GRCh37 Y chromosome
                    common_name: "chrY".to_string(),
                },
                GenomeAccession {
                    genome: self.clone(),
                    chromosome: "chrM".to_string(),
                    accession: "J01415.2".to_string(), // rCRS mitochondrial
                    common_name: "chrM".to_string(),
                },
            ],
        }
    }

    pub fn from_header(header: &rust_htslib::bam::HeaderView) -> Option<Self> {
        let header_text = String::from_utf8_lossy(header.as_bytes());

        // Look for assembly information in header text
        if header_text.contains("AS:GRCh38") || header_text.contains("GCA_000001405.15") {
            return Some(ReferenceGenome::GRCh38);
        } else if header_text.contains("AS:GRCh37") || header_text.contains("GCA_000001405.1") {
            return Some(ReferenceGenome::GRCh37);
        }

        // Check for CHM13v2 using multiple strategies
        if header_text.contains("AS:CHM13") || header_text.contains("GCA_009914755.4") {
            return Some(ReferenceGenome::CHM13v2);
        }

        // Check for CHM13v2 by characteristic patterns in UR field
        if header_text.contains("chm13") || header_text.contains("CHM13") || header_text.contains("t2t") || header_text.contains("T2T") {
            return Some(ReferenceGenome::CHM13v2);
        }

        // Check for CHM13v2 by characteristic chr1 length (248387328 bp)
        // and MD5 checksum
        if header_text.contains("SN:chr1") && header_text.contains("LN:248387328") {
            if header_text.contains("M5:e469247288ceb332aee524caec92bb22") {
                return Some(ReferenceGenome::CHM13v2);
            }
        }

        // Fallback: Check for characteristic GRCh38 chr1 length (248956422 bp)
        if header_text.contains("SN:chr1") && header_text.contains("LN:248956422") {
            return Some(ReferenceGenome::GRCh38);
        }

        // Fallback: Check for characteristic GRCh37 chr1 length (249250621 bp)
        if header_text.contains("SN:1") && header_text.contains("LN:249250621") {
            return Some(ReferenceGenome::GRCh37);
        }

        None
    }

    pub fn name(&self) -> &'static str {
        match self {
            ReferenceGenome::GRCh38 => "GRCh38",
            ReferenceGenome::CHM13v2 => "T2T-CHM13v2.0",
            ReferenceGenome::GRCh37 => "GRCh37",
        }
    }

    pub fn get_accession(&self, common_name: &str) -> Option<GenomeAccession> {
        self.accessions()
            .into_iter()
            .find(|acc| acc.common_name == common_name)
            .or_else(|| {
                self.accessions()
                    .into_iter()
                    .find(|acc| acc.accession == common_name)
            })
    }
}