use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize, Debug, Clone)]
pub enum LociType {
    SNP,
    INDEL,
}

#[derive(Deserialize, Debug, Clone)]
pub struct LociCoordinate {
    pub position: u32,
    pub chromosome: String,
    pub ancestral: String,
    pub derived: String,
}

#[derive(Deserialize, Debug, Clone)]
pub struct Locus {
    pub name: String,
    pub loci_type: LociType,
    pub coordinates: HashMap<String, LociCoordinate>, // genbank_id -> coordinate info
}

#[derive(Deserialize)]
pub struct ApiCoordinate {
    pub start: u32,
    pub stop: u32,
    pub anc: String,
    pub der: String,
}

#[derive(Deserialize)]
pub struct ApiVariant {
    pub name: String,
    pub coordinates: HashMap<String, ApiCoordinate>,
    #[serde(rename = "variantType")]
    pub variant_type: String,
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupNode {
    pub haplogroup_id: u32,
    pub parent_id: u32,
    pub name: String,
    pub is_root: bool,
    pub loci: Vec<Locus>,
    pub children: Vec<u32>,
}

#[derive(Deserialize)]
pub struct HaplogroupTree {
    #[serde(rename = "allNodes")]
    pub all_nodes: HashMap<String, HaplogroupNode>,
}

#[derive(Debug)]
pub struct Haplogroup {
    pub(crate) name: String,
    pub(crate) parent: Option<String>,
    pub(crate) loci: Vec<Locus>,
    pub(crate) children: Vec<Haplogroup>,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct HaplogroupScore {
    pub(crate) matches: usize,
    pub(crate) ancestral_matches: usize,
    pub(crate) no_calls: usize,
    pub(crate) total_snps: usize,
    pub(crate) score: f64,
    depth: usize, // Added depth field
}

impl Default for HaplogroupScore {
    fn default() -> Self {
        Self {
            matches: 0,
            ancestral_matches: 0,
            no_calls: 0,
            total_snps: 0,
            score: 0.0,
            depth: 0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct HaplogroupResult {
    pub(crate) name: String,
    pub(crate) score: f64,
    pub(crate) matching_snps: u32,
    pub(crate) mismatching_snps: u32,
    pub(crate) ancestral_matches: u32,
    pub(crate) no_calls: u32,
    pub(crate) total_snps: u32,
    pub(crate) cumulative_snps: u32, // Total unique SNPs from root to this branch
    pub(crate) depth: u32,
}

impl TryFrom<ApiVariant> for Locus {
    type Error = Box<dyn std::error::Error>;

    fn try_from(api_variant: ApiVariant) -> Result<Self, Self::Error> {
        let loci_type = match api_variant.variant_type.as_str() {
            "SNP" => LociType::SNP,
            "INDEL" => LociType::INDEL,
            unknown => return Err(format!("Unknown variant type: {}", unknown).into()),
        };

        let coordinates = api_variant
            .coordinates
            .into_iter()
            .map(|(genbank_id, coord)| {
                (
                    genbank_id,
                    LociCoordinate {
                        position: coord.start,
                        chromosome: "chrY".to_string(),
                        ancestral: coord.anc,
                        derived: coord.der,
                    },
                )
            })
            .collect();

        Ok(Locus {
            name: api_variant.name,
            loci_type,
            coordinates,
        })
    }
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
        for _tid in 0..header.target_count() {
            // Get SQ lines from header text
            let header_text = String::from_utf8_lossy(header.as_bytes());

            // Look for assembly information in header text
            if header_text.contains("AS:GRCh38") || header_text.contains("GCA_000001405.15") {
                return Some(ReferenceGenome::GRCh38);
            } else if header_text.contains("AS:CHM13") || header_text.contains("GCA_009914755.4") {
                return Some(ReferenceGenome::CHM13v2);
            } else if header_text.contains("AS:GRCh37") || header_text.contains("GCA_000001405.1") {
                return Some(ReferenceGenome::GRCh37);
            }
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
