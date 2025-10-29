use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct HaplogroupExport {
    pub analysis_type: HaplogroupType,
    pub top_result: HaplogroupResult,
    pub all_results: Vec<HaplogroupResult>,
    pub call_quality: CallQuality,
    pub snp_details: Option<SnpDetails>,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum HaplogroupType {
    #[serde(rename = "Y-DNA")]
    YDna,
    #[serde(rename = "mtDNA")]
    MtDna,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct HaplogroupResult {
    pub haplogroup: String,
    pub score: f64,
    pub confidence: String,  // "High", "Medium", "Low"
    pub matching_snps: usize,
    pub mismatching_snps: usize,
    pub total_snps_tested: usize,
    pub depth: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CallQuality {
    pub average_depth: f64,
    pub positions_tested: usize,
    pub positions_called: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SnpDetails {
    pub matching: Vec<SnpCall>,
    pub mismatching: Vec<SnpCall>,
    pub no_calls: Vec<SnpPosition>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SnpCall {
    pub name: String,
    pub position: u32,
    pub called_base: char,
    pub expected_base: char,
    pub depth: u32,
    pub frequency: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SnpPosition {
    pub name: String,
    pub position: u32,
}