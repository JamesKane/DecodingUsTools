use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct FingerprintExport {
    pub hexdigest: String,
    pub parameters: FingerprintParameters,
    pub statistics: FingerprintStatistics,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub hashes: Option<Vec<u64>>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FingerprintParameters {
    pub ksize: usize,
    pub scaled: usize,
    pub max_frequency: Option<u32>,
    pub region: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FingerprintStatistics {
    pub sequences_processed: u64,
    pub unique_kmers: usize,
}