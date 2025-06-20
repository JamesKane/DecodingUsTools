#[derive(Debug, Clone)]
pub struct Sequence {
    pub data: Vec<u8>,
    pub id: Option<String>,
    pub quality: Option<Vec<u8>>,
    pub metadata: SequenceMetadata,
}

#[derive(Debug, Clone, Default)]
pub struct SequenceMetadata {
    pub is_mapped: bool,
    pub chromosome: Option<String>,
    pub position: Option<u64>,
    pub mapping_quality: Option<u8>,
}