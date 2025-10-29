#[derive(Clone, Debug)]
pub struct CallableOptions {
    pub min_depth: u32,
    pub max_depth: u32,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub min_depth_for_low_mapq: u32,
    pub max_low_mapq: u8,
    pub max_low_mapq_fraction: f64,
    pub selected_contigs: Option<Vec<String>>,
}

impl CallableOptions {
    pub fn new(
        min_depth: u32,
        max_depth: u32,
        min_mapping_quality: u8,
        min_base_quality: u8,
        min_depth_for_low_mapq: u32,
        max_low_mapq: u8,
        max_low_mapq_fraction: f64,
    ) -> Self {
        Self {
            min_depth,
            max_depth,
            min_mapping_quality,
            min_base_quality,
            min_depth_for_low_mapq,
            max_low_mapq,
            max_low_mapq_fraction,
            selected_contigs: None,
        }
    }

    pub fn with_contigs(mut self, contigs: Option<Vec<String>>) -> Self {
        self.selected_contigs = contigs;
        self
    }
}