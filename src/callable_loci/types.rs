use serde::{Deserialize, Serialize};

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledState {
    REF_N,
    CALLABLE,
    NO_COVERAGE,
    LOW_COVERAGE,
    EXCESSIVE_COVERAGE,
    POOR_MAPPING_QUALITY,
}

#[derive(Debug)]
pub(crate) struct CoverageRange {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) state: CalledState,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BamAnalysisReport {
    pub metadata: BamMetadata,
    pub bam_stats: BamStatistics,
    pub contig_analyses: Vec<ContigAnalysis>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BamMetadata {
    pub reference_build: String,
    pub aligner: String,
    pub sample_count: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BamStatistics {
    pub average_read_length: f64,
    pub paired_percentage: f64,
    pub average_insert_size: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ContigAnalysis {
    pub name: String,
    pub length: usize,
    pub coverage_stats: ContigCoverageStats,
    pub quality_stats: ContigQualityStats,
    pub coverage_plot: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ContigCoverageStats {
    pub unique_reads: u64,
    pub state_counts: ContigStateCounts,
    pub coverage_percent: f64,
    pub average_depth: f64,
    pub covered_bases: u64,
    pub total_bases: u64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ContigQualityStats {
    pub average_mapq: f64,
    pub average_baseq: f64,
    pub q30_percentage: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ContigStateCounts {
    pub reference_n: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
    pub callable: u64,
}
