use crate::callable_loci::types::{ContigCoverageStats, ContigQualityStats};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Debug, Serialize, Deserialize)]
pub struct CoverageExport {
    pub summary: CoverageSummary,
    #[serde(with = "arc_serde")]
    pub contigs: Arc<Vec<ContigExport>>,
    pub quality_metrics: QualityMetrics,
}

// Implementation of Clone using Arc for the non-cloneable Vec<ContigExport>
impl Clone for CoverageExport {
    fn clone(&self) -> Self {
        CoverageExport {
            summary: self.summary.clone(),
            contigs: Arc::clone(&self.contigs),
            quality_metrics: self.quality_metrics.clone(),
        }
    }
}

// Helper module for serializing Arc<Vec<T>>
mod arc_serde {
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use std::sync::Arc;

    pub fn serialize<S, T>(val: &Arc<Vec<T>>, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
        T: Serialize,
    {
        (**val).serialize(serializer)
    }

    pub fn deserialize<'de, D, T>(deserializer: D) -> Result<Arc<Vec<T>>, D::Error>
    where
        D: Deserializer<'de>,
        T: Deserialize<'de>,
    {
        Vec::deserialize(deserializer).map(Arc::new)
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CoverageSummary {
    pub total_bases: u64,
    pub callable_bases: u64,
    pub callable_percentage: f64,
    pub average_depth: f64,
    pub contigs_analyzed: usize,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ContigExport {
    pub name: String,
    pub length: usize,
    pub coverage_stats: ContigCoverageStats,
    pub quality_stats: ContigQualityStats,
    pub state_distribution: StateDistribution,

    // Optional: include plot data for GUI rendering
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coverage_histogram: Option<Vec<HistogramBin>>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct StateDistribution {
    pub ref_n: u64,
    pub callable: u64,
    pub no_coverage: u64,
    pub low_coverage: u64,
    pub excessive_coverage: u64,
    pub poor_mapping_quality: u64,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct HistogramBin {
    pub start: u32,
    pub end: u32,
    pub state: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct QualityMetrics {
    pub average_mapq: f64,
    pub average_baseq: f64,
    pub q30_percentage: f64,
}