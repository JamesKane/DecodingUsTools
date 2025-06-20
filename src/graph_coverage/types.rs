use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeCoverage {
    pub coverage: u32,
    pub forward_coverage: u32,
    pub reverse_coverage: u32,
    pub total_quality: f64,
    pub unique_reads: HashSet<String>,
}

impl Default for NodeCoverage {
    fn default() -> Self {
        Self {
            coverage: 0,
            forward_coverage: 0,
            reverse_coverage: 0,
            total_quality: 0.0,
            unique_reads: HashSet::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeCoverage {
    pub from_node: String,
    pub to_node: String,
    pub coverage: u32,
    pub avg_quality: f64,
}

#[derive(Default, Debug, Serialize, Deserialize)]
#[derive(Clone)]
pub struct PathCoverage {
    pub path_name: String,
    pub positions: Vec<u32>,
    pub avg_coverage: f64,
    pub min_coverage: u32,
    pub max_coverage: u32,
}

#[derive(Default, Debug, Serialize, Deserialize)]
pub struct GraphStats {
    pub total_nodes: u64,
    pub covered_nodes: u64,
    pub avg_node_coverage: f64,
    pub path_coverages: HashMap<String, PathCoverage>,
    pub low_coverage_nodes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphCoverageStats {
    pub nodes: HashMap<String, NodeCoverage>,
    pub edges: Vec<EdgeCoverage>,
    pub total_reads: u64,
    pub mapped_reads: u64,
    pub unmapped_reads: u64,
    pub low_quality_reads: u64,
}

impl Default for GraphCoverageStats {
    fn default() -> Self {
        Self {
            nodes: HashMap::new(),
            edges: Vec::new(),
            total_reads: 0,
            mapped_reads: 0,
            unmapped_reads: 0,
            low_quality_reads: 0,
        }
    }
}
