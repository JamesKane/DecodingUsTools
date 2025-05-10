use serde::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize, Debug, Clone)]
pub struct Snp {
    pub(crate) position: u32,
    pub(crate) ancestral: String,
    pub(crate) derived: String,
    #[serde(skip)]
    pub(crate) chromosome: String,
    #[serde(skip)]
    pub(crate) build: String,
}

#[derive(Deserialize, Debug)]
pub struct Variant {
    pub variant: String,
    #[serde(rename = "position", default)]
    pub pos: u32,
    #[serde(default)]
    pub ancestral: String,
    #[serde(default)]
    pub derived: String,
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupNode {
    pub haplogroup_id: u32,
    pub parent_id: u32,
    pub name: String,
    pub is_root: bool,
    pub variants: Vec<Variant>,
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
    pub(crate) snps: Vec<Snp>,
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