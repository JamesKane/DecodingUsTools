use crate::haplogroup::types::{Haplogroup, HaplogroupTree};
use crate::utils::cache::TreeType;
use std::error::Error;

pub(crate) mod ftdna;

pub trait TreeProvider {
    fn url(&self, tree_type: TreeType) -> &str;
    fn cache_prefix(&self, tree_type: TreeType) -> &str;
    fn progress_message(&self, tree_type: TreeType) -> String;
    fn parse_tree(&self, data: &str) -> Result<HaplogroupTree, Box<dyn Error>>;
    fn build_tree(&self, tree: &HaplogroupTree, node_id: u32, tree_type: TreeType) -> Option<Haplogroup>;
}