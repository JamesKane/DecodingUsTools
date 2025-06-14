use super::TreeProvider;
use crate::haplogroup::types::{
    Haplogroup, HaplogroupNode, HaplogroupTree, LociCoordinate, LociType, Locus,
};
use crate::utils::cache::TreeType;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;

#[derive(Deserialize)]
pub struct FtdnaVariant {
    #[serde(default)]
    pub variant: String,
    #[serde(default)]
    pub position: Option<i32>,
    #[serde(default)]
    pub ancestral: String,
    #[serde(default)]
    pub derived: String,
    #[serde(default)]
    pub region: String,
    #[serde(default)]
    pub id: Option<u32>,
}

impl From<FtdnaVariant> for Locus {
    fn from(v: FtdnaVariant) -> Self {
        let mut coordinates = HashMap::new();
        if let Some(pos) = v.position {
            coordinates.insert(
                "GRCh38".to_string(),
                LociCoordinate {
                    position: pos.unsigned_abs(),
                    chromosome: "chrY".to_string(),
                    ancestral: v.ancestral,
                    derived: v.derived,
                }
            );
        }

        Locus {
            name: v.variant,
            loci_type: LociType::SNP,
            coordinates,
        }
    }
}

#[derive(Deserialize)]
struct FtdnaNode {
    #[serde(rename = "haplogroupId")]
    haplogroup_id: u32,
    #[serde(rename = "parentId", default)]
    parent_id: u32,
    name: String,
    #[serde(rename = "isRoot")]
    is_root: bool,
    root: String,
    #[serde(rename = "kitsCount")]
    kits_count: u32,
    #[serde(rename = "subBranches")]
    sub_branches: u32,
    #[serde(rename = "bigYCount")]
    big_y_count: u32,
    #[serde(default)]
    variants: Vec<FtdnaVariant>,
    #[serde(default)]
    children: Vec<u32>,
}

#[derive(Deserialize)]
struct FtdnaTreeJson {
    #[serde(rename = "allNodes")]
    nodes: HashMap<String, FtdnaNode>,
}

pub struct FtdnaTreeProvider;

impl FtdnaTreeProvider {
    pub fn new() -> Self {
        Self
    }
}

impl TreeProvider for FtdnaTreeProvider {
    fn url(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "https://www.familytreedna.com/public/y-dna-haplotree/get",
            TreeType::MTDNA => "https://www.familytreedna.com/public/mt-dna-haplotree/get",
        }
    }

    fn cache_prefix(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "ftdna-ytree",
            TreeType::MTDNA => "ftdna-mttree",
        }
    }

    fn progress_message(&self, tree_type: TreeType) -> String {
        format!(
            "Downloading FTDNA {} tree...",
            match tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        )
    }

    fn parse_tree(&self, data: &str) -> Result<HaplogroupTree, Box<dyn Error>> {
        let ftdna_tree: FtdnaTreeJson = serde_json::from_str(data)?;
        let all_nodes = ftdna_tree
            .nodes
            .into_iter()
            .map(|(k, n)| {
                let loci = n.variants.into_iter().map(Locus::from).collect();
                (
                    k,
                    HaplogroupNode {
                        haplogroup_id: n.haplogroup_id,
                        parent_id: n.parent_id,
                        name: n.name,
                        is_root: n.is_root,
                        loci,
                        children: n.children,
                    },
                )
            })
            .collect();

        Ok(HaplogroupTree { all_nodes })
    }

    fn build_tree(
        &self,
        tree: &HaplogroupTree,
        node_id: u32,
        tree_type: TreeType,
    ) -> Option<Haplogroup> {
        let node_str = node_id.to_string();
        let node = tree.all_nodes.get(&node_str)?;

        let children = node
            .children
            .iter()
            .filter_map(|&child_id| self.build_tree(tree, child_id, tree_type))
            .collect();

        Some(Haplogroup {
            name: node.name.clone(),
            parent: if node.parent_id == 0 {
                None
            } else {
                Some(
                    tree.all_nodes
                        .get(&node.parent_id.to_string())?
                        .name
                        .clone(),
                )
            },
            loci: node.loci.clone(),
            children,
        })
    }

    fn supported_builds(&self) -> Vec<String> {
        vec!["GRCh38".to_string(), "rCRS".to_string()]
    }
}