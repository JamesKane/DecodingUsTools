use super::TreeProvider;
use crate::haplogroup::types::{
    Haplogroup, HaplogroupNode, HaplogroupTree, LociCoordinate, LociType, Locus,
};
use crate::utils::cache::TreeType;
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;

#[derive(Deserialize)]
pub struct ApiCoordinate {
    pub start: u32,
    pub stop: u32,
    pub anc: String,
    pub der: String,
}

#[derive(Deserialize)]
pub struct ApiVariant {
    pub name: String,
    pub coordinates: HashMap<String, ApiCoordinate>,
    #[serde(rename = "variantType")]
    pub variant_type: String,
}

#[derive(Deserialize)]
pub struct ApiNode {
    pub name: String,
    #[serde(rename = "parentName")]
    pub parent_name: String,
    pub variants: Vec<ApiVariant>,
    #[serde(rename = "lastUpdated")]
    pub last_updated: String,
    #[serde(rename = "isBackbone")]
    pub is_backbone: bool,
}

pub struct DecodingUsTreeProvider;

impl DecodingUsTreeProvider {
    pub fn new() -> Self {
        Self
    }
}

impl TreeProvider for DecodingUsTreeProvider {
    fn url(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "https://decoding-us.com/api/v1/y-tree",
            TreeType::MTDNA => panic!("MT-DNA tree not yet supported by DecodingUs"),
        }
    }

    fn cache_prefix(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "decodingus-ytree",
            TreeType::MTDNA => panic!("MT-DNA tree not yet supported by DecodingUs"),
        }
    }

    fn progress_message(&self, tree_type: TreeType) -> String {
        format!(
            "Downloading DecodingUs {} tree...",
            match tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => panic!("MT-DNA tree not yet supported by DecodingUs"),
            }
        )
    }

    fn parse_tree(&self, data: &str) -> Result<HaplogroupTree, Box<dyn Error>> {
        let api_nodes: Vec<ApiNode> = serde_json::from_str(data)?;

        // First, create a mapping of names to IDs
        let mut name_to_id: HashMap<String, u32> = HashMap::new();
        for (idx, node) in api_nodes.iter().enumerate() {
            name_to_id.insert(node.name.clone(), idx as u32);
        }

        // Pre-calculate children relationships
        let mut children_map: HashMap<String, Vec<u32>> = HashMap::new();
        for (idx, node) in api_nodes.iter().enumerate() {
            if !node.parent_name.is_empty() {
                children_map
                    .entry(node.parent_name.clone())
                    .or_default()
                    .push(idx as u32);
            }
        }

        // Create the nodes map
        let mut all_nodes = HashMap::new();
        for (idx, node) in api_nodes.into_iter().enumerate() {
            let haplogroup_id = idx as u32;
            let parent_id = if node.parent_name.is_empty() {
                0 // Root node
            } else {
                *name_to_id.get(&node.parent_name).unwrap_or(&0)
            };

            // Convert ApiVariant to Locus
            let loci = node.variants.into_iter()
                .map(|v| {
                    let loci_type = match v.variant_type.as_str() {
                        "SNP" => LociType::SNP,
                        _ => LociType::INDEL,
                    };

                    let coordinates = v.coordinates.into_iter()
                        .map(|(build, coord)| {
                            (build, LociCoordinate {
                                position: coord.start,
                                chromosome: "chrY".to_string(),
                                ancestral: coord.anc,
                                derived: coord.der,
                            })
                        })
                        .collect();

                    Locus {
                        name: v.name,
                        loci_type,
                        coordinates,
                    }
                })
                .collect();

            // Get children from pre-calculated map
            let children = children_map.get(&node.name).cloned().unwrap_or_default();

            all_nodes.insert(
                haplogroup_id.to_string(),
                HaplogroupNode {
                    haplogroup_id,
                    parent_id,
                    name: node.name,
                    is_root: parent_id == 0,
                    loci,
                    children,
                },
            );
        }

        Ok(HaplogroupTree { all_nodes })
    }

    fn build_tree(&self, tree: &HaplogroupTree, node_id: u32, tree_type: TreeType) -> Option<Haplogroup> {
        let node_str = node_id.to_string();
        let node = tree.all_nodes.get(&node_str)?;

        let children = node.children
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
        vec![
            "GRCh38".to_string(),
            "GRCh37".to_string(),
            "T2T-CHM13v2.0".to_string(),
        ]
    }
}