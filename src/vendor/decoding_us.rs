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
    parent_name: Option<String>,
    pub variants: Vec<ApiVariant>,
    #[serde(rename = "lastUpdated")]
    pub last_updated: String,
    #[serde(rename = "isBackbone")]
    pub is_backbone: bool,
}

fn deserialize_null_string<'de, D>(deserializer: D) -> Result<String, D::Error>
where
    D: serde::Deserializer<'de>,
{
    Ok(Option::deserialize(deserializer)?.unwrap_or_default())
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
        let mut root_id = None;

        for (idx, node) in api_nodes.iter().enumerate() {
            name_to_id.insert(node.name.clone(), idx as u32);
            if node.parent_name.is_none() || node.parent_name.as_ref().map_or(true, |p| p.is_empty()) {
                if root_id.is_some() {
                    return Err("Multiple root nodes found in tree".into());
                }
                root_id = Some(idx as u32);
            }
        }

        let root_id = root_id.ok_or("No root node found")?;

        // Create the nodes map
        let mut all_nodes = HashMap::new();
        for (idx, node) in api_nodes.into_iter().enumerate() {
            let haplogroup_id = idx as u32;
            let is_root = haplogroup_id == root_id;

            // Set parent_id to root_id for nodes that should be direct children of root
            let parent_id = if is_root {
                0  // Only the root node has parent_id 0
            } else {
                match &node.parent_name {
                    Some(parent) if !parent.is_empty() => *name_to_id.get(parent).unwrap_or(&root_id),
                    _ => root_id  // If no parent is specified, attach to root
                }
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
                            // Map accession numbers to build IDs
                            let build_id = match build.as_str() {
                                "CM000686.2" => "GRCh38",
                                "NC_000024.10" => "GRCh38",
                                "NC_060948.1" => "T2T-CHM13v2.0",
                                "CP086569.2" => "T2T-CHM13v2.0",
                                _ => build.as_str(),
                            };

                            (build_id.to_string(), LociCoordinate {
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

            all_nodes.insert(
                haplogroup_id.to_string(),
                HaplogroupNode {
                    haplogroup_id,
                    parent_id,
                    name: node.name,
                    is_root,
                    loci,
                    children: vec![],  // Will be populated later
                },
            );
        }

        // After creating all_nodes, collect the child relationships separately
        let mut children_updates: Vec<(u32, u32)> = Vec::new(); // (parent_id, child_id)
        for node in all_nodes.values() {
            if !node.is_root {
                children_updates.push((node.parent_id, node.haplogroup_id));
            }
        }

        // Now update the children vectors
        for (parent_id, child_id) in children_updates {
            if let Some(parent) = all_nodes.get_mut(&parent_id.to_string()) {
                parent.children.push(child_id);
            }
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