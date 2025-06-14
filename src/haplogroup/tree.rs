use crate::haplogroup::types::{Haplogroup, HaplogroupResult, LociType, Locus};
use crate::utils::cache::{TreeCache, TreeType};
use std::collections::HashMap;

pub(crate) fn load_tree(tree_type: TreeType) -> Result<Haplogroup, Box<dyn std::error::Error>> {
    // Get tree from cache
    let tree_cache = TreeCache::new(tree_type)?;
    let tree = tree_cache.get_tree()?;

    let root_count = tree
        .all_nodes
        .values()
        .filter(|node| node.parent_id == 0)
        .count();
    if root_count > 1 {
        return Err("Multiple root nodes found in tree".into());
    }

    let root_node = tree
        .all_nodes
        .values()
        .find(|node| node.parent_id == 0)
        .ok_or("No root node found")?;

    let haplogroup_tree = tree_cache
        .provider
        .build_tree(&tree, root_node.haplogroup_id, tree_type)
        .ok_or("Failed to build tree")?;
    Ok(haplogroup_tree)
}

pub(crate) fn collect_snps<'a>(
    haplogroup: &'a Haplogroup,
    positions: &mut HashMap<u32, Vec<(&'a str, &'a Locus)>>,
    build_id: &str,
) {
    for locus in &haplogroup.loci {
        if let Some(coord) = locus.coordinates.get(build_id) {
            if matches!(locus.loci_type, LociType::SNP) {
                positions
                    .entry(coord.position)
                    .or_default()
                    .push((&haplogroup.name, locus));
            }
        }
    }

    for child in &haplogroup.children {
        collect_snps(child, positions, build_id);
    }
}

pub(crate) fn find_path_to_root(
    haplogroup: &Haplogroup,
    target_name: &str,
    scores: &[HaplogroupResult],
) -> Option<Vec<String>> {
    if haplogroup.name == target_name {
        return Some(vec![haplogroup.name.clone()]);
    }

    for child in &haplogroup.children {
        if let Some(mut path) = find_path_to_root(child, target_name, scores) {
            path.push(haplogroup.name.clone());
            return Some(path);
        }
    }

    None
}
