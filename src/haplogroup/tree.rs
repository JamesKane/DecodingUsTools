use crate::haplogroup::types::{Haplogroup, HaplogroupResult, LociType, Locus};
use crate::utils::cache::{TreeCache, TreeType};
use std::collections::HashMap;
use indicatif::{ProgressBar, ProgressStyle};
use crate::cli;

pub(crate) fn load_tree(tree_type: TreeType, provider: crate::cli::TreeProvider) -> Result<Haplogroup, Box<dyn std::error::Error>> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    progress.set_message("Initializing tree cache...");

    // Get tree from cache
    let tree_cache = TreeCache::new(tree_type, provider.clone())?;

    progress.set_message("Fetching haplogroup tree...");
    let tree = tree_cache.get_tree().map_err(|e| {
        format!("Failed to get haplogroup tree: {}", e)
    })?;

    progress.set_message("Building tree structure...");
    println!("Total nodes in tree: {}", tree.all_nodes.len());

    // Find root node based on provider
    let root_node = match provider {
        cli::TreeProvider::DecodingUs => {
            // For DecodingUs, use is_root flag
            tree.all_nodes.values()
                .find(|node| node.is_root)
                .ok_or("No node marked as root found in DecodingUs tree")?
        },
        cli::TreeProvider::FTDNA => {
            // For FTDNA, use parent_id == 0
            let root = tree.all_nodes.values()
                .find(|node| node.parent_id == 0)
                .ok_or("No root node found in FTDNA tree")?;
            if tree.all_nodes.values().filter(|node| node.parent_id == 0).count() > 1 {
                return Err("Multiple root nodes found in FTDNA tree".into());
            }
            root
        }
    };

    println!("Found root node: id={}, name='{}'",
             root_node.haplogroup_id, root_node.name);

    let haplogroup_tree = tree_cache
        .provider
        .build_tree(&tree, root_node.haplogroup_id, tree_type)
        .ok_or("Failed to build tree")?;

    progress.finish_with_message("Tree loaded successfully");
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
