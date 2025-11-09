use crate::cli;
use crate::haplogroup::types::{Haplogroup, HaplogroupResult, LociCoordinate, LociType, Locus};
use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashMap;

/// Loads a haplogroup tree based on the specified `TreeType` and `TreeProvider`.
///
/// This function initializes a progress bar to provide feedback on the process of loading
/// the haplogroup tree. It retrieves the tree either from an existing cache or by fetching
/// and building the tree structure. The root node of the tree is determined based on the
/// specific `TreeProvider` implementation (e.g., `DecodingUs` or `FTDNA`). Once the root
/// node is identified, the tree is built and returned as a `Haplogroup` structure.
///
/// # Arguments
///
/// * `tree_type` - Specifies the type of tree to be loaded (`TreeType` enum).
/// * `provider` - Defines the source or method of loading the tree (via `TreeProvider`).
///
/// # Returns
///
/// This function returns a `Result`:
/// * `Ok(Haplogroup)` - A successfully loaded haplogroup tree.
/// * `Err(Box<dyn std::error::Error>)` - An error in case of failure during any step
///   (e.g., cache fetching, tree building, or invalid tree structure).
///
/// # Errors
///
/// Errors can occur due to the following reasons:
/// * Issues with initializing or using the tree cache.
/// * Failure to fetch the tree from the specified provider.
/// * Missing or invalid root node based on the `TreeProvider`.
/// * Errors during the tree-building phase (e.g., malformed or incomplete data).
///
/// # Example
///
/// ```no_run
/// use crate::cli::{TreeType, TreeProvider};
///
/// let tree_type = TreeType::Binary;
/// let provider = TreeProvider::DecodingUs;
///
/// match load_tree(tree_type, provider) {
///     Ok(tree) => println!("Successfully loaded tree with root: {}", tree.root),
///     Err(e) => eprintln!("Failed to load tree: {}", e),
/// }
/// ```
///
/// # Notes
///
/// * Progress feedback is displayed using a spinner for user interaction.
/// * For `DecodingUs`, the root node is determined by an `is_root` flag in the data.
/// * For `FTDNA`, the root node is identified by a `parent_id` of `0`.
/// * If multiple root nodes are found for the `FTDNA` provider, an error is returned.
///
/// # Dependencies
///
/// This function relies on the `TreeCache` for caching mechanisms and assumes a valid
/// `TreeProvider` interface for building the tree.
///
/// # See Also
///
/// * `TreeCache` - Provides caching functionality for trees.
/// * `TreeProvider` - Specifies logic for fetching or building tree data based on source.
pub(crate) fn load_tree(
    tree_type: TreeType,
    provider: crate::cli::TreeProvider,
) -> Result<Haplogroup, Box<dyn std::error::Error>> {
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
    let tree = tree_cache
        .get_tree()
        .map_err(|e| format!("Failed to get haplogroup tree: {}", e))?;

    progress.set_message("Building tree structure...");
    println!("Total nodes in tree: {}", tree.all_nodes.len());

    // Find root node based on provider
    let root_node = match provider {
        cli::TreeProvider::DecodingUs => {
            // For DecodingUs, use is_root flag
            tree.all_nodes
                .values()
                .find(|node| node.is_root)
                .ok_or("No node marked as root found in DecodingUs tree")?
        }
        cli::TreeProvider::FTDNA => {
            // For FTDNA, use parent_id == 0
            let root = tree
                .all_nodes
                .values()
                .find(|node| node.parent_id == 0)
                .ok_or("No root node found in FTDNA tree")?;
            if tree
                .all_nodes
                .values()
                .filter(|node| node.parent_id == 0)
                .count()
                > 1
            {
                return Err("Multiple root nodes found in FTDNA tree".into());
            }
            root
        }
    };

    println!(
        "Found root node: id={}, name='{}'",
        root_node.haplogroup_id, root_node.name
    );

    let haplogroup_tree = tree_cache
        .provider
        .build_tree(&tree, root_node.haplogroup_id, tree_type)
        .ok_or("Failed to build tree")?;

    progress.finish_with_message("Tree loaded successfully");
    Ok(haplogroup_tree)
}

/// Collects SNPs (Single Nucleotide Polymorphisms) from a phylogenetic tree of haplogroups and
/// organizes them into a hashmap indexed by position.
///
/// This function traverses the specified haplogroup tree and its descendants, collecting
/// information about loci of type `SNP` for a specified genomic build. The output is stored
/// in a hashmap where the key is the genomic position, and the value is a vector of tuples
/// containing the name of the haplogroup and a reference to the corresponding locus.
///
/// # Parameters
///
/// * `tree`: A reference to the root of the `Haplogroup` tree to traverse.
/// * `positions`: A mutable reference to a `HashMap` where the collected SNP data will be stored.  
///    The key is the genomic position (`u32`), and the value is a `Vec` of tuples. Each tuple
///    contains:
///   - A reference to the name of the haplogroup (`&'a str`).
///   - A reference to the corresponding `Locus` object (`&'a Locus`).
/// * `build_id`: A string slice (`&str`) that specifies the genomic build ID. Only loci with
///    coordinates corresponding to this build will be collected.
///
/// # Debugging
///
/// The function allows debugging using the following environment variables:
///
/// - `DECODINGUS_DEBUG_TREE`: If set to `1`, enables debug output for all processing steps in the tree.
/// - `DECODINGUS_DEBUG_NODE`: If set to a non-empty string, allows debugging specific nodes in the tree.
/// - `DECODINGUS_DEBUG_SITE`: If set to a specific locus name, enables debug output for that locus across all haplogroups.
///
/// When debugging is enabled, diagnostic messages are printed using `eprintln!`, providing details
/// about the haplogroup, locus, position, and alleles for each SNP.
///
/// # Notes
///
/// - Only loci of type `LociType::SNP` are collected.
/// - If a locus does not have coordinates for the specified `build_id`, it will be omitted from the
///   results. However, if debugging is enabled, a warning message will be logged.
///
/// # Example
///
/// ```rust
/// use std::collections::HashMap;
///
/// let mut positions = HashMap::new();
/// let root_haplogroup = get_sample_haplogroup_tree(); // Assume this function provides your tree.
/// collect_snps(&root_haplogroup, &mut positions, "GRCh37");
///
/// // Now `positions` contains SNP data, indexed by position.
/// for (pos, loci) in &positions {
///     println!("Position: {}", pos);
///     for (haplogroup, locus) in loci {
///         println!("  Haplogroup: {}, Locus: {}", haplogroup, locus.name);
///     }
/// }
/// ```
///
/// # Recursion
///
/// This function works recursively, processing all child `Haplogroup` nodes in the tree, ensuring
/// that all descendant SNPs are collected into the `positions` map.
///
/// # Performance
///
/// Depending on the size of the tree and the number of loci, this function can be computationally
/// intensive. Custom filters using the environment variables can help limit the scope for debugging.
///
/// # Panics
///
/// This function does not explicitly panic, but care should be taken to ensure:
/// - The `tree` structure is valid and properly initialized.
/// - The `build_id` corresponds to valid locus coordinate data within `tree`.
///
pub(crate) fn collect_snps<'a>(
    tree: &'a Haplogroup,
    positions: &mut HashMap<u32, Vec<(&'a str, &'a Locus)>>,
    build_id: &str,
) {
    let debug = std::env::var("DECODINGUS_DEBUG_TREE").ok().as_deref() == Some("1");
    let dbg_node = std::env::var("DECODINGUS_DEBUG_NODE")
        .ok()
        .filter(|v| !v.is_empty());

    // Debug specific locus name across all haplogroups
    let debug_locus_name = std::env::var("DECODINGUS_DEBUG_SITE").ok();

    for locus in &tree.loci {
        if matches!(locus.loci_type, LociType::SNP) {
            if let Some(coord) = locus.coordinates.get(build_id) {
                // Check if this locus matches the debug filter
                let is_debug_locus = debug_locus_name
                    .as_ref()
                    .map(|name| locus.name.contains(name))
                    .unwrap_or(false);

                if debug || is_debug_locus {
                    eprintln!("[tree.collect] haplogroup='{}' locus='{}' at pos={} (build={}) anc='{}' der='{}'",
                              tree.name, locus.name, coord.position, build_id,
                              coord.ancestral, coord.derived);
                }

                positions
                    .entry(coord.position)
                    .or_insert_with(Vec::new)
                    .push((&tree.name, locus));
            } else if debug {
                eprintln!(
                    "[tree.collect] {} locus {} MISSING coordinate for build {}",
                    tree.name, locus.name, build_id
                );
            }
        }
    }

    for child in &tree.children {
        collect_snps(child, positions, build_id);
    }
}

/// Finds the path from a specific haplogroup to the root of a hierarchy.
///
/// This function recursively traverses the haplogroup tree starting from the given haplogroup,
/// searching for the target haplogroup identified by `target_name`. If the target haplogroup
/// is found, the function constructs and returns a path (vector of haplogroup names) from
/// the target haplogroup to the root. The `scores` parameter is passed but not utilized in
/// the current implementation, suggesting it might be for future functionality or ignored.
///
/// # Arguments
///
/// * `haplogroup` - A reference to the current `Haplogroup` node being traversed.
/// * `target_name` - The name of the target haplogroup for which the path to the root is sought.
/// * `scores` - A slice of `HaplogroupResult` values, currently unused.
///
/// # Returns
///
/// * `Some(Vec<String>)` - If the target haplogroup is found, returns a vector of strings
///   representing the path of haplogroup names from the target to the root, with the target
///   being the first element and the root being the last.
/// * `None` - If the target haplogroup is not found within the hierarchy.
///
/// # Example
///
/// ```
/// // Sample data setup (assuming the appropriate structures and implementations are in place)
/// let haplogroup = Haplogroup {
///     name: "root".to_string(),
///     children: vec![
///         Haplogroup {
///             name: "child1".to_string(),
///             children: vec![
///                 Haplogroup {
///                     name: "target".to_string(),
///                     children: vec![],
///                 },
///             ],
///         },
///         Haplogroup {
///             name: "child2".to_string(),
///             children: vec![],
///         },
///     ],
/// };
///
/// let path = find_path_to_root(&haplogroup, "target", &[]);
/// assert_eq!(path, Some(vec!["target".to_string(), "child1".to_string(), "root".to_string()]));
/// ```
///
/// # Note
///
/// This function assumes that the haplogroup hierarchy is a tree (acyclic structure),
/// meaning there are no circular relationships between haplogroups. If such relationships exist,
/// this function may result in a stack overflow due to infinite recursion.
///
/// Additionally, the `scores` parameter is currently unused and could be removed or utilized
/// in the future to affect the logic of the function.
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

use crate::types::ReferenceGenome;
use crate::utils::liftover::Liftover;

/// Project a haplogroup tree's SNP coordinates from `src_build` into `dst_genome` build.
/// Adds/overwrites entries keyed by `dst_genome.name()` while preserving existing entries,
/// and omits unmappable loci in the destination build.
pub(crate) fn project_tree_to_build(
    tree: &Haplogroup,
    src_build: &str,
    dst_genome: ReferenceGenome,
    lifter: &Liftover,
) -> Haplogroup {
    let debug = std::env::var("DECODINGUS_DEBUG_TREE").ok().as_deref() == Some("1");

    fn project_locus(
        locus: &Locus,
        src_build: &str,
        dst_build_id: &str,
        lifter: &Liftover,
        debug: bool,
        parent_name: &str,
    ) -> Locus {
        let mut new_locus = locus.clone();
        if let Some(coord) = locus.coordinates.get(src_build) {
            // Only handle SNPs here; INDELs remain unmodified.
            // Use the new atomic liftover API to map position + alleles with strand awareness.
            if let Some(mapped) = lifter.map_snp_alleles(
                &coord.chromosome,
                coord.position,
                &coord.ancestral,
                &coord.derived,
            ) {
                let dbg_node = std::env::var("DECODINGUS_DEBUG_NODE")
                    .ok()
                    .filter(|v| !v.is_empty());
                if debug
                    && dbg_node
                        .as_ref()
                        .map(|f| parent_name.contains(f))
                        .unwrap_or(false)
                {
                    eprintln!(
                        "[tree.debug] {} locus {}: {} {}:{} ({}->{}) -> {} {}:{} ({}->{}){}",
                        parent_name,
                        locus.name,
                        src_build,
                        coord.chromosome,
                        coord.position,
                        coord.ancestral,
                        coord.derived,
                        dst_build_id,
                        mapped.contig,
                        mapped.position,
                        mapped.ancestral,
                        mapped.derived,
                        if mapped.was_reverse { " [RC]" } else { "" }
                    );
                }

                let projected = LociCoordinate {
                    position: mapped.position,
                    chromosome: mapped.contig,
                    ancestral: mapped.ancestral,
                    derived: mapped.derived,
                };
                let mut coords = new_locus.coordinates.clone();
                coords.insert(dst_build_id.to_string(), projected);
                new_locus.coordinates = coords;
            } else if debug {
                let dbg_node = std::env::var("DECODINGUS_DEBUG_NODE")
                    .ok()
                    .filter(|v| !v.is_empty());
                if dbg_node
                    .as_ref()
                    .map(|f| parent_name.contains(f))
                    .unwrap_or(false)
                {
                    eprintln!(
                        "[tree.debug] {} locus {}: FAILED to liftover from {} {}:{}",
                        parent_name, locus.name, src_build, coord.chromosome, coord.position
                    );
                }
            }
        }
        new_locus
    }

    let dst_build_id = dst_genome.name();

    let mut projected = Haplogroup {
        name: tree.name.clone(),
        parent: tree.parent.clone(),
        loci: tree
            .loci
            .iter()
            .map(|l| project_locus(l, src_build, dst_build_id, lifter, debug, &tree.name))
            .collect(),
        children: Vec::new(),
    };

    projected.children = tree
        .children
        .iter()
        .map(|c| project_tree_to_build(c, src_build, dst_genome.clone(), lifter))
        .collect();

    projected
}
