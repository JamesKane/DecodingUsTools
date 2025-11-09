pub mod caller;
mod scoring;
mod tree;
pub(crate) mod types;
mod validation;

use crate::haplogroup::types::{Haplogroup, HaplogroupResult, LociType, Locus};
use crate::types::ReferenceGenome;
use crate::utils::cache::TreeType;
use crate::utils::liftover::Liftover;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::{self, Read as BamRead};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Analyzes the haplogroup of genomic data using a given BAM/CRAM file and tree structure, and produces a scored output.
///
/// # Parameters
/// - `bam_file`: Path to the input BAM/CRAM file containing the sequencing data.
/// - `reference_file`: Path to the reference genome file required for analyses, particularly for CRAM file validation.
/// - `output_file`: Path to the desired output file where the results will be stored.
/// - `min_depth`: Minimum read depth required in the BAM/CRAM file to consider a valid SNP.
/// - `min_quality`: Minimum base quality for SNPs used in scoring the haplogroups.
/// - `tree_type`: Type of tree to use for analyzing haplogroups. This could be Y-DNA (`TreeType::YDNA`) or mitochondrial DNA (`TreeType::MTDNA`).
/// - `provider`: A provider for the haplogroup tree (e.g., FTDNA).
/// - `show_snps`: Boolean flag to indicate whether to output the list of SNPs used during the analysis.
/// - `include_off_path`: Boolean flag indicating whether to include haplogroup loci off the main phylogenetic path while scoring haplogroups.
///
/// # Return Value
/// - Returns a `Result<()>` on successful completion. Errors are wrapped in a `Box<dyn Error>` for broad compatibility.
///
/// # Workflow
/// 1. **Reference Validation**: Validates the BAM file against the selected reference genome.
///    - If the input file is a CRAM file, the reference genome location is set using the `REF_PATH` environment variable.
/// 2. **Tree Loading**: Loads the haplogroup tree structure according to the specified type (`YDNA` or `MTDNA`) and extracts the build information from both the tree and BAM file.
/// 3. **Coordinate Liftover**: Converts tree SNP coordinates to the BAM's genomic build if necessary, creating a "projected tree" compatible with the data being analyzed.
/// 4. **SNP Collection**: Collects relevant SNP positions from the projected tree and evaluates corresponding SNPs from the BAM file, filtering them by the minimum depth and quality thresholds.
/// 5. **Debugging (Optional)**:
///    - Allows optional debug output for tree projection, specific haplogroups, and tracked SNP sites controlled by environment variables.
///    - Debug output shows detailed mappings, failed projections, and tracked loci for deeper introspection.
/// 6. **Haplogroup Scoring**: Computes scores for haplogroups based on the filtered SNP calls and constructs an ordered list of scores.
/// 7. **Output**: The results, including haplogroup scores and optional SNP details, are saved to the specified `output_file`.
///
/// # Notes
/// - Environmental variables `DECODINGUS_DEBUG_PROJECTION`, `DECODINGUS_DEBUG_NODE`, and `DECODINGUS_DEBUG_SITE` can be used for debugging purposes:
///   - `DECODINGUS_DEBUG_PROJECTION=1`: Enables verbose debug output of tree projection details.
///   - `DECODINGUS_DEBUG_NODE`: Filters debug output for specific haplogroups.
///   - `DECODINGUS_DEBUG_SITE`: Tracks specific genomic positions in SNP debugging.
/// - Haplogroup scoring and projections are normalized to a consistent coordinate system. Liftover operations between builds are performed when necessary.
/// - The function supports multiple output customizations through its flags, simplifying downstream analysis.
///
/// # Errors
/// - Errors can occur due to:
///   - Missing or invalid input files (`bam_file`, `reference_file`).
///   - Failure to load or project the haplogroup tree.
///   - Invalid SNP data (e.g., insufficient read depth or base quality).
///   - Issues during scoring or output generation.
///
/// # Example
/// ```rust
/// use crate::{analyze_haplogroup, TreeType};
///
/// let result = analyze_haplogroup(
///     "input.bam".to_string(),
///     "reference.fa".to_string(),
///     "output.txt".to_string(),
///     10,               // Minimum depth
///     30,               // Minimum quality
///     TreeType::YDNA,   // Tree type
///     TreeProvider::FTDNA, // Tree source provider
///     true,             // Show SNPs
///     true              // Include off-path loci
/// );
///
/// match result {
///     Ok(()) => println!("Haplogroup analysis completed successfully."),
///     Err(e) => eprintln!("Error during haplogroup analysis: {}", e),
/// }
/// ```
pub fn analyze_haplogroup(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    tree_type: TreeType,
    provider: crate::cli::TreeProvider,
    show_snps: bool,
    include_off_path: bool,
) -> Result<(), Box<dyn Error>> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );

    progress.set_message("Validating BAM reference genome...");

    // Set reference for CRAM files
    if bam_file.ends_with(".cram") {
        std::env::set_var("REF_PATH", &reference_file);
    }

    // Open BAM/CRAM. For BAM we don't require an index here.
    let mut bam = bam::Reader::from_path(&bam_file)?;
    let header = bam.header().clone();
    let (genome, chromosome) = validation::validate_reference(&bam, tree_type)?;

    let haplogroup_tree = tree::load_tree(tree_type, provider)?;

    // Determine the build used by the tree and the build of the BAM
    let bam_build = genome; // ReferenceGenome
    let src_build_id = match tree_type {
        TreeType::YDNA => "GRCh38", // FTDNA provides GRCh38
        TreeType::MTDNA => "rCRS",
    };

    // Build liftover from tree's build to BAM build
    let lifter = if tree_type == TreeType::MTDNA {
        // rCRS identity mapping with chromosome aliasing handled downstream
        Liftover::identity()
    } else {
        // Parse src_build_id to ReferenceGenome
        let source = match src_build_id {
            "GRCh38" => ReferenceGenome::GRCh38,
            "GRCh37" => ReferenceGenome::GRCh37,
            _ => ReferenceGenome::GRCh38,
        };
        Liftover::load_or_fetch(source.clone(), bam_build.clone())?
    };

    // Project the tree to the BAM's build; this normalizes positions and alleles (strand-aware)
    let projected_tree =
        tree::project_tree_to_build(&haplogroup_tree, src_build_id, bam_build.clone(), &lifter);

    // Debug: verify projection worked for key haplogroups
    let debug_projection =
        std::env::var("DECODINGUS_DEBUG_PROJECTION").ok().as_deref() == Some("1");
    if debug_projection {
        eprintln!(
            "[debug] Tree projection: {} -> {}",
            src_build_id,
            bam_build.name()
        );

        // If a specific haplogroup filter is provided, show projection details for it
        if let Some(filter) = std::env::var("DECODINGUS_DEBUG_NODE")
            .ok()
            .filter(|s| !s.is_empty())
        {
            if let Some(original_hg) = find_haplogroup(&haplogroup_tree, &filter) {
                if let Some(projected_hg) = find_haplogroup(&projected_tree, &filter) {
                    eprintln!("[debug] Haplogroup {} projection check:", filter);
                    for locus in &original_hg.loci {
                        if let Some(orig_coord) = locus.coordinates.get(src_build_id) {
                            eprintln!(
                                "[debug]   Original ({}): {} pos={} anc={} der={}",
                                src_build_id,
                                locus.name,
                                orig_coord.position,
                                orig_coord.ancestral,
                                orig_coord.derived
                            );
                        }
                    }
                    for locus in &projected_hg.loci {
                        if let Some(proj_coord) = locus.coordinates.get(bam_build.name()) {
                            eprintln!(
                                "[debug]   Projected ({}): {} pos={} anc={} der={}",
                                bam_build.name(),
                                locus.name,
                                proj_coord.position,
                                proj_coord.ancestral,
                                proj_coord.derived
                            );
                        } else {
                            eprintln!(
                                "[debug]   Projected ({}): {} MISSING (liftover failed)",
                                bam_build.name(),
                                locus.name
                            );
                        }
                    }
                }
            }
        }
    }

    // Collect positions directly from the projected tree in BAM build coordinates
    let build_id = bam_build.name();
    let mut positions_bam: HashMap<u32, Vec<(&str, &Locus)>> = HashMap::new();
    tree::collect_snps(&projected_tree, &mut positions_bam, build_id);

    if debug_projection {
        eprintln!(
            "[debug] Collected {} unique positions from projected tree for build '{}'",
            positions_bam.len(),
            build_id
        );
        // Show a sample of collected positions
        let mut sample_positions: Vec<u32> = positions_bam.keys().copied().collect();
        sample_positions.sort();
        eprintln!(
            "[debug] Sample positions collected: {:?}",
            sample_positions.iter().take(10).collect::<Vec<_>>()
        );
    }

    // Optional: focus debug on a specific site by locus name (e.g., M526)
    let debug_site = std::env::var("DECODINGUS_DEBUG_SITE").ok();
    let mut debug_tree_sites: Vec<u32> = Vec::new();
    if let Some(ref site_tag) = debug_site {
        for (pos, entries) in positions_bam.iter() {
            if entries
                .iter()
                .any(|(_, locus)| locus.name.contains(site_tag))
            {
                debug_tree_sites.push(*pos);
            }
        }
        if !debug_tree_sites.is_empty() {
            eprintln!(
                "[debug] Tracking sites by locus tag '{}': {:?}",
                site_tag, debug_tree_sites
            );
        }
    }

    let mut snp_calls_bam: HashMap<u32, (char, u32, f64)> = HashMap::new();

    // Collect SNP calls directly at BAM coordinates
    caller::collect_snp_calls(
        min_depth,
        min_quality,
        &mut bam,
        &header,
        build_id.to_string(),
        chromosome.clone(),
        &mut positions_bam,
        &mut snp_calls_bam,
    )?;

    // No remapping needed; scoring consumes BAM-build coordinates against the projected tree
    let snp_calls = snp_calls_bam;

    // Debug: dump snp_calls for any sites selected via DECODINGUS_DEBUG_SITE
    if (debug_projection || debug_site.is_some()) && !debug_tree_sites.is_empty() {
        for pos in debug_tree_sites.iter() {
            if let Some((base, depth, freq)) = snp_calls.get(pos) {
                eprintln!(
                    "[debug] snp_calls at tracked position {}: base='{}' depth={} freq={:.3}",
                    pos, base, depth, freq
                );
            } else {
                eprintln!(
                    "[debug] snp_calls at tracked position {}: NO CALL FOUND",
                    pos
                );
            }
            if let Some(entries) = positions_bam.get(pos) {
                eprintln!(
                    "[debug] positions_bam at tracked position {} has {} entries:",
                    pos,
                    entries.len()
                );
                for (hg_name, locus) in entries {
                    eprintln!(
                        "[debug]   - haplogroup='{}' locus='{}'",
                        hg_name, locus.name
                    );
                }
            } else {
                eprintln!(
                    "[debug] positions_bam at tracked position {}: NOT TRACKED",
                    pos
                );
            }
        }
    }

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups with projected tree and BAM build id
    let mut scores = Vec::new();
    scoring::calculate_haplogroup_score(
        &projected_tree,
        &snp_calls,
        &mut scores,
        None,
        0,
        build_id,
    );

    let mut ordered_scores = collect_scored_paths(scores, &projected_tree, include_off_path);

    // Optional debug dump of top haplogroups (env-gated)
    if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
        eprintln!("[debug] Top haplogroups (name, score, matches, ancestral, mismatching, no_calls, total, cumulative, depth)");
        for r in ordered_scores.iter().take(25) {
            eprintln!(
                "[debug] {}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.name,
                r.score,
                r.matching_snps,
                r.ancestral_matches,
                r.mismatching_snps,
                r.no_calls,
                r.total_snps,
                r.cumulative_snps,
                r.depth
            );
        }
        // Optional focused view filtered by DECODINGUS_DEBUG_NODE
        if let Some(filter) = std::env::var("DECODINGUS_DEBUG_NODE")
            .ok()
            .filter(|s| !s.is_empty())
        {
            let entries: Vec<_> = ordered_scores
                .iter()
                .filter(|r| r.name.contains(&filter))
                .take(25)
                .collect();
            if !entries.is_empty() {
                eprintln!("[debug] Filtered snapshot (filter='{}'):", filter);
                for r in entries {
                    eprintln!(
                        "[debug] {}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        r.name,
                        r.score,
                        r.matching_snps,
                        r.ancestral_matches,
                        r.mismatching_snps,
                        r.no_calls,
                        r.total_snps,
                        r.cumulative_snps,
                        r.depth
                    );
                }
            }
        }
    }

    // Use ordered_scores for writing the report
    let mut writer = BufWriter::new(File::create(output_file)?);

    if show_snps {
        writeln!(
            writer,
            "Haplogroup\tScore\tMatching_SNPs\tMismatching_SNPs\tAncestral_Matches\tNo_Calls\tTotal_SNPs\tCumulative_SNPs\tDepth\tMatching_SNP_Details\tMismatching_SNP_Details\tNo_Call_Details"
        )?;
    } else {
        writeln!(
            writer,
            "Haplogroup\tScore\tMatching_SNPs\tMismatching_SNPs\tAncestral_Matches\tNo_Calls\tTotal_SNPs\tCumulative_SNPs\tDepth"
        )?;
    }

    for result in ordered_scores {
        // Recompute branch-local metrics for reporting to avoid cumulative leakage
        let (b_match, b_mismatch, b_anc, b_nocall, b_total) =
            compute_branch_local_metrics(&result.name, &projected_tree, &snp_calls, build_id);

        if show_snps {
            let (matching_snps, mismatching_snps, no_call_snps) =
                get_snp_details(&result.name, &projected_tree, &snp_calls, build_id);
            writeln!(
                writer,
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.name,
                result.score,
                b_match,
                b_mismatch,
                b_anc,
                b_nocall,
                b_total,
                result.cumulative_snps,
                result.depth,
                matching_snps,
                mismatching_snps,
                no_call_snps
            )?;
        } else {
            writeln!(
                writer,
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.name,
                result.score,
                b_match,
                b_mismatch,
                b_anc,
                b_nocall,
                b_total,
                result.cumulative_snps,
                result.depth
            )?;
        }
    }

    progress.finish_with_message("Analysis complete!");
    Ok(())
}

/// Retrieves detailed information about the matching, mismatching, and uncalled SNPs (Single Nucleotide Polymorphisms)
/// for a given haplogroup in a phylogenetic tree.
///
/// # Arguments
///
/// * `haplogroup_name` - A `&str` containing the name of the haplogroup to investigate.
/// * `tree` - A reference to the phylogenetic tree of type `Haplogroup`, where haplogroups and their loci are managed.
/// * `snp_calls` - A reference to a `HashMap` where the keys are SNP positions (u32) and the values are tuples
///   containing:
///   - `char`: The base that was called at this position.
///   - `u32`: The read depth for the SNP call.
///   - `f64`: The allele frequency for the SNP call.
/// * `build_id` - A `&str` identifying the genome build/reference version to use when determining SNP positions.
///
/// # Returns
///
/// A tuple of three `String`s:
/// 1. `matching`: Concatenated pairs of locus names and positions (in the format `name:position`) where the
///     base in the SNP call matches the derived base at that position, separated by semicolons.
/// 2. `mismatching`: Concatenated pairs of locus names and positions (in the format `name:position`) where the
///     base in the SNP call does not match the derived base, separated by semicolons.
/// 3. `no_calls`: Concatenated pairs of locus names and positions (in the format `name:position`) where no SNP
///     call was made for the corresponding locus, separated by semicolons.
///
/// # Behavior
///
/// - If the specified `haplogroup_name` is found within the `tree`, the loci for the haplogroup are iterated over.
/// - For each locus, its genomic coordinate in the specified `build_id` is retrieved.
/// - It checks the SNP call for that specific genomic position in the `snp_calls` map and categorizes the locus
///   as `matching`, `mismatching`, or `no_calls` depending on whether the called base matches the derived base.
///
/// # Notes
///
/// - This function assumes that the input data structures are properly initialized and contain valid information.
/// - The derived base is extracted from the `derived` string of the `coordinates` field of a `locus` using its
///   first character (`chars().next().unwrap()`).
/// - If a matching haplogroup is not found in the tree, empty strings are returned for all three components of the tuple.
///
/// # Example
///
/// ```
/// use std::collections::HashMap;
///
/// // Assuming the existence of properly initialized `tree` and `snp_calls`
/// let haplogroup_name = "H1";
/// let build_id = "GRCh37";
/// let (matching, mismatching, no_calls) = get_snp_details(haplogroup_name, &tree, &snp_calls, build_id);
///
/// println!("Matching SNPs: {}", matching);
/// println!("Mismatching SNPs: {}", mismatching);
/// println!("No call SNPs: {}", no_calls);
/// ```
///
/// # Panics
///
/// - The function will panic if the derived base string (`locus.coordinates[build_id].derived`) is empty,
///   as the `unwrap()` is directly applied to retrieve the first character.
///
/// # Dependencies
///
/// The function relies on an external structure named `Haplogroup` and the `find_haplogroup` function, which are
/// expected to define how to navigate the phylogenetic tree and retrieve data for a given haplogroup.
fn get_snp_details(
    haplogroup_name: &str,
    tree: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    build_id: &str,
) -> (String, String, String) {
    let mut matching = Vec::new();
    let mut mismatching = Vec::new();
    let mut no_calls = Vec::new();

    if let Some(node) = find_haplogroup(tree, haplogroup_name) {
        for locus in &node.loci {
            if let Some(coord) = locus.coordinates.get(build_id) {
                match snp_calls.get(&coord.position) {
                    Some((called_base, depth, freq)) => {
                        let derived_base = coord.derived.chars().next().unwrap();
                        if *called_base == derived_base {
                            matching.push(format!("{}:{}", locus.name, coord.position));
                        } else {
                            mismatching.push(format!("{}:{}", locus.name, coord.position));
                        }
                    }
                    None => {
                        no_calls.push(format!("{}:{}", locus.name, coord.position));
                    }
                }
            }
        }
    }

    (
        matching.join(";"),
        mismatching.join(";"),
        no_calls.join(";"),
    )
}

/// Recursively searches for a haplogroup by name within a hierarchical haplogroup tree.
///
/// # Parameters
/// - `tree`: A reference to the root `Haplogroup` node from where the search begins.
/// - `name`: A string slice representing the name of the haplogroup to search for.
///
/// # Returns
/// - `Option<&'a Haplogroup>`:
///   - Returns `Some(&Haplogroup)` if a haplogroup with the specified name is found in the tree.
///   - Returns `None` if no haplogroup with the specified name exists within the tree.
///
/// # Lifetimes
/// - `'a`: The lifetime of the `tree` reference ensures the returned reference, if any, is valid as long as the input tree is valid.
///
/// # Example
/// ```
/// let tree = Haplogroup {
///     name: "Root".to_string(),
///     children: vec![
///         Haplogroup {
///             name: "A".to_string(),
///             children: vec![],
///         },
///         Haplogroup {
///             name: "B".to_string(),
///             children: vec![
///                 Haplogroup {
///                     name: "B1".to_string(),
///                     children: vec![],
///                 }
///             ],
///         }
///     ],
/// };
///
/// if let Some(haplogroup) = find_haplogroup(&tree, "B1") {
///     println!("Found haplogroup: {}", haplogroup.name);
/// } else {
///     println!("Haplogroup not found.");
/// }
/// ```
///
/// # Notes
/// - This function performs a depth-first search through the tree structure.
/// - It assumes the `Haplogroup` structure has fields `name` (a `String`) and `children` (a `Vec<Haplogroup>`).
fn find_haplogroup<'a>(tree: &'a Haplogroup, name: &str) -> Option<&'a Haplogroup> {
    if tree.name == name {
        return Some(tree);
    }
    for child in &tree.children {
        if let Some(found) = find_haplogroup(child, name) {
            return Some(found);
        }
    }
    None
}

/// Collects and ranks scored haplogroup paths based on SNP data and tree structure.
///
/// This function processes a list of scored haplogroups (`scores`) along with a haplogroup tree (`tree`),
/// and returns a deduplicated, ranked list of haplogroups. The ranking prioritizes terminal haplogroups
/// (leaves) as well as paths supported by derived evidence.
///
/// # Parameters
///
/// - `scores`: A vector of `HaplogroupResult` representing haplogroups and their associated scores.
/// - `tree`: A reference to the `Haplogroup` tree structure which defines the hierarchy of haplogroups.
/// - `include_off_path`: A boolean flag indicating whether to include haplogroups outside the primary path
/// to the terminal node for ranking purposes.
///
/// # Returns
///
/// - A vector of `HaplogroupResult` objects, ordered by their relevance based on scoring and tree path rules.
///
/// The function performs the following steps:
/// 1. Deduplicates the `scores` list by retaining the result with the highest score for each haplogroup.
/// 2. Filters to retain only terminal haplogroups (leaves) from the tree.
/// 3. Computes derived and ancestral path-based evidence for leaves.
/// 4. Removes leaves that lack sufficient derived evidence along their ancestral paths.
/// 5. Sorts the remaining results based on various criteria, including:
///    - Consistency scores (difference between matching and ancestral SNPs).
///    - Overall score.
///    - Cumulative SNP match counts.
/// 6. Applies further rules and heuristics to rank deeper leaves or nodes with fewer ancestral conflicts.
///
/// # Helper Functions
///
/// - `collect_leaf_names`: Recursively collects the names of all leaf nodes in the haplogroup tree.
/// - `find_path_to_root`: A utility function (assumed to be defined elsewhere) that returns the path from
///   a given haplogroup to the root of the tree.
/// - `pick_best_by_rules`: Traverses the tree using evidence-based rules and selects the most reliable terminal
///   haplogroup.
///
/// # Example
///
/// ```
/// // Assuming Haplogroup and HaplogroupResult structures are defined as:
/// struct Haplogroup {
///     name: String,
///     children: Vec<Haplogroup>,
/// }
///
/// struct HaplogroupResult {
///     name: String,
///     score: f64,
///     matching_snps: u32,
///     ancestral_matches: u32,
///     cumulative_snps: u32,
///     depth: u32,
/// }
///
/// let tree = Haplogroup {
///     name: "Root".to_string(),
///     children: vec![
///         Haplogroup {
///             name: "A".to_string(),
///             children: vec![],
///         },
///         Haplogroup {
///             name: "B".to_string(),
///             children: vec![
///                 Haplogroup {
///                     name: "B1".to_string(),
///                     children: vec![],
///                 },
///             ],
///         },
///     ],
/// };
///
/// let scores = vec![
///     HaplogroupResult {
///         name: "A".to_string(),
///         score: 95.0,
///         matching_snps: 50,
///         ancestral_matches: 10,
///         cumulative_snps: 60,
///         depth: 1,
///     },
///     HaplogroupResult {
///         name: "B1".to_string(),
///         score: 90.0,
///         matching_snps: 60,
///         ancestral_matches: 5,
///         cumulative_snps: 65,
///         depth: 2,
///     },
/// ];
///
/// let include_off_path = true;
/// let result = collect_scored_paths(scores, &tree, include_off_path);
///
/// assert_eq!(result.len(), 1);
/// assert_eq!(result[0].name, "B1");
/// ```
///
/// # Notes
///
/// - The function assumes the presence of additional helper functions like `find_path_to_root` for traversing
/// the tree. Ensure all dependencies are implemented or imported.
/// - Sorting and scoring criteria can be adjusted based on the specific requirements or datasets used.
fn collect_scored_paths(
    scores: Vec<HaplogroupResult>,
    tree: &Haplogroup,
    include_off_path: bool,
) -> Vec<HaplogroupResult> {
    let mut ordered_results = Vec::new();

    // Helper: collect leaf names
    fn collect_leaf_names(node: &Haplogroup, acc: &mut std::collections::HashSet<String>) {
        if node.children.is_empty() {
            acc.insert(node.name.clone());
        } else {
            for c in &node.children {
                collect_leaf_names(c, acc);
            }
        }
    }

    let mut leaf_names = std::collections::HashSet::new();
    collect_leaf_names(tree, &mut leaf_names);

    // First, deduplicate the results by keeping only the entry with the highest score for each haplogroup
    let mut unique_results: HashMap<String, HaplogroupResult> = HashMap::new();
    for result in scores {
        unique_results
            .entry(result.name.clone())
            .and_modify(|existing| {
                if result.score > existing.score {
                    *existing = result.clone();
                }
            })
            .or_insert(result);
    }

    // Convert to vec and filter with stricter criteria; keep primarily leaves for ranking
    let mut remaining: Vec<HaplogroupResult> = unique_results
        .values()
        .cloned()
        .filter(|result| {
            // Consider all leaves; do not pre-filter by clean/consistent to allow downstream evidence to drive selection
            leaf_names.contains(&result.name)
        })
        .collect();

    // Compute path-derived counts for each leaf to ensure we only consider leaves supported by any derived evidence along the path
    let mut path_derived: HashMap<String, u32> = HashMap::new();
    let mut path_ancestral: HashMap<String, u32> = HashMap::new();
    for leaf in &remaining {
        if let Some(path) = tree::find_path_to_root(tree, &leaf.name, &[]) {
            let mut d: u32 = 0;
            let mut a_m: u32 = 0;
            for name in path.iter() {
                if let Some(entry) = unique_results.get(name) {
                    d = d.saturating_add(entry.matching_snps);
                    a_m = a_m.saturating_add(entry.ancestral_matches);
                }
            }
            path_derived.insert(leaf.name.clone(), d);
            path_ancestral.insert(leaf.name.clone(), a_m);
        }
    }

    // Keep only leaves that have at least one derived call along their ancestral path
    remaining.retain(|leaf| path_derived.get(&leaf.name).copied().unwrap_or(0) >= 1);

    // Sort by consistency first (no-call tolerant), then by score, then by cumulative SNPs
    remaining.sort_by(|a, b| {
        let a_consistency = (a.matching_snps as i64) - 2 * (a.ancestral_matches as i64);
        let b_consistency = (b.matching_snps as i64) - 2 * (b.ancestral_matches as i64);
        b_consistency
            .cmp(&a_consistency)
            .then_with(|| {
                b.score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    });

    // Helper: comparator prioritizing deeper, consistent leaves with fewer ancestral conflicts
    let cmp = |a: &HaplogroupResult, b: &HaplogroupResult| {
        let a_consistency = (a.matching_snps as i64) - 2 * (a.ancestral_matches as i64);
        let b_consistency = (b.matching_snps as i64) - 2 * (b.ancestral_matches as i64);
        b.depth
            .cmp(&a.depth) // prefer deeper (closer to terminal)
            .then_with(|| b.ancestral_matches.cmp(&a.ancestral_matches)) // fewer ancestral preferred
            .then_with(|| b_consistency.cmp(&a_consistency))
            .then_with(|| {
                b.score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    };

    /// Picks the best node in a haplogroup tree based on specific rules and criteria.
    ///
    /// This function processes a node in a haplogroup tree and decides whether to
    /// stop at the current node or continue traversing its children based on
    /// SNP (Single-Nucleotide Polymorphism) matching rules, ancestral evidence, and
    /// conflict tracking. The goal is to determine the most appropriate haplogroup
    /// based on the provided data.
    ///
    /// # Arguments
    ///
    /// - `node`: A reference to the current haplogroup node being evaluated.
    /// - `map`: A `HashMap` containing haplogroup names as keys and their corresponding
    ///   `HaplogroupResult` as values. This provides SNP and ancestral evidence for nodes.
    /// - `leaves`: A vector containing all `HaplogroupResult` for leaf nodes in the tree.
    /// - `tree`: A reference to the root of the haplogroup tree, allowing access to its structure.
    /// - `consec_conflicts`: The number of consecutive ancestor blocks with conflicts encountered so far.
    ///
    /// # Returns
    ///
    /// Returns a reference to the most appropriate haplogroup node according to the
    /// defined decision rules.
    ///
    /// # Decision Rules
    ///
    /// 1. **Stopping Conditions**:
    ///    - If no entry exists in the `map` for the current node, the traversal stops here.
    ///    - If the node is a leaf (i.e., has no children), the traversal stops here.
    ///    - If three consecutive conflicting ancestor blocks are encountered, the traversal stops here.
    ///
    /// 2. **Special Case Handling**:
    ///    - If a SNP block contains more than five SNPs but only one derived call, the branch is treated
    ///      as ancestral.
    ///    - If the parent node and one of its direct children are strictly all-ancestral, no further
    ///      descent is allowed.
    ///
    /// 3. **Descent and Comparison**:
    ///    - Only consider descent into a child when there is sufficient derived evidence or any
    ///      downstream derived evidence from descendants.
    ///    - For each child node:
    ///        - Evaluate its suitability as the best candidate node based on various SNP metrics
    ///          and consistency rules.
    ///        - Consider grandchildren nodes only up to one level below the current node (e.g.,
    ///          immediate children of the current node's children).
    ///
    /// 4. **Leaf Node Selection**:
    ///    - Use a comparison function `cmp_leaf` to select the best leaf node in the subtree based on:
    ///        - Node depth.
    ///        - Ancestral matches.
    ///        - Derived SNP consistency metrics.
    ///        - Cumulative SNP scores.
    ///
    /// 5. **Traversal Limitation**:
    ///    - Only allow traversal one level deeper into the tree from the current node.
    ///
    /// # Comparison Function
    ///
    /// A helper comparison function `cmp_leaf` is used to rank leaf nodes based on:
    /// - Depth of the node in the tree.
    /// - Number of ancestral matches.
    /// - SNP consistency metrics, calculated as `(matching_snps - 2 * ancestral_matches)`.
    /// - Score of the node, considering SNP evidence.
    /// - Cumulative SNP count.
    ///
    /// # Behavior
    ///
    /// The function ensures that decisions are made using clear, evidence-based rules to
    /// avoid unnecessary conflicts or traversals. It balances strict ancestral evidence
    /// validation with the flexibility to explore better-fitting haplogroups deeper in the tree.
    ///
    /// # Example
    ///
    /// Suppose you have a haplogroup tree and a map of SNP evidence for various nodes. You want to
    /// determine the best haplogroup for a particular node based on the provided rules:
    ///
    /// ```rust
    /// let result = pick_best_by_rules(node, &haplogroup_map, &leaf_nodes, tree, 0);
    /// println!("Best haplogroup: {}", result.name);
    /// ```
    ///
    /// This will return the best haplogroup node that satisfies the rules, ensuring a consistent
    /// and reliable result.
    fn pick_best_by_rules<'a>(
        node: &'a Haplogroup,
        map: &HashMap<String, HaplogroupResult>,
        leaves: &Vec<HaplogroupResult>,
        tree: &'a Haplogroup,
        consec_conflicts: u32,
    ) -> &'a Haplogroup {
        let name = &node.name;
        let entry = map.get(name);
        // No entry -> stop here
        if entry.is_none() {
            return node;
        }
        let e = entry.unwrap();
        let total = e.total_snps;
        let mut derived = e.matching_snps;
        let ancestral = e.ancestral_matches;
        // When there is only one derived call in a block of more than five, treat the branch as ancestral
        if total > 5 && derived == 1 {
            derived = 0;
        }
        // Define conflict for succession tracking only for larger blocks (>4) with half-or-more ancestral
        let block_conflict = total > 4 && (ancestral as u32) >= ((total + 1) / 2);
        let consec = if block_conflict {
            consec_conflicts + 1
        } else {
            0
        };
        // If leaf, stop
        if node.children.is_empty() {
            return node;
        }
        // Terminate if three successive conflicting ancestor blocks
        if consec >= 3 {
            return node;
        }
        // Special case: only one phylo variant and it is ancestral -> inspect children rather than stopping
        let single_ancestral = total == 1 && derived == 0 && ancestral == 1;
        // Early termination: stop when parent and a direct child are strictly all-ancestral
        let parent_all_ancestral = total > 0 && derived == 0 && ancestral == total;
        if parent_all_ancestral {
            for child in &node.children {
                if let Some(ce) = map.get(&child.name) {
                    let c_total = ce.total_snps;
                    let c_derived = ce.matching_snps;
                    let c_ancestral = ce.ancestral_matches;
                    let child_all_ancestral =
                        c_total > 0 && c_derived == 0 && c_ancestral == c_total;
                    if child_all_ancestral {
                        return node;
                    }
                }
            }
        }
        // Compare leaves helper
        let cmp_leaf = |a: &HaplogroupResult, b: &HaplogroupResult| {
            let a_consistency = (a.matching_snps as i64) - 2 * (a.ancestral_matches as i64);
            let b_consistency = (b.matching_snps as i64) - 2 * (b.ancestral_matches as i64);
            b.depth
                .cmp(&a.depth)
                .then_with(|| b.ancestral_matches.cmp(&a.ancestral_matches))
                .then_with(|| b_consistency.cmp(&a_consistency))
                .then_with(|| {
                    b.score
                        .partial_cmp(&a.score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
        };
        // For each child, find best downstream leaf and whether any derived evidence exists along that subtree
        let mut best_child: Option<&Haplogroup> = None;
        let mut best_child_leaf: Option<HaplogroupResult> = None;
        for child in &node.children {
            let child_entry = map.get(&child.name);
            let (cd, ca, ct) = child_entry
                .map(|r| {
                    (
                        r.matching_snps as i32,
                        r.ancestral_matches as i32,
                        r.total_snps as i32,
                    )
                })
                .unwrap_or((0, 0, 0));
            // Revised rule: allow descent when there is sufficient non-trivial derived support and not worse than ancestral.
            // This helps nodes with mixed evidence (e.g., 3 derived / 3 ancestral) proceed over siblings defined by INDELs.
            let clear_derived_here = (cd >= 2 && cd >= ca) || (cd >= 1 && ca == 0);
            // Hard gate: do not consider children that have zero SNP loci (e.g., INDEL-only definitions)
            let has_snp_loci = ct > 0;
            // Determine candidate leaves limited to ONE extra level below this child (grandchildren that are leaves), or the child if it is a leaf
            let mut best_leaf: Option<HaplogroupResult> = None;
            let mut any_downstream_derived = false;
            // Include the child itself if it is a leaf and has SNP loci
            if child.children.is_empty() && has_snp_loci {
                if let Some(entry) = map.get(&child.name) {
                    any_downstream_derived = any_downstream_derived || entry.matching_snps > 0;
                    best_leaf = Some(entry.clone());
                }
            }
            // Consider direct grandchildren that are leaves
            for gc in &child.children {
                if gc.children.is_empty() {
                    if let Some(entry) = map.get(&gc.name) {
                        if entry.matching_snps > 0 {
                            any_downstream_derived = true;
                        }
                        match &best_leaf {
                            None => best_leaf = Some(entry.clone()),
                            Some(curr) => {
                                if cmp_leaf(entry, curr).is_lt() {
                                    best_leaf = Some(entry.clone());
                                }
                            }
                        }
                    }
                }
            }
            // Parent gate: if parent is strictly all-ancestral, do not descend
            let can_descend_parent_gate = !(parent_all_ancestral);
            // Decide if we can descend into this child
            let can_descend = can_descend_parent_gate
                && (clear_derived_here || (single_ancestral && any_downstream_derived));
            if can_descend {
                if let Some(bl) = best_leaf {
                    match (&best_child_leaf, best_child) {
                        (None, _) => {
                            best_child_leaf = Some(bl);
                            best_child = Some(child);
                        }
                        (Some(curr_leaf), Some(_)) => {
                            if cmp_leaf(&bl, curr_leaf).is_lt() {
                                best_child_leaf = Some(bl);
                                best_child = Some(child);
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
        // Only allow ONE extra level of descent from the current node
        if let Some(ch) = best_child {
            return ch;
        }
        // Otherwise, stop here per rules
        node
    }

    let mut best_leaf_overall: Option<HaplogroupResult> = None;
    // Build a leaf based on the rule-based picker
    let chosen = pick_best_by_rules(tree, &unique_results, &remaining, tree, 0);
    if chosen.children.is_empty() {
        // Chosen is already a leaf
        if let Some(entry) = unique_results.get(&chosen.name).cloned() {
            best_leaf_overall = Some(entry);
        }
    } else {
        // Choose the best descendant LEAF under the chosen subtree using the same comparator
        let mut best_under_chosen: Option<HaplogroupResult> = None;
        for leaf in &remaining {
            if let Some(path) = tree::find_path_to_root(tree, &leaf.name, &[]) {
                if path.contains(&chosen.name) {
                    match &best_under_chosen {
                        None => best_under_chosen = Some(leaf.clone()),
                        Some(curr) => {
                            if cmp(leaf, curr).is_lt() {
                                best_under_chosen = Some(leaf.clone());
                            }
                        }
                    }
                }
            }
        }
        if let Some(best) = best_under_chosen {
            best_leaf_overall = Some(best);
        } else {
            // Fallback: From the root, choose the child subtree whose best leaf satisfies the comparator
            for child in &tree.children {
                // find best leaf under this child
                let mut best_in_child: Option<HaplogroupResult> = None;
                for leaf in &remaining {
                    if let Some(path) = tree::find_path_to_root(tree, &leaf.name, &[]) {
                        if path.contains(&child.name) {
                            match &best_in_child {
                                None => best_in_child = Some(leaf.clone()),
                                Some(curr) => {
                                    if cmp(leaf, curr).is_lt() {
                                        best_in_child = Some(leaf.clone());
                                    }
                                }
                            }
                        }
                    }
                }
                if let Some(candidate) = best_in_child {
                    match &best_leaf_overall {
                        None => best_leaf_overall = Some(candidate),
                        Some(curr) => {
                            if cmp(&candidate, curr).is_lt() {
                                best_leaf_overall = Some(candidate);
                            }
                        }
                    }
                }
            }
        }
    }

    // Optional debug on key splits (IJK -> IJ vs K, K -> R)
    if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
        // Helper: find first node matching a name prefix or exact known markers
        fn find_by_pred<'a, F: Fn(&str) -> bool>(
            node: &'a Haplogroup,
            pred: &F,
        ) -> Option<&'a Haplogroup> {
            if pred(&node.name) {
                return Some(node);
            }
            for c in &node.children {
                if let Some(found) = find_by_pred(c, pred) {
                    return Some(found);
                }
            }
            None
        }
        // Helper: pretty print tallies and decisions
        let mut print_children = |label: &str, parent: &Haplogroup| {
            eprintln!(
                "[debug] {} split at {}: {} children",
                label,
                parent.name,
                parent.children.len()
            );
            // Determine which child would be chosen by rules from this parent
            let chosen = pick_best_by_rules(parent, &unique_results, &remaining, tree, 0);
            for ch in &parent.children {
                if let Some(entry) = unique_results.get(&ch.name) {
                    let total = entry.total_snps;
                    let mut derived = entry.matching_snps;
                    let ancestral = entry.ancestral_matches;
                    let no_calls = entry.no_calls;
                    let norm_single = total > 5 && derived == 1;
                    if norm_single {
                        derived = 0;
                    }
                    let conflict = total > 4 && (ancestral as u32) >= ((total + 1) / 2);
                    let clear_here = ((entry.matching_snps >= 1 && entry.ancestral_matches == 0)
                        || (entry.matching_snps >= 2
                            && entry.matching_snps > entry.ancestral_matches));
                    let mark = if ch.name == chosen.name {
                        " <chosen>"
                    } else {
                        ""
                    };
                    eprintln!(
                        "[debug]   child {}{}\tbranch_score={:.4}\t(der={}, anc={}, noc={}, total={})\tnormalized_single={}\tconflict={}\tclear_here={}",
                        entry.name,
                        mark,
                        entry.score,
                        derived,
                        ancestral,
                        no_calls,
                        total,
                        norm_single,
                        conflict,
                        clear_here
                    );
                } else {
                    eprintln!("[debug]   child {} has no score entry", ch.name);
                }
            }
            eprintln!("[debug] {} decision -> {}", label, chosen.name);
        };
        // IJK split
        if let Some(ijk) = find_by_pred(tree, &|n: &str| n.starts_with("IJK-") || n.contains("IJK"))
        {
            print_children("IJK", ijk);
        }
        // K split
        if let Some(k_node) = find_by_pred(tree, &|n: &str| {
            n.starts_with("K-") || n == "K" || n.contains("K-M9")
        }) {
            print_children("K", k_node);
        }
    }

    // Take the top LEAF from the best root child path and output the LEAF first, then its ancestral path
    if let Some(top_result) = best_leaf_overall.or_else(|| remaining.first().cloned()) {
        // Output the top leaf first
        if let Some(entry) = unique_results.get(&top_result.name).cloned() {
            if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
                eprintln!(
                    "[debug] PATH {}\tbranch={:.4}\tcum={:.4}",
                    entry.name, entry.score, entry.score
                );
            }
            ordered_results.push(entry);
        }
        // Then include its ancestral path (excluding the leaf itself), from root down to parent
        if let Some(path) = tree::find_path_to_root(tree, &top_result.name, &[]) {
            let mut cumulative = 0.0;
            for name in path.into_iter().rev().filter(|n| n != &top_result.name) {
                if let Some(entry) = unique_results.get(&name).cloned() {
                    cumulative += entry.score;
                    if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
                        eprintln!(
                            "[debug] PATH {}\tbranch={:.4}\tcum={:.4}",
                            entry.name, entry.score, cumulative
                        );
                    }
                    ordered_results.push(entry);
                }
            }
        }
        // Remove the top leaf from remaining
        if let Some(pos) = remaining.iter().position(|r| r.name == top_result.name) {
            remaining.remove(pos);
        }
    }

    // Optionally include off-path leaves for diagnostics with evidence gating
    if include_off_path {
        let eps = 1e-9f64;
        // Evidence gate: (matches >= 2) || (matches >= 1 && ancestral == 0), and score > EPS
        let mut gated: Vec<HaplogroupResult> = remaining
            .into_iter()
            .filter(|r| {
                let matches = r.matching_snps as i32;
                let ancestral = r.ancestral_matches as i32;
                (matches >= 2 || (matches >= 1 && ancestral == 0)) && r.score > eps
            })
            .collect();
        // Order by score desc, then cumulative SNPs desc
        gated.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
        });
        ordered_results.extend(gated);
    }

    ordered_results
}

// Helper: compute branch-local metrics for a single haplogroup name for reporting purposes
// Returns (matching_derived, mismatching_conflicts, ancestral_matches, no_calls, total_snps)
fn compute_branch_local_metrics(
    haplogroup_name: &str,
    tree: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    build_id: &str,
) -> (u32, u32, u32, u32, u32) {
    let mut matching: u32 = 0;
    let mut mismatching: u32 = 0;
    let mut ancestral: u32 = 0;
    let mut no_calls: u32 = 0;
    let mut total: u32 = 0;

    if let Some(node) = find_haplogroup(tree, haplogroup_name) {
        for locus in &node.loci {
            if matches!(locus.loci_type, LociType::SNP) {
                if let Some(coord) = locus.coordinates.get(build_id) {
                    total += 1;
                    match snp_calls.get(&coord.position) {
                        Some((called_base, depth_called, _freq)) => {
                            if *depth_called < 1 {
                                no_calls += 1;
                                continue;
                            }
                            let base = called_base.to_ascii_uppercase();
                            let der = coord
                                .derived
                                .chars()
                                .next()
                                .unwrap_or('N')
                                .to_ascii_uppercase();
                            let anc = coord
                                .ancestral
                                .chars()
                                .next()
                                .unwrap_or('N')
                                .to_ascii_uppercase();
                            if base == der {
                                matching += 1;
                            } else if base == anc {
                                ancestral += 1;
                            } else {
                                mismatching += 1;
                            }
                        }
                        None => {
                            no_calls += 1;
                        }
                    }
                }
            }
        }
    }

    (matching, mismatching, ancestral, no_calls, total)
}
