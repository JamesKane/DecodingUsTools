use crate::haplogroup::types::{Haplogroup, HaplogroupResult, HaplogroupScore, LociType};
use std::collections::{HashMap, HashSet};

const HIGH_QUALITY_THRESHOLD: f64 = 0.7;
const MEDIUM_QUALITY_THRESHOLD: f64 = 0.5;
// Respect caller-level min depth by keeping scoring permissive; caller already filters by depth.
const MIN_DEPTH: u32 = 1;

// --- Debug helpers ---
// Diagnostics for INDEL-only branches (e.g., DF13 -> Z39589/FGC11134 split):
// - Set DECODINGUS_DEBUG_SCORES=1 to enable debug output.
// - Optionally set DECODINGUS_DEBUG_NODE to a substring of the haplogroup name to narrow logs
//   (default if unset is to match "Z39589").
// - Optionally set DECODINGUS_DEBUG_CHILD_LOCI=1 to print per-locus details for the
//   immediate children of an INDEL-only node.
//
// Acceptance criteria: when a node has 0 SNP loci (bm.total == 0) and matches the filter,
// we print a header and one-line summaries for each immediate child including metrics and
// raw/adjusted scores, without affecting scoring behavior.
fn debug_scores_enabled() -> bool {
    std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1")
}

fn matches_debug_node(name: &str) -> bool {
    match std::env::var("DECODINGUS_DEBUG_NODE").ok() {
        Some(filter) if !filter.is_empty() => name.contains(&filter),
        _ => false, // no default hard-coded target
    }
}

fn child_loci_debug_enabled() -> bool {
    std::env::var("DECODINGUS_DEBUG_CHILD_LOCI").ok().as_deref() == Some("1")
}

#[derive(Clone, Debug, Default)]
struct BranchMetrics {
    derived: u32,
    ancestral: u32,
    conflicts: u32, // true mismatches (neither ancestral nor derived)
    no_calls: u32,
    low_quality: u32,
    total: u32,
    raw_branch_score: f64,
}

/// Computes the branch metrics for a given haplogroup based on observed SNP (Single Nucleotide Polymorphism)
/// calls, considering read depth, allele frequency, and ancestral/derived states.
///
/// # Parameters
///
/// * `haplogroup` - A reference to a `Haplogroup` object that defines loci (positions) of interest and
///                  contains associated metadata like locus type and coordinate mappings.
/// * `snp_calls` - A `HashMap` mapping genomic positions (u32) to observed SNP calls, represented as a tuple
///                 containing the called base (`char`), the read depth (`u32`), and the allele frequency (`f64`).
/// * `build_id` - A string slice representing the genomic build ID, used to match coordinate mappings in the
///                haplogroup.
/// * `depth` - A `u32` indicating the sequencing depth for the current context, used to adjust computed
///             branch scores.
/// * `debug_per_locus` - A boolean flag enabling per-locus debug output for diagnostic purposes, providing
///                       detailed information about SNP evaluation at each position.
///
/// # Returns
///
/// * `BranchMetrics` - A struct summarizing SNP-based metrics for the haplogroup, including counts of defining
///                     loci, tallied allele classifications (e.g., derived, ancestral, conflicts), as well as
///                     computed confidence scores (raw_branch_score).
///
/// # Process
///
/// The function follows these steps:
/// 1. Filters loci in the haplogroup to consider only SNP loci with coordinates specified for the given `build_id`.
/// 2. Initializes a default `BranchMetrics` object to keep track of metrics.
/// 3. Iterates over the filtered loci and analyzes associated SNP calls from the `snp_calls` map. Depending on
///    SNP call availability and quality thresholds:
///     - Tallies counts for derived, ancestral, conflict, or no-call loci.
///     - Records loci as low-quality if allele frequency falls below `HIGH_QUALITY_THRESHOLD`.
/// 4. When debugging is enabled, emits detailed per-locus diagnostic output.
/// 5. Computes an overall raw branch score based on derived allele counts, their frequency, read depth, and
///    certain thresholds and scaling factors.
/// 6. Returns the populated `BranchMetrics` object.
///
/// # Note
///
/// - Debug output will appear in stderr and is controlled through `debug_per_locus` or specific environment
///   variables (`DECODINGUS_DEBUG_SCORES` and `DECODINGUS_DEBUG_NODE`).
/// - The function assumes a predefined constant `MIN_DEPTH` for minimum read depth and `HIGH_QUALITY_THRESHOLD`
///   for allele frequency to consider high confidence.
/// - The scoring logic penalizes branches with fewer callable loci or an overwhelming presence of ancestral
///   alleles, compared to derived ones.
///
/// # Example
///
/// ```
/// # use std::collections::HashMap;
/// # let haplogroup = /* construct Haplogroup */;
/// # let snp_calls = HashMap::new(); // Mapping of positions to SNP calls
/// # let build_id = "GRCh38";
/// # let depth = 100;
/// # let debug_per_locus = false;
///
/// let metrics = compute_branch_metrics(&haplogroup, &snp_calls, build_id, depth, debug_per_locus);
///
/// println!("Derived loci: {}", metrics.derived);
/// println!("Raw branch score: {:.2}", metrics.raw_branch_score);
/// ```
fn compute_branch_metrics(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    build_id: &str,
    depth: u32,
    debug_per_locus: bool,
) -> BranchMetrics {
    let defining_loci: Vec<_> = haplogroup
        .loci
        .iter()
        .filter(|locus| {
            matches!(locus.loci_type, LociType::SNP) && locus.coordinates.contains_key(build_id)
        })
        .collect();

    let mut bm = BranchMetrics::default();

    // Debug gating:
    // - When debug_per_locus=true, always emit per-locus lines (caller already gated by env vars and node match)
    // - Otherwise, emit only when DECODINGUS_DEBUG_SCORES=1 and the node matches DECODINGUS_DEBUG_NODE filter
    let debug_scores = debug_scores_enabled();
    let debug = debug_per_locus || (debug_scores && matches_debug_node(&haplogroup.name));

    bm.total = defining_loci.len() as u32;

    // Iterate loci and tally metrics from snp_calls
    for locus in defining_loci {
        if let Some(coord) = locus.coordinates.get(build_id) {
            let pos = coord.position;
            match snp_calls.get(&pos) {
                None => {
                    bm.no_calls += 1;
                    if debug {
                        if debug_per_locus {
                            eprintln!("[debug.loc] {:<8} pos={} no_call", locus.name, pos);
                        }
                    }
                }
                Some((call_base, depth_called, freq)) => {
                    if *depth_called < MIN_DEPTH {
                        bm.no_calls += 1;
                        if debug && debug_per_locus {
                            eprintln!(
                                "[debug.loc] {:<8} pos={} depth={} < MIN_DEPTH no_call",
                                locus.name, pos, depth_called
                            );
                        }
                        continue;
                    }

                    // Do not early-exit on medium frequency; classify allele first and mark low quality afterwards

                    // Compare against expected alleles
                    let der = coord
                        .derived
                        .as_bytes()
                        .get(0)
                        .map(|b| (*b as char).to_ascii_uppercase());
                    let anc = coord
                        .ancestral
                        .as_bytes()
                        .get(0)
                        .map(|b| (*b as char).to_ascii_uppercase());
                    let base = call_base.to_ascii_uppercase();

                    if Some(base) == der {
                        bm.derived += 1;
                        // downgrade to low_quality bucket if not high confidence
                        if *freq < HIGH_QUALITY_THRESHOLD {
                            bm.low_quality += 1;
                        }
                        if debug && debug_per_locus {
                            eprintln!(
                                "[debug.loc] {:<8} pos={} base={} depth={} freq={:.3} DER",
                                locus.name, pos, call_base, depth_called, freq
                            );
                        }
                    } else if Some(base) == anc {
                        bm.ancestral += 1;
                        if *freq < HIGH_QUALITY_THRESHOLD {
                            bm.low_quality += 1;
                        }
                        if debug && debug_per_locus {
                            eprintln!(
                                "[debug.loc] {:<8} pos={} base={} depth={} freq={:.3} ANC",
                                locus.name, pos, call_base, depth_called, freq
                            );
                        }
                    } else {
                        bm.conflicts += 1;
                        // mark as low quality if not high confidence
                        if *freq < HIGH_QUALITY_THRESHOLD {
                            bm.low_quality += 1;
                        }
                        if debug && debug_per_locus {
                            eprintln!(
                                "[debug.loc] {:<8} pos={} base={} depth={} freq={:.3} CONFLICT anc='{}' der='{}'",
                                locus.name,
                                pos,
                                call_base,
                                depth_called,
                                freq,
                                coord.ancestral,
                                coord.derived
                            );
                        }
                    }
                }
            }
        }
    }

    // Base score (pre-experimental multipliers). Keep aligned with legacy branch calc below.
    let callable = bm.derived + bm.ancestral + bm.low_quality;
    if callable > 0 {
        let ratio = (bm.derived as f64) / (callable as f64);
        let evidence = (1.0 + bm.derived as f64).ln();
        let mut branch_score = ratio * evidence;
        let depth_factor = 1.0 + 0.05 * (depth as f64);
        branch_score *= depth_factor;
        if bm.derived > 0 {
            if bm.ancestral >= bm.derived {
                branch_score *= 0.3;
            }
            if bm.ancestral > bm.derived * 2 {
                branch_score *= 0.15;
            }
        }
        if bm.total > 0 {
            let callable_frac = callable as f64 / bm.total as f64;
            if callable_frac < 0.3 {
                branch_score *= 0.6;
            } else if callable_frac < 0.6 {
                branch_score *= 0.85;
            }
        }
        if bm.derived < 2 {
            branch_score *= 0.9;
        }
        let quality_factor = if bm.low_quality == 0 { 1.1 } else { 0.92 };
        bm.raw_branch_score = branch_score * quality_factor;
    } else {
        bm.raw_branch_score = 0.0;
    }

    bm
}

/**
 * Computes the score for a given haplogroup branch based on its associated SNP calls, parent information,
 * and child branches, storing the results in a scores vector for further evaluation.
 *
 * This method processes the haplogroup tree recursively, utilizing branch metrics and configurable weight factors
 * to compute adjusted scores for each haplogroup node. The computation considers derived SNPs, ancestral SNPs,
 * low-quality data, and other relevant factors to score and evaluate the haplogroup.
 *
 * # Arguments
 *
 * * `haplogroup` - The current haplogroup node being evaluated. Contains information about the name, loci, and child nodes.
 * * `snp_calls` - A map containing SNP calls for various loci. Each SNP call is represented by a tuple of `(base, position, quality)`.
 * * `scores` - A mutable vector to store the resulting scores for all haplogroups evaluated during the recursion.
 * * `parent_info` - Optional information about the parent score and cumulative SNPs for the parent haplogroup node.
 * * `depth` - The depth of the current haplogroup node in the haplogroup tree structure.
 * * `build_id` - A string identifier representing the genome build used for coordinate system (e.g., "hg19").
 *
 * # Returns
 *
 * A tuple containing:
 * * `HaplogroupScore` - The computed score for the current haplogroup node, encapsulating matches, mismatches, and other metrics.
 * * `HashSet<u32>` - A cumulative set of SNP positions that combine the current node's loci and parent cumulative SNPs.
 *
 * # Algorithm Details
 *
 * 1. **Score Adjustments**:
 *    - The raw branch score gets adjusted based on a set of conditions, such as:
 *      - Shallow node constraints (e.g., penalizing singleton derived calls or excessive ancestral calls).
 *      - Derived vs ancestral contrast.
 *      - Callable SNP fraction scaling.
 *
 * 2. **Sibling Contrast and Softmax**:
 *    - Contrast metrics for child haplogroups are computed using sibling bonus adjustments. Child nodes are evaluated
 *      with softmax-derived scaling.
 *
 * 3. **Pruning Heuristics**:
 *    - Prevents descending into child nodes that show no derived SNP evidence or whose immediate children also lack derived evidence.
 *    - Applies this pruning rule to reduce propagation into haplogroups that do not align statistically.
 *
 * 4. **Debugging Support**:
 *    - Prints detailed information about the current node at key splits if debugging is enabled via the `DECODINGUS_DEBUG_SCORES` environment variable.
 *    - Includes options to show child relationships and immediate metrics (`DECODINGUS_DEBUG_TREE`).
 *
 * # Panics
 *
 * Panics if an invalid haplogroup tree state is encountered where a child node shares the same name as its parent node.
 *
 * # Example Usage
 *
 * ```rust
 * let mut scores = Vec::new();
 * let haplogroup = Haplogroup::new("H", ...); // Simplified example
 * let snp_calls = HashMap::new();
 * let result = calculate_branch_score(&haplogroup, &snp_calls, &mut scores, None, 0, "hg19");
 *
 * println!("Result score: {:?}", result.0.score);
 * println!("Cumulative SNPs: {:?}", result.1);
 * ```
 *
 * # Notes
 *
 * - The function aims to balance sensitivity and specificity when evaluating haplogroups by applying configurable
 *   weights and thresholds for scoring derived vs ancestral SNPs.
 * - Depth-based adjustments help account for the hierarchical nature of haplogroup trees.
 */
fn calculate_branch_score(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<HaplogroupResult>,
    parent_info: Option<&(HaplogroupScore, HashSet<u32>)>,
    depth: u32,
    build_id: &str,
) -> (HaplogroupScore, HashSet<u32>) {
    // Configurable multipliers (conservative defaults)
    let dampen_majority_anc = 0.35f64;
    let singleton_downweight = 0.6f64;
    let alpha: f64 = 1.7; // ancestral weight in contrast
    let bonus_budget: f64 = 0.3; // sibling bonus budget
    let coherence_bonus: f64 = 1.05; // path-coherence

    // Prepare cumulative SNP set
    let mut cumulative_snps: HashSet<u32> = parent_info
        .map(|(_, snps)| snps.clone())
        .unwrap_or_default();
    // Add current node loci positions
    for locus in &haplogroup.loci {
        if matches!(locus.loci_type, LociType::SNP) {
            if let Some(coord) = locus.coordinates.get(build_id) {
                cumulative_snps.insert(coord.position);
            }
        }
    }

    // Compute branch metrics and adjusted score for current node
    let bm = compute_branch_metrics(haplogroup, snp_calls, build_id, depth, false);
    let mut adjusted_branch = bm.raw_branch_score;
    // Stronger ancestral penalty and singleton guard at shallow nodes
    if bm.total > 4 && bm.ancestral as u32 >= ((bm.total + 1) / 2) {
        adjusted_branch *= dampen_majority_anc;
    }
    if bm.total > 5 && bm.derived == 1 {
        adjusted_branch *= singleton_downweight;
    }
    // Minimum support guard for shallow nodes: near-zero unless (derived >=2) or (derived>=1 and ancestral==0)
    let shallow = depth <= 6;
    if shallow {
        let ok_support = (bm.derived >= 2) || (bm.derived >= 1 && bm.ancestral == 0);
        if !ok_support {
            adjusted_branch *= 0.1;
        }
    }
    // Strengthen ancestral vs derived contrast
    if bm.derived > 0 {
        if bm.ancestral >= bm.derived {
            adjusted_branch *= 0.25;
        }
        if bm.ancestral > bm.derived * 2 {
            adjusted_branch *= 0.10;
        }
    }
    // Callable fraction scaling a bit stronger
    if bm.total > 0 {
        let callable = bm.derived + bm.ancestral + bm.low_quality;
        let callable_frac = callable as f64 / bm.total as f64;
        if callable_frac < 0.3 {
            adjusted_branch *= 0.5;
        } else if callable_frac < 0.6 {
            adjusted_branch *= 0.8;
        }
    }

    let mut current_score = HaplogroupScore::default();
    current_score.matches = bm.derived as usize;
    current_score.ancestral_matches = bm.ancestral as usize;
    current_score.no_calls = bm.no_calls as usize;
    current_score.total_snps = bm.total as usize;
    current_score.score = adjusted_branch;

    // Push branch-only result for the current node
    scores.push(HaplogroupResult {
        name: haplogroup.name.clone(),
        score: current_score.score,
        matching_snps: bm.derived,
        mismatching_snps: bm.conflicts,
        low_quality_snps: bm.low_quality,
        ancestral_matches: bm.ancestral,
        no_calls: bm.no_calls,
        total_snps: bm.total,
        cumulative_snps: cumulative_snps.len() as u32,
        depth,
    });

    if haplogroup.children.is_empty() {
        return (current_score, cumulative_snps);
    }

    // Precompute child metrics to derive sibling bonuses
    let mut child_metrics: Vec<(&Haplogroup, BranchMetrics)> = Vec::new();
    for child in &haplogroup.children {
        if child.name == haplogroup.name {
            panic!(
                "Invalid haplogroup tree state: Child node '{}' has same name as its parent",
                child.name
            );
        }
        let cbm = compute_branch_metrics(child, snp_calls, build_id, depth + 1, false);
        child_metrics.push((child, cbm));
    }

    let parent_all_ancestral = bm.total > 0 && bm.derived == 0 && bm.ancestral == bm.total;

    // Sibling contrast -> softmax
    let mut contrasts: Vec<f64> = child_metrics
        .iter()
        .map(|(_, cbm)| (cbm.derived as f64) - alpha * (cbm.ancestral as f64))
        .collect();
    let max_c = contrasts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exp_vals: Vec<f64> = contrasts.iter().map(|c| (*c - max_c).exp()).collect();
    let sum_exp: f64 = exp_vals.iter().sum::<f64>().max(1e-9);

    // Optional debug at key splits
    let debug = std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1");
    let debug_node = debug && matches_debug_node(&haplogroup.name);
    if debug_node {
        eprintln!("[debug] Split at {} (depth {}):", haplogroup.name, depth);
        eprintln!(
            "[debug]   parent bm: der={} anc={} noc={} lowq={} total={} raw={:.4} adj={:.4}",
            bm.derived,
            bm.ancestral,
            bm.no_calls,
            bm.low_quality,
            bm.total,
            bm.raw_branch_score,
            adjusted_branch
        );

        // Optional tree view of immediate children names
        if std::env::var("DECODINGUS_DEBUG_TREE").ok().as_deref() == Some("1") {
            let child_names: Vec<String> =
                haplogroup.children.iter().map(|c| c.name.clone()).collect();
            eprintln!(
                "[debug.tree] immediate children: [{}]",
                child_names.join(", ")
            );
        }

        // Immediate child branch-local metrics (pre-recursion)
        for (child, cbm) in child_metrics.iter() {
            let mut child_adj = cbm.raw_branch_score;
            if cbm.total > 4 && cbm.ancestral as u32 >= ((cbm.total + 1) / 2) {
                child_adj *= dampen_majority_anc;
            }
            if cbm.total > 5 && cbm.derived == 1 {
                child_adj *= singleton_downweight;
            }
            eprintln!(
                "[debug]   immediate child {:<20} der={} anc={} noc={} lowq={} total={} raw={:.4} adj={:.4}",
                child.name,
                cbm.derived,
                cbm.ancestral,
                cbm.no_calls,
                cbm.low_quality,
                cbm.total,
                cbm.raw_branch_score,
                child_adj
            );
        }
    }

    for (idx, (child, cbm)) in child_metrics.iter().enumerate() {
        // Adjust child branch score locally
        let mut child_adj = cbm.raw_branch_score;
        if cbm.total > 4 && cbm.ancestral as u32 >= ((cbm.total + 1) / 2) {
            child_adj *= dampen_majority_anc;
        }
        if cbm.total > 5 && cbm.derived == 1 {
            child_adj *= singleton_downweight;
        }
        let soft = exp_vals[idx] / sum_exp;
        let sibling_bonus = 1.0 + bonus_budget * soft;
        let clear_local = (cbm.derived > cbm.ancestral && cbm.derived >= 2)
            || (cbm.derived >= 1 && cbm.ancestral == 0);
        let coh = if clear_local && !parent_all_ancestral {
            coherence_bonus
        } else {
            1.0
        };

        // Pruning rule: do not descend into a child that shows no derived evidence,
        // and whose immediate children (if any) also show no derived evidence while
        // showing ancestral calls. This prevents walking down the wrong sibling when
        // another sibling (e.g., FGC11134) has the only derived support.
        let mut prune_child = false;
        if cbm.total > 0 {
            // Child has SNP loci but none derived and at least one ancestral → prune
            if cbm.derived == 0 && cbm.ancestral > 0 {
                prune_child = true;
            }
        } else {
            // INDEL-only child: look at immediate grandchildren metrics
            if !child.children.is_empty() {
                let mut any_gc_with_loci = false;
                let mut all_gc_no_derived = true;
                let mut any_gc_ancestral = false;
                for gc in &child.children {
                    let gcbm = compute_branch_metrics(gc, snp_calls, build_id, depth + 2, false);
                    if gcbm.total > 0 {
                        any_gc_with_loci = true;
                        if gcbm.derived > 0 {
                            all_gc_no_derived = false;
                        }
                        if gcbm.ancestral > 0 {
                            any_gc_ancestral = true;
                        }
                    }
                }
                if any_gc_with_loci && all_gc_no_derived && any_gc_ancestral {
                    prune_child = true;
                }
            }
        }
        if prune_child {
            if debug_node {
                eprintln!(
                    "[debug]   pruning descent into {:<20} (der={} anc={} total={})",
                    child.name, cbm.derived, cbm.ancestral, cbm.total
                );
            }
            continue;
        }

        // Recurse down the child to get its cumulative score
        let (child_score, child_cumulative) = calculate_branch_score(
            child,
            snp_calls,
            scores,
            Some(&(current_score.clone(), cumulative_snps.clone())),
            depth + 1,
            build_id,
        );

        // Combine scores (apply sibling/path multipliers to the whole child subtree to influence ordering)
        let combined_matches = current_score.matches + child_score.matches;
        let combined_ancestral = current_score.ancestral_matches + child_score.ancestral_matches;
        let combined_no_calls = current_score.no_calls + child_score.no_calls;
        let combined_total = current_score.total_snps + child_score.total_snps;
        let combined_score = current_score.score + (child_score.score * sibling_bonus * coh);

        if debug_node {
            eprintln!(
                "[debug]   subtree via {:<20} final={:.4} matches={} total_snps={} depth={}",
                child.name,
                child_score.score * sibling_bonus * coh,
                child_score.matches,
                child_score.total_snps,
                depth + 1
            );
        }

        // Push cumulative result for this child (use child-side metrics for conflicts/lowq)
        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: combined_score,
            matching_snps: combined_matches as u32,
            mismatching_snps: cbm.conflicts,
            low_quality_snps: cbm.low_quality,
            ancestral_matches: combined_ancestral as u32,
            no_calls: combined_no_calls as u32,
            total_snps: combined_total as u32,
            cumulative_snps: child_cumulative.len() as u32,
            depth: depth + 1,
        });
    }

    (current_score, cumulative_snps)
}

pub(crate) fn calculate_haplogroup_score(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<HaplogroupResult>,
    parent_info: Option<&(HaplogroupScore, HashSet<u32>)>,
    depth: u32,
    build_id: &str,
) -> (HaplogroupScore, HashSet<u32>) {
    // Experimental scoring is now the default implementation
    calculate_branch_score(haplogroup, snp_calls, scores, parent_info, depth, build_id)
}
