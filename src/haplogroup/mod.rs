pub mod caller;
mod scoring;
mod tree;
pub(crate) mod types;
mod validation;

use crate::haplogroup::types::{Haplogroup, HaplogroupResult, Locus};
use crate::utils::cache::TreeType;
use crate::utils::liftover::Liftover;
use crate::types::ReferenceGenome;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::{self, Read as BamRead, Reader};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn analyze_haplogroup(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    tree_type: TreeType,
    provider: crate::cli::TreeProvider,
    show_snps: bool,
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
    let projected_tree = tree::project_tree_to_build(&haplogroup_tree, src_build_id, bam_build.clone(), &lifter);

    // Debug: verify projection worked for key haplogroups
    let debug_projection = std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1");
    if debug_projection {
        eprintln!("[debug] Tree projection: {} -> {}", src_build_id, bam_build.name());

        // Check a few key haplogroups to verify their coordinates changed
        let check_names = vec!["R-P312", "R-L21", "R-M269"];
        for check_name in check_names {
            if let Some(original_hg) = find_haplogroup(&haplogroup_tree, check_name) {
                if let Some(projected_hg) = find_haplogroup(&projected_tree, check_name) {
                    eprintln!("[debug] Haplogroup {} projection check:", check_name);
                    for locus in &original_hg.loci {
                        if let Some(orig_coord) = locus.coordinates.get(src_build_id) {
                            eprintln!("[debug]   Original ({}): {} pos={} anc={} der={}",
                                      src_build_id, locus.name, orig_coord.position,
                                      orig_coord.ancestral, orig_coord.derived);
                        }
                    }
                    for locus in &projected_hg.loci {
                        if let Some(proj_coord) = locus.coordinates.get(bam_build.name()) {
                            eprintln!("[debug]   Projected ({}): {} pos={} anc={} der={}",
                                      bam_build.name(), locus.name, proj_coord.position,
                                      proj_coord.ancestral, proj_coord.derived);
                        } else {
                            eprintln!("[debug]   Projected ({}): {} MISSING (liftover failed)",
                                      bam_build.name(), locus.name);
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
        eprintln!("[debug] Collected {} unique positions from projected tree for build '{}'",
                  positions_bam.len(), build_id);
        // Show a sample of collected positions
        let mut sample_positions: Vec<u32> = positions_bam.keys().copied().collect();
        sample_positions.sort();
        eprintln!("[debug] Sample positions collected: {:?}",
                  sample_positions.iter().take(10).collect::<Vec<_>>());
    }

    // Optional: focus debug on a specific site by locus name (e.g., M526)
    let debug_site = std::env::var("DECODINGUS_DEBUG_SITE").ok();
    let mut debug_tree_sites: Vec<u32> = Vec::new();
    if let Some(ref site_tag) = debug_site {
        for (pos, entries) in positions_bam.iter() {
            if entries.iter().any(|(_, locus)| locus.name.contains(site_tag)) {
                debug_tree_sites.push(*pos);
            }
        }
        if !debug_tree_sites.is_empty() {
            eprintln!("[debug] Tracking sites by locus tag '{}': {:?}", site_tag, debug_tree_sites);
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

    // Debug: dump snp_calls for P312 position if debugging
    if debug_projection || debug_site.is_some() {
        let p312_pos = 20_901_962u32;
        if let Some((base, depth, freq)) = snp_calls.get(&p312_pos) {
            eprintln!("[debug] snp_calls at P312 position {}: base='{}' depth={} freq={:.3}",
                      p312_pos, base, depth, freq);
        } else {
            eprintln!("[debug] snp_calls at P312 position {}: NO CALL FOUND", p312_pos);
        }

        // Also check if positions_bam has P312
        if let Some(entries) = positions_bam.get(&p312_pos) {
            eprintln!("[debug] positions_bam at P312 position {} has {} entries:", p312_pos, entries.len());
            for (hg_name, locus) in entries {
                eprintln!("[debug]   - haplogroup='{}' locus='{}'", hg_name, locus.name);
            }
        } else {
            eprintln!("[debug] positions_bam at P312 position {}: NOT TRACKED", p312_pos);
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

    let mut ordered_scores = collect_scored_paths(scores, &projected_tree);
    // Post-process: filter zero-score artifacts and sort by depth descending (leaf to root)
    ordered_scores.retain(|r| r.score > 0.0);
    ordered_scores.sort_by(|a, b| b.depth.cmp(&a.depth).then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal)));

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
        // Focused view for R-* lineage entries
        let r_entries: Vec<_> = ordered_scores
            .iter()
            .filter(|r| r.name.starts_with("R-"))
            .take(25)
            .collect();
        if !r_entries.is_empty() {
            eprintln!("[debug] R-lineage snapshot:");
            for r in r_entries {
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
        if show_snps {
            let (matching_snps, mismatching_snps, no_call_snps) =
                get_snp_details(&result.name, &projected_tree, &snp_calls, build_id);
            writeln!(
                writer,
                "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.name,
                result.score,
                result.matching_snps,
                result.mismatching_snps,
                result.ancestral_matches,
                result.no_calls,
                result.total_snps,
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
                result.matching_snps,
                result.mismatching_snps,
                result.ancestral_matches,
                result.no_calls,
                result.total_snps,
                result.cumulative_snps,
                result.depth
            )?;
        }
    }

    progress.finish_with_message("Analysis complete!");
    Ok(())
}

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
        no_calls.join(";")
    )
}

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

fn collect_scored_paths(scores: Vec<HaplogroupResult>, tree: &Haplogroup) -> Vec<HaplogroupResult> {
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
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
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
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    };

    // Rule-based descent to pick the best leaf closest to a terminal node per guidelines
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
        if entry.is_none() { return node; }
        let e = entry.unwrap();
        let total = e.total_snps;
        let mut derived = e.matching_snps;
        let ancestral = e.ancestral_matches;
        // When there is only one derived call in a block of more than five, treat the branch as ancestral
        if total > 5 && derived == 1 { derived = 0; }
        // Define conflict for succession tracking only for larger blocks (>4) with half-or-more ancestral
        let block_conflict = total > 4 && (ancestral as u32) >= ((total + 1) / 2);
        let consec = if block_conflict { consec_conflicts + 1 } else { 0 };
        // If leaf, stop
        if node.children.is_empty() { return node; }
        // Terminate if three successive conflicting ancestor blocks
        if consec >= 3 { return node; }
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
                    let child_all_ancestral = c_total > 0 && c_derived == 0 && c_ancestral == c_total;
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
                .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
                .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
        };
        // For each child, find best downstream leaf and whether any derived evidence exists along that subtree
        let mut best_child: Option<&Haplogroup> = None;
        let mut best_child_leaf: Option<HaplogroupResult> = None;
        for child in &node.children {
            let child_entry = map.get(&child.name);
            let (cd, ca, ct) = child_entry.map(|r| (r.matching_snps as i32, r.ancestral_matches as i32, r.total_snps as i32)).unwrap_or((0, 0, 0));
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
                        if entry.matching_snps > 0 { any_downstream_derived = true; }
                        match &best_leaf {
                            None => best_leaf = Some(entry.clone()),
                            Some(curr) => { if cmp_leaf(entry, curr).is_lt() { best_leaf = Some(entry.clone()); } }
                        }
                    }
                }
            }
            // Parent gate: if parent is strictly all-ancestral, do not descend
            let can_descend_parent_gate = !(parent_all_ancestral);
            // Decide if we can descend into this child
            let can_descend = can_descend_parent_gate && (clear_derived_here || (single_ancestral && any_downstream_derived));
            if can_descend {
                if let Some(bl) = best_leaf {
                    match (&best_child_leaf, best_child) {
                        (None, _) => { best_child_leaf = Some(bl); best_child = Some(child); },
                        (Some(curr_leaf), Some(_)) => { if cmp_leaf(&bl, curr_leaf).is_lt() { best_child_leaf = Some(bl); best_child = Some(child); } },
                        _ => {}
                    }
                }
            }
        }
        // Only allow ONE extra level of descent from the current node
        if let Some(ch) = best_child { return ch; }
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
                        Some(curr) => { if cmp(leaf, curr).is_lt() { best_under_chosen = Some(leaf.clone()); } }
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
                                    if cmp(leaf, curr).is_lt() { best_in_child = Some(leaf.clone()); }
                                }
                            }
                        }
                    }
                }
                if let Some(candidate) = best_in_child {
                    match &best_leaf_overall {
                        None => best_leaf_overall = Some(candidate),
                        Some(curr) => {
                            if cmp(&candidate, curr).is_lt() { best_leaf_overall = Some(candidate); }
                        }
                    }
                }
            }
        }
    }

    // Optional debug on key splits (IJK -> IJ vs K, K -> R)
    if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
        // Helper: find first node matching a name prefix or exact known markers
        fn find_by_pred<'a, F: Fn(&str) -> bool>(node: &'a Haplogroup, pred: &F) -> Option<&'a Haplogroup> {
            if pred(&node.name) { return Some(node); }
            for c in &node.children { if let Some(found) = find_by_pred(c, pred) { return Some(found); } }
            None
        }
        // Helper: pretty print tallies and decisions
        let mut print_children = |label: &str, parent: &Haplogroup| {
            eprintln!("[debug] {} split at {}: {} children", label, parent.name, parent.children.len());
            // Determine which child would be chosen by rules from this parent
            let chosen = pick_best_by_rules(parent, &unique_results, &remaining, tree, 0);
            for ch in &parent.children {
                if let Some(entry) = unique_results.get(&ch.name) {
                    let total = entry.total_snps;
                    let mut derived = entry.matching_snps;
                    let ancestral = entry.ancestral_matches;
                    let no_calls = entry.no_calls;
                    let norm_single = total > 5 && derived == 1;
                    if norm_single { derived = 0; }
                    let conflict = total > 4 && (ancestral as u32) >= ((total + 1) / 2);
                    let clear_here = ((entry.matching_snps >= 1 && entry.ancestral_matches == 0)
                        || (entry.matching_snps >= 2 && entry.matching_snps > entry.ancestral_matches));
                    let mark = if ch.name == chosen.name { " <chosen>" } else { "" };
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
        if let Some(ijk) = find_by_pred(tree, &|n: &str| n.starts_with("IJK-") || n.contains("IJK")) {
            print_children("IJK", ijk);
        }
        // K split
        if let Some(k_node) = find_by_pred(tree, &|n: &str| n.starts_with("K-") || n == "K" || n.contains("K-M9")) {
            print_children("K", k_node);
        }
    }

    // Take the top LEAF from the best root child path and output the LEAF first, then its ancestral path
    if let Some(top_result) = best_leaf_overall.or_else(|| remaining.first().cloned()) {
        // Output the top leaf first
        if let Some(entry) = unique_results.get(&top_result.name).cloned() {
            if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
                eprintln!("[debug] PATH {}\tbranch={:.4}\tcum={:.4}", entry.name, entry.score, entry.score);
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
                        eprintln!("[debug] PATH {}\tbranch={:.4}\tcum={:.4}", entry.name, entry.score, cumulative);
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

    // Sort remaining LEAVES by score descending first, cumulative SNPs as tiebreaker
    remaining.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    });

    // Add remaining results
    ordered_results.extend(remaining);

    ordered_results
}
