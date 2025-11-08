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
    let tree_build_id = match tree_type {
        TreeType::YDNA => "GRCh38", // FTDNA provides GRCh38; DecodingUs also supports others but GRCh38 is common
        TreeType::MTDNA => "rCRS",
    };

    // Create map of positions to check in TREE coordinates
    let mut tree_positions: HashMap<u32, Vec<(&str, &Locus)>> = HashMap::new();
    tree::collect_snps(&haplogroup_tree, &mut tree_positions, tree_build_id);

    // Optional: focus debug on a specific site by locus name (e.g., M526)
    let debug_site = std::env::var("DECODINGUS_DEBUG_SITE").ok();
    let mut debug_tree_sites: Vec<u32> = Vec::new();
    if let Some(ref site_tag) = debug_site {
        for (pos, entries) in tree_positions.iter() {
            if entries.iter().any(|(_, locus)| locus.name.contains(site_tag)) {
                debug_tree_sites.push(*pos);
            }
        }
        if !debug_tree_sites.is_empty() {
            eprintln!("[debug] Tracking sites by locus tag '{}': {:?}", site_tag, debug_tree_sites);
        }
    }

    // If BAM build differs from tree build, liftover positions to BAM coordinates; also build reverse map
    let mut positions_bam: HashMap<u32, Vec<(&str, &Locus)>> = HashMap::new();
    let mut bam_to_tree_pos: HashMap<u32, u32> = HashMap::new();
    let mut tree_to_bam_pos: HashMap<u32, u32> = HashMap::new();

    if bam_build.name() != tree_build_id {
        // Parse tree_build_id -> ReferenceGenome
        let source = match tree_build_id {
            "GRCh38" => ReferenceGenome::GRCh38,
            "GRCh37" => ReferenceGenome::GRCh37,
            _ => ReferenceGenome::GRCh38, // fallback
        };
        // For mitochondrial DNA, UCSC chain files are not reliable across accessions
        // (e.g., rCRS NC_012920 vs. UCSC NC_001807 label). Coordinates are the same length
        // and typically align 1:1, so use identity liftover while allowing chromosome aliasing.
        let lifter = if tree_type == TreeType::MTDNA {
            Liftover::identity()
        } else {
            Liftover::load_or_fetch(source.clone(), bam_build.clone())?
        };
        // Determine source chromosome name to use for liftover
        let src_chrom = if tree_type == TreeType::YDNA { "chrY" } else { "chrM" };

        // Optional verbose liftover diagnostics for Y/MT mapping.
        let verbose_liftover = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().as_deref() == Some("1");
        let mut mapped_count: u32 = 0;
        let mut unmapped_count: u32 = 0;
        let mut sample_hits: Vec<(u32, Option<u32>)> = Vec::new();

        if tree_type == TreeType::MTDNA {
            // Original MT logic with relaxed window around base position
            for (pos_tree, entries) in tree_positions.iter() {
                let base = lifter.map_pos(src_chrom, *pos_tree).unwrap_or(*pos_tree);
                let start = base.saturating_sub(5);
                let end = base + 5;
                for pos_bam in start..=end {
                    positions_bam.entry(pos_bam).or_default().extend(entries.iter().copied());
                    // Only map back to the same tree position; first write wins
                    bam_to_tree_pos.entry(pos_bam).or_insert(*pos_tree);
                }
                if verbose_liftover && sample_hits.len() < 10 {
                    sample_hits.push((*pos_tree, Some(base)));
                }
                mapped_count += 1; // treat MT as mapped via relaxed window
            }
        } else {
            // YDNA: batch liftover for performance
            let mut positions: Vec<u32> = tree_positions.keys().copied().collect();
            positions.sort_unstable();
            let mapped_vec = lifter.map_many(src_chrom, &positions);
            for (pos_tree, mapped) in positions.into_iter().zip(mapped_vec.into_iter()) {
                if let Some(pos_bam) = mapped {
                    if let Some(entries) = tree_positions.get(&pos_tree) {
                        positions_bam.entry(pos_bam).or_default().extend(entries.iter().copied());
                    }
                    bam_to_tree_pos.entry(pos_bam).or_insert(pos_tree);
                    tree_to_bam_pos.entry(pos_tree).or_insert(pos_bam);
                    mapped_count += 1;
                } else {
                    unmapped_count += 1;
                }
                if verbose_liftover && sample_hits.len() < 10 {
                    sample_hits.push((pos_tree, mapped));
                }
            }
        }

        if verbose_liftover {
            eprintln!(
                "[liftover] {} -> {} on {}: mapped={}, unmapped={}, examples={:?}",
                tree_build_id,
                bam_build.name(),
                src_chrom,
                mapped_count,
                unmapped_count,
                sample_hits
            );
        }
    } else {
        positions_bam = tree_positions.clone();
        for (p, _) in &positions_bam { bam_to_tree_pos.insert(*p, *p); }
    }

    let mut snp_calls_bam: HashMap<u32, (char, u32, f64)> = HashMap::new();

    // Collect SNP calls against BAM coordinates
    caller::collect_snp_calls(
        min_depth,
        min_quality,
        &mut bam,
        &header,
        bam_build.name().to_string(),
        chromosome.clone(),
        &mut positions_bam,
        &mut snp_calls_bam,
    )?;

    // Remap SNP calls back to tree coordinates so scoring can use tree_build_id
    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();
    for (bam_pos, call) in snp_calls_bam.into_iter() {
        if let Some(tree_pos) = bam_to_tree_pos.get(&bam_pos) {
            snp_calls.insert(*tree_pos, call);
        }
    }

    // Targeted debug for a specific locus tag (e.g., M526) without probes
    if let Some(ref site_tag) = debug_site {
        // Reopen BAM for a fresh read over specific positions
        if !debug_tree_sites.is_empty() {
            let mut bam_dbg = bam::Reader::from_path(&bam_file)?;
            let header_dbg = bam_dbg.header().clone();
            for tp in &debug_tree_sites {
                if let Some(bp) = tree_to_bam_pos.get(tp) {
                    // Tally base counts at the mapped BAM position
                    let counts = caller::tally_bases_at_positions(
                        &mut bam_dbg,
                        &header_dbg,
                        &chromosome,
                        std::iter::once(*bp),
                        min_quality,
                    )?;
                    let counts_str = if let Some(c) = counts.get(bp) {
                        let mut v: Vec<(char, u32)> = c.iter().map(|(b, n)| (*b, *n)).collect();
                        v.sort_by_key(|(b, _)| *b);
                        format!("{:?}", v)
                    } else { "None".to_string() };
                    let call_bam = snp_calls.get(tp).cloned();
                    eprintln!(
                        "[debug] Site {} tree_pos={} -> bam_pos={} counts={} call={:?}",
                        site_tag, tp, bp, counts_str, call_bam
                    );
                    // Print tree allele annotations at this site
                    if let Some(entries) = tree_positions.get(tp) {
                        for (_, locus) in entries {
                            if let Some(coord) = locus.coordinates.get(tree_build_id) {
                                eprintln!(
                                    "[debug] Locus {} alleles (tree={}): ancestral={}, derived={}",
                                    locus.name, tree_build_id, coord.ancestral, coord.derived
                                );
                            }
                        }
                    }
                } else {
                    eprintln!("[debug] Site {} tree_pos={} did not liftover to BAM", site_tag, tp);
                }
            }
        }
    }

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups with tree coordinates
    let mut scores = Vec::new();
    scoring::calculate_haplogroup_score(
        &haplogroup_tree,
        &snp_calls,
        &mut scores,
        None,
        0,
        tree_build_id,
    );

    let ordered_scores = collect_scored_paths(scores, &haplogroup_tree);

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
                get_snp_details(&result.name, &haplogroup_tree, &snp_calls, tree_build_id);
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
            let is_leaf = leaf_names.contains(&result.name);
            // Keep only results that:
            // 1. Are leaves
            // 2. Either: (a) clean sparse signal (ancestral == 0 and derived >= 1), or (b) consistent signal (derived - 2*ancestral >= 1 and derived >= 2)
            let consistency = (result.matching_snps as i64) - 2 * (result.ancestral_matches as i64);
            let clean_sparse = result.ancestral_matches == 0 && result.matching_snps >= 1;
            let consistent = consistency >= 1 && result.matching_snps >= 2;
            is_leaf && (clean_sparse || consistent)
        })
        .collect();

    // Sort by consistency first (no-call tolerant), then by score, then by cumulative SNPs
    remaining.sort_by(|a, b| {
        let a_consistency = (a.matching_snps as i64) - 2 * (a.ancestral_matches as i64);
        let b_consistency = (b.matching_snps as i64) - 2 * (b.ancestral_matches as i64);
        b_consistency
            .cmp(&a_consistency)
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    });

    // Helper: comparator prioritizing no-call tolerance and consistency
    let cmp = |a: &HaplogroupResult, b: &HaplogroupResult| {
        let a_consistency = (a.matching_snps as i64) - 2 * (a.ancestral_matches as i64);
        let b_consistency = (b.matching_snps as i64) - 2 * (b.ancestral_matches as i64);
        b.ancestral_matches
            .cmp(&a.ancestral_matches) // fewer ancestral preferred
            .then_with(|| b_consistency.cmp(&a_consistency))
            .then_with(|| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal))
            .then_with(|| b.cumulative_snps.cmp(&a.cumulative_snps))
    };

    // From the root, choose the child subtree whose best leaf satisfies the comparator
    let mut best_leaf_overall: Option<HaplogroupResult> = None;
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

    // Optional debug on key splits (IJK -> IJ vs K, K -> R)
    if std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") {
        // Helper: find first node matching a name prefix or exact known markers
        fn find_by_pred<'a, F: Fn(&str) -> bool>(node: &'a Haplogroup, pred: &F) -> Option<&'a Haplogroup> {
            if pred(&node.name) { return Some(node); }
            for c in &node.children { if let Some(found) = find_by_pred(c, pred) { return Some(found); } }
            None
        }
        // IJK split
        if let Some(ijk) = find_by_pred(tree, &|n: &str| n.starts_with("IJK-") || n.contains("IJK")) {
            eprintln!("[debug] IJK split at {}: {} children", ijk.name, ijk.children.len());
            for ch in &ijk.children {
                if let Some(entry) = unique_results.get(&ch.name) {
                    let consistency = (entry.matching_snps as i64) - 2 * (entry.ancestral_matches as i64);
                    eprintln!("[debug] IJK child {}\t{:.4}\t{}\t{}\tconsistency={}", entry.name, entry.score, entry.matching_snps, entry.ancestral_matches, consistency);
                } else {
                    eprintln!("[debug] IJK child {} has no score entry", ch.name);
                }
            }
        }
        // K split
        if let Some(k_node) = find_by_pred(tree, &|n: &str| n.starts_with("K-") || n == "K" || n.contains("K-M9")) {
            eprintln!("[debug] K split at {}: {} children", k_node.name, k_node.children.len());
            for ch in &k_node.children {
                if let Some(entry) = unique_results.get(&ch.name) {
                    let consistency = (entry.matching_snps as i64) - 2 * (entry.ancestral_matches as i64);
                    eprintln!("[debug] K child {}\t{:.4}\t{}\t{}\tconsistency={}", entry.name, entry.score, entry.matching_snps, entry.ancestral_matches, consistency);
                } else {
                    eprintln!("[debug] K child {} has no score entry", ch.name);
                }
            }
        }
    }

    // Take the top LEAF from the best root child path and ensure its ancestral path is included first
    if let Some(top_result) = best_leaf_overall.or_else(|| remaining.first().cloned()) {
        if let Some(path) = tree::find_path_to_root(tree, &top_result.name, &[]) {
            // Insert path elements (root -> leaf) first
            let mut cumulative = 0.0;
            for name in path.clone().into_iter().rev() {
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
