use crate::haplogroup::types::{Haplogroup, HaplogroupResult, HaplogroupScore, LociType};
use std::collections::{HashMap, HashSet};

const HIGH_QUALITY_THRESHOLD: f64 = 0.7;
const MEDIUM_QUALITY_THRESHOLD: f64 = 0.5;
// Respect caller-level min depth by keeping scoring permissive; caller already filters by depth.
const MIN_DEPTH: u32 = 1;

pub(crate) fn calculate_haplogroup_score(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<HaplogroupResult>,
    parent_info: Option<&(HaplogroupScore, HashSet<u32>)>,
    depth: u32,
    build_id: &str,
) -> (HaplogroupScore, HashSet<u32>) {
    let mut current_score = HaplogroupScore::default();
    let mut cumulative_snps = parent_info
        .map(|(_, snps)| snps.clone())
        .unwrap_or_default();

    // Filter valid loci and add their positions to cumulative_snps
    let defining_loci: Vec<_> = haplogroup
        .loci
        .iter()
        .filter(|locus| {
            matches!(locus.loci_type, LociType::SNP) &&
                locus.coordinates.contains_key(build_id)
        })
        .collect();

    // Add positions to cumulative_snps
    for locus in &defining_loci {
        if let Some(coord) = locus.coordinates.get(build_id) {
            cumulative_snps.insert(coord.position);
        }
    }

    let mut branch_derived = 0;
    let mut branch_ancestral = 0;
    let mut branch_no_calls = 0;
    let mut branch_low_quality = 0;

    // Process defining loci
    for locus in &defining_loci {
        if let Some(coord) = locus.coordinates.get(build_id) {
            if let Some((called_base, depth, freq)) = snp_calls.get(&coord.position) {
                if *depth >= MIN_DEPTH {
                    let maybe_derived = coord.derived.chars().next();
                    let maybe_ancestral = coord.ancestral.chars().next();
                    if let (Some(derived_base), Some(ancestral_base)) = (maybe_derived, maybe_ancestral) {
                        let dbg_site = std::env::var("DECODINGUS_DEBUG_SITE").ok();
                        if let Some(ref tag) = dbg_site {
                            if locus.name.contains(tag) {
                                eprintln!(
                                    "[debug] scoring locus {} pos={} called={} depth={} freq={:.2} anc={} der={}",
                                    locus.name,
                                    coord.position,
                                    called_base,
                                    depth,
                                    freq,
                                    ancestral_base,
                                    derived_base
                                );
                            }
                        }
                        if *called_base == derived_base {
                            if *freq >= HIGH_QUALITY_THRESHOLD && *depth >= 1 {
                                branch_derived += 1;
                            } else if *freq >= MEDIUM_QUALITY_THRESHOLD && *depth >= 1 {
                                branch_derived += 1;
                            } else {
                                branch_low_quality += 1;
                            }
                        } else if *called_base == ancestral_base {
                            if *freq >= HIGH_QUALITY_THRESHOLD && *depth >= 1 {
                                branch_ancestral += 1;
                            } else {
                                branch_low_quality += 1;
                            }
                        } else {
                            // Called base is neither derived nor ancestral; treat as mismatch/low-quality for this branch
                            branch_low_quality += 1;
                        }
                    } else {
                        // Missing allele information; cannot evaluate this locus reliably
                        branch_no_calls += 1;
                    }
                } else {
                    branch_no_calls += 1;
                }
            } else {
                branch_no_calls += 1;
            }
        }
    }

    // Calculate score with ratio-based logic to prefer consistent evidence
    let callable = branch_derived + branch_ancestral + branch_low_quality;
    if callable > 0 {
        let ratio = (branch_derived as f64) / (callable as f64);
        // Evidence size factor: ln(1 + derived)
        let evidence = (1.0 + branch_derived as f64).ln();
        let mut branch_score = ratio * evidence;
        // Slightly prefer coherent deeper paths
        let depth_factor = 1.0 + 0.05 * (depth as f64);
        branch_score *= depth_factor;
        // Penalize conflicting ancestral evidence more strongly
        if branch_derived > 0 {
            if branch_ancestral >= branch_derived { branch_score *= 0.3; }
            if branch_ancestral > branch_derived * 2 { branch_score *= 0.15; }
        }
        // Penalize if most loci are uncallable at this branch
        let total = defining_loci.len() as u32;
        if total > 0 {
            let callable_frac = callable as f64 / total as f64;
            if callable_frac < 0.3 { branch_score *= 0.6; }
            else if callable_frac < 0.6 { branch_score *= 0.85; }
        }
        // Down-weight singleton derived hits on a branch, but keep it mostly intact for ancient/sparse data
        if branch_derived < 2 { branch_score *= 0.9; }
        // Low-quality penalty (slightly reward clean calls)
        let quality_factor = if branch_low_quality == 0 { 1.1 } else { 0.92 };
        current_score.score = branch_score * quality_factor;
    }

    current_score.matches += branch_derived;
    current_score.ancestral_matches += branch_ancestral;
    current_score.no_calls += branch_no_calls;
    current_score.total_snps += defining_loci.len();

    if branch_ancestral > branch_derived * 10 {
        scores.push(HaplogroupResult {
            name: haplogroup.name.clone(),
            score: 0.0,
            matching_snps: branch_derived.try_into().unwrap_or(0),
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: branch_ancestral.try_into().unwrap_or(0),
            no_calls: branch_no_calls.try_into().unwrap_or(0),
            total_snps: defining_loci.len() as u32,
            cumulative_snps: cumulative_snps.len() as u32,
            depth,
        });
        return (current_score, cumulative_snps);
    }

    // Push branch-only result for the current node to enable path-level aggregation
    scores.push(HaplogroupResult {
        name: haplogroup.name.clone(),
        score: current_score.score,
        matching_snps: current_score.matches as u32,
        mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
        ancestral_matches: current_score.ancestral_matches as u32,
        no_calls: current_score.no_calls as u32,
        total_snps: current_score.total_snps as u32,
        cumulative_snps: cumulative_snps.len() as u32,
        depth,
    });

    // Process children recursively
    for child in &haplogroup.children {
        let (child_score, child_cumulative) = calculate_haplogroup_score(
            child,
            snp_calls,
            scores,
            Some(&(current_score.clone(), cumulative_snps.clone())),
            depth + 1,
            build_id,
        );

        // Combine current branch evidence with child for cumulative scoring down the path
        let combined_matches = current_score.matches + child_score.matches;
        let combined_ancestral = current_score.ancestral_matches + child_score.ancestral_matches;
        let combined_no_calls = current_score.no_calls + child_score.no_calls;
        let combined_total = current_score.total_snps + child_score.total_snps;
        let combined_score = current_score.score + child_score.score;

        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: combined_score,
            matching_snps: combined_matches.try_into().unwrap_or(0),
            // We don't separately track true mismatches; approximate with low-quality from current branch
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: combined_ancestral.try_into().unwrap_or(0),
            no_calls: combined_no_calls.try_into().unwrap_or(0),
            total_snps: combined_total as u32,
            cumulative_snps: child_cumulative.len() as u32,
            depth: depth + 1,
        });
    }

    (current_score, cumulative_snps)
}