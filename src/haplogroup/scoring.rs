use crate::haplogroup::types::{Haplogroup, HaplogroupResult, HaplogroupScore, LociType};
use std::collections::{HashMap, HashSet};

const HIGH_QUALITY_THRESHOLD: f64 = 0.7;
const MEDIUM_QUALITY_THRESHOLD: f64 = 0.5;
const MIN_DEPTH: u32 = 4;

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
                    let derived_base = coord.derived.chars().next().unwrap();
                    let ancestral_base = coord.ancestral.chars().next().unwrap();

                    if *called_base == derived_base {
                        if *freq >= HIGH_QUALITY_THRESHOLD {
                            branch_derived += 1;
                        } else if *freq >= MEDIUM_QUALITY_THRESHOLD {
                            branch_derived += 1;
                        } else {
                            branch_low_quality += 1;
                        }
                    } else if *called_base == ancestral_base {
                        if *freq >= HIGH_QUALITY_THRESHOLD {
                            branch_ancestral += 1;
                        } else {
                            branch_low_quality += 1;
                        }
                    } else if *freq >= HIGH_QUALITY_THRESHOLD {
                        branch_derived += 1;
                    } else {
                        branch_low_quality += 1;
                    }
                } else {
                    branch_no_calls += 1;
                }
            } else {
                branch_no_calls += 1;
            }
        }
    }

    // Calculate score with improved logic
    let total_calls = branch_derived + branch_ancestral + branch_low_quality;
    if total_calls > 0 {
        // Scoring logic remains the same
        let branch_score = if branch_ancestral == 0 {
            match branch_derived {
                d if d >= 1 => 3.08,
                _ => 1.0,
            }
        } else {
            match (branch_derived, branch_ancestral) {
                (d, a) if d >= 3 && a <= d / 2 => 2.8,
                (d, a) if d >= 2 && a <= d => 2.5,
                (d, _) if d >= 2 => 2.0,
                (1, a) if a <= 2 => 1.5,
                (d, a) if a > d * 3 => 0.0,
                _ => 1.0,
            }
        };

        let quality_factor = if branch_low_quality == 0 { 1.1 } else { 0.9 };
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

        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: child_score.score,
            matching_snps: child_score.matches.try_into().unwrap_or(0),
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: child_score.ancestral_matches.try_into().unwrap_or(0),
            no_calls: child_score.no_calls.try_into().unwrap_or(0),
            total_snps: defining_loci.len() as u32,
            cumulative_snps: child_cumulative.len() as u32,
            depth,
        });
    }

    (current_score, cumulative_snps)
}