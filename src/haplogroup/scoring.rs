use crate::haplogroup::types::{Haplogroup, HaplogroupResult, HaplogroupScore, LociType};
use std::collections::{HashMap, HashSet};

const HIGH_QUALITY_THRESHOLD: f64 = 0.7;
const MEDIUM_QUALITY_THRESHOLD: f64 = 0.5;
// Respect caller-level min depth by keeping scoring permissive; caller already filters by depth.
const MIN_DEPTH: u32 = 1;

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

fn compute_branch_metrics(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    build_id: &str,
    depth: u32,
) -> BranchMetrics {
    let defining_loci: Vec<_> = haplogroup
        .loci
        .iter()
        .filter(|locus| matches!(locus.loci_type, LociType::SNP) && locus.coordinates.contains_key(build_id))
        .collect();

    let mut bm = BranchMetrics::default();
    
    // Enable debug for P312 specifically
    let debug = std::env::var("DECODINGUS_DEBUG_SCORES").ok().as_deref() == Some("1") 
        && haplogroup.name.contains("P312");

    for locus in &defining_loci {
        if let Some(coord) = locus.coordinates.get(build_id) {
            if debug {
                // Always print position being checked for P312
                eprintln!(
                    "[scoring.debug] {} locus={} checking position={} (chrom={})",
                    haplogroup.name, locus.name, coord.position, coord.chromosome
                );
            }
            
            if let Some((called_base, depth_call, freq)) = snp_calls.get(&coord.position) {
                if debug {
                    eprintln!(
                        "[scoring.debug] {} locus={} pos={} called_base='{}' depth={} freq={:.3}",
                        haplogroup.name, locus.name, coord.position, called_base, depth_call, freq
                    );
                    eprintln!(
                        "[scoring.debug]   coord.derived='{}' coord.ancestral='{}'",
                        coord.derived, coord.ancestral
                    );
                }
                
                if *depth_call >= MIN_DEPTH {
                    let maybe_derived = coord.derived.chars().next();
                    let maybe_ancestral = coord.ancestral.chars().next();
                    
                    if debug {
                        eprintln!(
                            "[scoring.debug]   derived_char={:?} ancestral_char={:?}",
                            maybe_derived, maybe_ancestral
                        );
                    }
                    
                    if let (Some(derived_base), Some(ancestral_base)) = (maybe_derived, maybe_ancestral) {
                        if *called_base == derived_base {
                            if *freq >= HIGH_QUALITY_THRESHOLD && *depth_call >= 1 {
                                bm.derived += 1;
                                if debug { eprintln!("[scoring.debug]   -> counted as DERIVED (high qual)"); }
                            } else if *freq >= MEDIUM_QUALITY_THRESHOLD && *depth_call >= 1 {
                                bm.derived += 1;
                                if debug { eprintln!("[scoring.debug]   -> counted as DERIVED (medium qual)"); }
                            } else {
                                bm.low_quality += 1;
                                if debug { eprintln!("[scoring.debug]   -> counted as LOW_QUALITY"); }
                            }
                        } else if *called_base == ancestral_base {
                            if *freq >= HIGH_QUALITY_THRESHOLD && *depth_call >= 1 {
                                bm.ancestral += 1;
                                if debug { eprintln!("[scoring.debug]   -> counted as ANCESTRAL"); }
                            } else {
                                bm.low_quality += 1;
                                if debug { eprintln!("[scoring.debug]   -> counted as LOW_QUALITY (anc, low qual)"); }
                            }
                        } else {
                            // true conflict: called base is neither ancestral nor derived
                            bm.conflicts += 1;
                            if debug {
                                eprintln!(
                                    "[scoring.debug]   -> CONFLICT: called '{}' != derived '{}' != ancestral '{}'",
                                    called_base, derived_base, ancestral_base
                                );
                            }
                        }
                    } else {
                        bm.no_calls += 1;
                        if debug { eprintln!("[scoring.debug]   -> NO_CALL (missing derived/ancestral in coord)"); }
                    }
                } else {
                    bm.no_calls += 1;
                    if debug { eprintln!("[scoring.debug]   -> NO_CALL (depth < MIN_DEPTH)"); }
                }
            } else {
                bm.no_calls += 1;
                if debug { 
                    eprintln!("[scoring.debug] {} locus={} pos={} -> NO_CALL (no snp_call at this position)", 
                              haplogroup.name, locus.name, coord.position); 
                }
            }
        }
    }

    bm.total = defining_loci.len() as u32;

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

fn experimental_calculate(
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
    let bm = compute_branch_metrics(haplogroup, snp_calls, build_id, depth);
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
        if callable_frac < 0.3 { adjusted_branch *= 0.5; }
        else if callable_frac < 0.6 { adjusted_branch *= 0.8; }
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
        let cbm = compute_branch_metrics(child, snp_calls, build_id, depth + 1);
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
    let is_key_split = haplogroup.name.contains("IJK")
        || haplogroup.name.starts_with("K-")
        || haplogroup.name == "K"
        || haplogroup.name.contains("K-M9")
        || haplogroup.name.contains("R-L151");
    if debug && is_key_split {
        eprintln!("[debug] Split at {} (depth {}):", haplogroup.name, depth);
        eprintln!("[debug]   parent bm: der={} anc={} noc={} lowq={} total={} raw={:.4} adj={:.4}", bm.derived, bm.ancestral, bm.no_calls, bm.low_quality, bm.total, bm.raw_branch_score, adjusted_branch);
    }

    for (idx, (child, cbm)) in child_metrics.iter().enumerate() {
        // Adjust child branch score locally
        let mut child_adj = cbm.raw_branch_score;
        if cbm.total > 4 && cbm.ancestral as u32 >= ((cbm.total + 1) / 2) { child_adj *= dampen_majority_anc; }
        if cbm.total > 5 && cbm.derived == 1 { child_adj *= singleton_downweight; }
        let soft = exp_vals[idx] / sum_exp;
        let sibling_bonus = 1.0 + bonus_budget * soft;
        let clear_local = (cbm.derived > cbm.ancestral && cbm.derived >= 2) || (cbm.derived >= 1 && cbm.ancestral == 0);
        let coh = if clear_local && !parent_all_ancestral { coherence_bonus } else { 1.0 };

        // Recurse down the child to get its cumulative score
        let (child_score, child_cumulative) = experimental_calculate(
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

        if debug && is_key_split {
            eprintln!(
                "[debug]   child {:<20} der={} anc={} noc={} total={} raw={:.4} adj={:.4} soft={:.3} sib_bonus={:.3} coh={:.3} final_child={:.4}",
                child.name,
                cbm.derived,
                cbm.ancestral,
                cbm.no_calls,
                cbm.total,
                cbm.raw_branch_score,
                child_adj,
                soft,
                sibling_bonus,
                coh,
                child_score.score * sibling_bonus * coh
            );
        }

        // Push cumulative result for this child
        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: combined_score,
            matching_snps: combined_matches as u32,
            mismatching_snps: bm.conflicts,
            low_quality_snps: bm.low_quality,
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
    experimental_calculate(haplogroup, snp_calls, scores, parent_info, depth, build_id)
}