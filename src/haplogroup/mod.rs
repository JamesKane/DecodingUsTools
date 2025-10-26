mod caller;
mod scoring;
mod tree;
pub(crate) mod types;
mod validation;

use crate::haplogroup::types::{Haplogroup, HaplogroupResult, Locus};
use crate::utils::cache::TreeType;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::IndexedReader;
use rust_htslib::faidx;
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
    // Initialize the FASTA reader
    let mut fasta_reader = faidx::Reader::from_path(&reference_file)?;

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

    // Use IndexedReader instead of Reader
    let bam = IndexedReader::from_path(&bam_file)?;
    let (genome, chromosome) = validation::validate_reference(&bam, tree_type)?;

    let haplogroup_tree = tree::load_tree(tree_type, provider)?;

    // Determine build ID and chromosome based on genome and tree type
    let build_id = match tree_type {
        TreeType::YDNA => genome.name(),
        TreeType::MTDNA => "rCRS",
    };

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Locus)>> = HashMap::new();
    tree::collect_snps(&haplogroup_tree, &mut positions, build_id);

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    // Pass both build_id and actual chromosome accession to collect_snp_calls
    caller::collect_snp_calls(
        min_depth,
        min_quality,
        &mut fasta_reader,
        bam,
        build_id.parse().unwrap(),
        chromosome,
        &mut positions,
        &mut snp_calls,
    )?;

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups
    let mut scores = Vec::new();
    scoring::calculate_haplogroup_score(
        &haplogroup_tree,
        &snp_calls,
        &mut scores,
        None,
        0,
        build_id,
    );

    let ordered_scores = collect_scored_paths(scores, &haplogroup_tree);

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
                get_snp_details(&result.name, &haplogroup_tree, &snp_calls, build_id);
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

    // Convert to vec and filter with stricter criteria
    let mut remaining: Vec<HaplogroupResult> = unique_results
        .into_values()
        .filter(|result| {
            // Keep only results that:
            // 1. Have a non-zero score
            // 2. Don't have overwhelmingly more ancestral than derived matches
            // 3. Have at least some derived matches
            result.score > 0.0
                && result.ancestral_matches <= result.matching_snps * 3
                && result.matching_snps > 0
        })
        .collect();

    // Sort by cumulative SNPs descending, then by score descending for ties
    remaining.sort_by(|a, b| {
        b.cumulative_snps.cmp(&a.cumulative_snps).then_with(|| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    // Take the top result and ensure its ancestral path is included first
    if let Some(top_result) = remaining.first().cloned() {
        if let Some(path) = tree::find_path_to_root(tree, &top_result.name, &remaining) {
            // Remove all path elements from remaining
            for name in &path {
                if let Some(pos) = remaining.iter().position(|r| &r.name == name) {
                    ordered_results.push(remaining.remove(pos));
                }
            }
        }
    }

    // Sort remaining results by cumulative SNPs descending, score as tiebreaker
    remaining.sort_by(|a, b| {
        b.cumulative_snps.cmp(&a.cumulative_snps).then_with(|| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    });

    // Add remaining results
    ordered_results.extend(remaining);

    ordered_results
}
