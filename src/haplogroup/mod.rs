mod caller;
mod scoring;
mod tree;
pub(crate) mod types;
mod validation;

use crate::haplogroup::types::{Haplogroup, HaplogroupResult, Snp};
use crate::utils::cache::TreeType;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
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
) -> Result<(), Box<dyn Error>> {
    // Initialize the FASTA reader
    let mut fasta_reader = bio::io::fasta::IndexedReader::from_file(&reference_file)?;

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );

    progress.set_message("Validating BAM reference genome...");

    // Use IndexedReader instead of Reader
    let bam = IndexedReader::from_path(&bam_file)?;
    validation::validate_hg38_reference(&bam)?;

    let haplogroup_tree = tree::load_tree(tree_type)?;

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Snp)>> = HashMap::new();
    tree::collect_snps(&haplogroup_tree, &mut positions);

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    caller::collect_snp_calls(
        min_depth,
        min_quality,
        &mut fasta_reader,
        bam,
        &mut positions,
        &mut snp_calls,
    )?;

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups
    let mut scores = Vec::new();
    scoring::calculate_haplogroup_score(&haplogroup_tree, &snp_calls, &mut scores, None, 0);

    let ordered_scores = collect_scored_paths(scores, &haplogroup_tree);

    // Use ordered_scores for writing the report
    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(
        writer,
        "Haplogroup\tScore\tMatching_SNPs\tMismatching_SNPs\tAncestral_Matches\tNo_Calls\tTotal_SNPs\tCumulative_SNPs\tDepth"
    )?;
    for result in ordered_scores {
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

    progress.finish_with_message("Analysis complete!");
    Ok(())
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
