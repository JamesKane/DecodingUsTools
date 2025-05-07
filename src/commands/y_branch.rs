use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::utils::cache::TreeCache;

#[derive(Deserialize, Debug)]
struct SNP {
    position: u32,
    ancestral: String,
    derived: String,
    chromosome: String,
    build: String,
}

#[derive(Deserialize, Debug)]
struct Haplogroup {
    name: String,
    parent: Option<String>,
    snps: Vec<SNP>,
    children: Vec<Haplogroup>,
}

struct BranchResult {
    has_excess_ancestral: bool,
    previous_had_excess: bool,
}

pub fn run(
    bam_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );
    let tree_cache = TreeCache::new()?;
    let tree = tree_cache.get_tree()?;

    progress.set_message("Processing BAM file...");

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &SNP)>> = HashMap::new();
    collect_snps(&tree, &mut positions);

    // Process BAM file
    let mut bam = bam::Reader::from_path(&bam_file)?;
    let mut pileup = bam.pileup();
    pileup.set_max_depth(1000);

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    for p in pileup {
        let pileup = p?;
        let pos = pileup.pos() as u32;

        if let Some(snp_info) = positions.get(&pos) {
            if pileup.depth() >= min_depth {
                let base_counts = count_bases(&pileup, min_quality);
                if let Some((base, count, frequency)) =
                    get_dominant_base(base_counts, pileup.depth())
                {
                    snp_calls.insert(pos, (base, count, frequency));
                }
            }
        }
    }

    progress.set_message("Scoring haplogroups...");
    let scores = score_haplogroups(&tree, &snp_calls)?;

    // Write results
    let output_file = File::create(output_file)?;
    let mut writer = BufWriter::new(output_file);

    writeln!(writer, "Haplogroup\tScore\tMatching_SNPs\tTotal_SNPs")?;
    for (haplogroup, score, matching, total) in scores.iter().take(10) {
        writeln!(
            writer,
            "{}\t{:.4}\t{}\t{}",
            haplogroup, score, matching, total
        )?;
    }

    progress.finish_with_message("Analysis complete!");
    Ok(())
}

fn collect_snps<'a>(
    haplogroup: &'a Haplogroup,
    positions: &mut HashMap<u32, Vec<(&'a str, &'a SNP)>>,
) {
    for snp in &haplogroup.snps {
        if snp.chromosome == "Y" && snp.build == "hg38" {
            positions
                .entry(snp.position)
                .or_default()
                .push((&haplogroup.name, snp));
        }
    }

    for child in &haplogroup.children {
        collect_snps(child, positions);
    }
}

fn count_bases(pileup: &bam::pileup::Pileup, min_quality: u8) -> HashMap<char, u32> {
    let mut counts: HashMap<char, u32> = HashMap::new();

    for alignment in pileup.alignments() {
        if alignment.record().mapq() >= min_quality {
            if let Some(qpos) = alignment.qpos() {
                let seq = alignment.record().seq().as_bytes();
                let base = seq[qpos] as char;
                *counts.entry(base).or_insert(0) += 1;
            }
        }
    }

    counts
}

fn get_dominant_base(counts: HashMap<char, u32>, total_depth: u32) -> Option<(char, u32, f64)> {
    counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(base, count)| (base, count, count as f64 / total_depth as f64))
}

fn calculate_haplogroup_score(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<(String, f64, u32, u32)>,
    previous_had_excess: bool,
) -> BranchResult {
    let mut matching_snps = 0;
    let mut ancestral_snps = 0;
    let mut called_snps = 0;
    let total_snps = haplogroup.snps.len();

    // First pass: count matching and ancestral SNPs
    for snp in &haplogroup.snps {
        if let Some((called_base, _, _)) = snp_calls.get(&snp.position) {
            called_snps += 1;
            if *called_base == snp.derived.chars().next().unwrap() {
                matching_snps += 1;
            } else if *called_base == snp.ancestral.chars().next().unwrap() {
                ancestral_snps += 1;
            }
        }
    }

    // Calculate score and check ancestral ratio
    let score = if total_snps > 0 {
        matching_snps as f64 / total_snps as f64
    } else {
        0.0
    };

    scores.push((
        haplogroup.name.clone(),
        score,
        matching_snps,
        total_snps as u32,
    ));

    // Check if we have too many ancestral alleles
    let has_excess_ancestral = if called_snps > 0 {
        let ancestral_ratio = ancestral_snps as f64 / called_snps as f64;
        ancestral_ratio > 0.2 // More than 20% ancestral
    } else {
        false
    };

    // Stop processing children only if we have excess ancestral SNPs in two successive branches
    if has_excess_ancestral && previous_had_excess {
        return BranchResult {
            has_excess_ancestral: true,
            previous_had_excess: true,
        };
    }

    // Process children
    let mut child_result = BranchResult {
        has_excess_ancestral: false,
        previous_had_excess: false,
    };

    for child in &haplogroup.children {
        child_result = calculate_haplogroup_score(child, snp_calls, scores, has_excess_ancestral);
        if child_result.has_excess_ancestral && child_result.previous_had_excess {
            break;
        }
    }

    BranchResult {
        has_excess_ancestral: child_result.has_excess_ancestral,
        previous_had_excess: has_excess_ancestral,
    }
}


fn score_haplogroups(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
) -> Result<Vec<(String, f64, u32, u32)>, Box<dyn std::error::Error>> {
    let mut scores = Vec::new();
    calculate_haplogroup_score(haplogroup, snp_calls, &mut scores, false);

    // Sort by score descending but prefer haplogroups with more total SNPs when scores are close
    scores.sort_by(|a, b| {
        let score_diff = b.1.partial_cmp(&a.1).unwrap();
        if (b.1 - a.1).abs() < 0.1 {
            // If scores are within 10% of each other, prefer the one with more total SNPs
            b.3.cmp(&a.3)
        } else {
            score_diff
        }
    });

    Ok(scores)
}
