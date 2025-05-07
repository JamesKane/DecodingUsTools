use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Deserialize, Debug)]
struct MtSnp {
    pos: u32,
    ancestral: String,
    derived: String,
}

#[derive(Deserialize, Debug)]
struct MtHaplogroup {
    name: String,
    snps: Vec<MtSnp>,
    #[serde(default)]
    children: Vec<MtHaplogroup>,
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
    // Get MT tree from cache or download
    let tree_cache = TreeCache::new(TreeType::MTDNA)?;
    let mt_tree: MtHaplogroup = tree_cache.get_tree()?;

    // Open BAM file
    let mut bam = bam::Reader::from_path(&bam_file)?;

    // Process MT chromosome positions
    let mut pileup = bam.pileup();
    pileup.set_max_depth(1000000);

    let progress = ProgressBar::new_spinner();
    progress.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} {msg}")
        .unwrap());
    progress.set_message("Processing MT-DNA positions...");

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    for p in pileup {
        let pileup = p?;
        let pos = pileup.pos();
        let depth = pileup.depth();

        if depth >= min_depth {
            let bases: Vec<_> = pileup
                .alignments()
                .filter(|aln| aln.record().mapq() >= min_quality)
                .filter_map(|aln| Option::from(aln.record().seq()[aln.qpos()?] as char))
                .collect();

            if !bases.is_empty() {
                // Count base frequencies
                let mut base_counts: HashMap<char, u32> = HashMap::new();
                for base in bases {
                    *base_counts.entry(base).or_insert(0) += 1;
                }

                // Find most common base
                if let Some((base, count)) = base_counts.iter()
                    .max_by_key(|&(_, count)| count)
                {
                    let freq = *count as f64 / depth as f64;
                    if freq >= 0.8 {  // 80% threshold for calling a variant
                        snp_calls.insert(pos, (*base, depth, freq));
                    }
                }
            }
        }
    }

    progress.finish_with_message("MT-DNA positions processed");

    // Score haplogroups
    let mut scores = Vec::new();
    calculate_haplogroup_score(&mt_tree, &snp_calls, &mut scores, false);

    // Sort and write results
    scores.sort_by(|a, b| {
        let score_diff = b.1.partial_cmp(&a.1).unwrap();
        if (b.1 - a.1).abs() < 0.1 {
            b.3.cmp(&a.3)  // Prefer haplogroups with more total SNPs when scores are close
        } else {
            score_diff
        }
    });

    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(writer, "Haplogroup\tScore\tMatching_SNPs\tTotal_SNPs")?;
    for (name, score, matching, total) in scores {
        writeln!(writer, "{}\t{:.3}\t{}\t{}", name, score, matching, total)?;
    }

    Ok(())
}

fn calculate_haplogroup_score(
    haplogroup: &MtHaplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<(String, f64, u32, u32)>,
    previous_had_excess: bool,
) -> BranchResult {
    let mut matching_snps = 0;
    let mut ancestral_snps = 0;
    let mut called_snps = 0;
    let total_snps = haplogroup.snps.len();

    for snp in &haplogroup.snps {
        if let Some((called_base, _, _)) = snp_calls.get(&snp.pos) {
            called_snps += 1;
            if *called_base == snp.derived.chars().next().unwrap() {
                matching_snps += 1;
            } else if *called_base == snp.ancestral.chars().next().unwrap() {
                ancestral_snps += 1;
            }
        }
    }

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

    let has_excess_ancestral = if called_snps > 0 {
        let ancestral_ratio = ancestral_snps as f64 / called_snps as f64;
        ancestral_ratio > 0.2  // More than 20% ancestral
    } else {
        false
    };

    if has_excess_ancestral && previous_had_excess {
        return BranchResult {
            has_excess_ancestral: true,
            previous_had_excess: true,
        };
    }

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