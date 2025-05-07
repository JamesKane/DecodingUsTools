use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Deserialize, Debug)]
pub struct Snp {
    #[serde(alias = "pos")]
    position: u32,
    ancestral: String,
    derived: String,
    #[serde(default = "default_chromosome")]
    chromosome: String,
    #[serde(default = "default_build")]
    build: String,
}

fn default_chromosome() -> String {
    "MT".to_string()
}

fn default_build() -> String {
    "rCRS".to_string()
}

#[derive(Deserialize, Debug)]
pub struct Haplogroup {
    name: String,
    #[serde(default)]
    parent: Option<String>,
    snps: Vec<Snp>,
    #[serde(default)]
    children: Vec<Haplogroup>,
}

struct BranchResult {
    has_excess_ancestral: bool,
    previous_had_excess: bool,
}

pub fn analyze_haplogroup(
    bam_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    tree_type: TreeType,
) -> Result<(), Box<dyn std::error::Error>> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap(),
    );

    let mut bam = bam::Reader::from_path(&bam_file)?;
    validate_hg38_reference(&bam)?;

    // Get tree from cache
    let tree_cache = TreeCache::new(tree_type)?;
    let tree: Haplogroup = tree_cache.get_tree()?;

    progress.set_message("Processing BAM file...");

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Snp)>> = HashMap::new();
    collect_snps(&tree, &mut positions);

    // Process BAM file
    let mut pileup = bam.pileup();
    pileup.set_max_depth(1000000);

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    for p in pileup {
        let pileup = p?;
        let pos = pileup.pos() as u32;

        if positions.contains_key(&pos) && pileup.depth() >= min_depth {
            let bases: Vec<_> = pileup
                .alignments()
                .filter(|aln| aln.record().mapq() >= min_quality)
                .filter_map(|aln| Option::from(aln.record().seq()[aln.qpos()?] as char))
                .collect();

            if !bases.is_empty() {
                let mut base_counts: HashMap<char, u32> = HashMap::new();
                for base in bases {
                    *base_counts.entry(base).or_insert(0) += 1;
                }

                if let Some((base, count)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                    let freq = *count as f64 / pileup.depth() as f64;
                    if freq >= 0.8 {
                        snp_calls.insert(pos, (*base, pileup.depth(), freq));
                    }
                }
            }
        }
    }

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups
    let mut scores = Vec::new();
    calculate_haplogroup_score(&tree, &snp_calls, &mut scores, false);

    // Sort results
    scores.sort_by(|a, b| {
        let score_diff = b.1.partial_cmp(&a.1).unwrap();
        if (b.1 - a.1).abs() < 0.1 {
            b.3.cmp(&a.3)
        } else {
            score_diff
        }
    });

    // Write results
    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(writer, "Haplogroup\tScore\tMatching_SNPs\tTotal_SNPs")?;
    for (name, score, matching, total) in scores {
        writeln!(writer, "{}\t{:.4}\t{}\t{}", name, score, matching, total)?;
    }

    progress.finish_with_message("Analysis complete!");
    Ok(())
}

fn collect_snps<'a>(
    haplogroup: &'a Haplogroup,
    positions: &mut HashMap<u32, Vec<(&'a str, &'a Snp)>>,
) {
    for snp in &haplogroup.snps {
        if is_valid_snp(snp) {
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

fn is_valid_snp(snp: &Snp) -> bool {
    match snp.chromosome.as_str() {
        "MT" => true,
        "Y" => snp.build == "hg38",
        _ => false,
    }
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
    let valid_snps: Vec<_> = haplogroup
        .snps
        .iter()
        .filter(|snp| is_valid_snp(snp))
        .collect();
    let total_snps = valid_snps.len();

    for snp in valid_snps {
        if let Some((called_base, _, _)) = snp_calls.get(&snp.position) {
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
        ancestral_ratio > 0.2
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

fn validate_hg38_reference(bam: &bam::Reader) -> Result<(), Box<dyn std::error::Error>> {
    let header_view = bam.header();

    // Get sequence dictionary
    let sequences: HashMap<&[u8], u64> = header_view
        .target_names()
        .iter()
        .enumerate()
        .filter_map(|(tid, name)| header_view.target_len(tid as u32).map(|len| (*name, len)))
        .collect();

    // Check chrY length
    let chr_y_len = sequences
        .get(&b"chrY"[..])
        .ok_or("chrY not found in BAM header")?;
    if *chr_y_len != 57227415 {
        return Err("BAM file appears to not be aligned to hg38: chrY length mismatch".into());
    }

    // Check chrM length
    let chr_m_len = sequences
        .get(&b"chrM"[..])
        .ok_or("chrM not found in BAM header")?;
    if *chr_m_len != 16569 {
        return Err("BAM file appears to not be aligned to hg38: chrM length mismatch".into());
    }

    // Check for presence of chrY_KI270740v1_random
    if !sequences.contains_key(&b"chrY_KI270740v1_random"[..]) {
        return Err(
            "BAM file appears to not be aligned to hg38: missing chrY_KI270740v1_random".into(),
        );
    }

    Ok(())
}
