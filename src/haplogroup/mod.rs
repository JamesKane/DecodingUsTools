use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Deserialize, Debug)]
pub struct Variant {
    variant: String,
    #[serde(default)]
    pos: u32,
    #[serde(default)]
    ancestral: String,
    #[serde(default)]
    derived: String,
    #[serde(default)]
    region: String,
    #[serde(rename = "snpId", default)]
    snp_id: u32,
}

impl Variant {
    fn to_snp(&self) -> Option<Snp> {
        // Only convert variants that have all required fields
        if self.variant.is_empty() || self.ancestral.is_empty() || self.derived.is_empty() {
            return None;
        }

        Some(Snp {
            position: self.pos,
            ancestral: self.ancestral.clone(),
            derived: self.derived.clone(),
            chromosome: default_chromosome(),
            build: default_build(),
        })
    }
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupNode {
    #[serde(rename = "haplogroupId")]
    haplogroup_id: u32,
    #[serde(rename = "parentId")]
    parent_id: u32,
    name: String,
    #[serde(rename = "isRoot")]
    is_root: bool,
    root: String,
    #[serde(rename = "kitsCount")]
    kits_count: u32,
    #[serde(rename = "subBranches")]
    sub_branches: u32,
    #[serde(rename = "bigYCount")]
    big_y_count: u32,
    #[serde(default)]
    variants: Vec<Variant>,
    #[serde(default)]
    children: Vec<u32>,
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupTree {
    #[serde(rename = "allNodes")]
    all_nodes: HashMap<String, HaplogroupNode>,
}

#[derive(Debug)]
pub struct Haplogroup {
    name: String,
    parent: Option<String>,
    snps: Vec<Snp>,
    children: Vec<Haplogroup>,
}

impl HaplogroupTree {
    fn build_tree(&self, node_id: u32) -> Option<Haplogroup> {
        let node = self.all_nodes.get(&node_id.to_string())?;

        let parent = if node.parent_id != 0 {
            self.all_nodes
                .get(&node.parent_id.to_string())
                .map(|p| p.name.clone())
        } else {
            None
        };

        let snps = node.variants.iter().filter_map(|v| v.to_snp()).collect();

        let children = node
            .children
            .iter()
            .filter_map(|&child_id| self.build_tree(child_id))
            .collect();

        Some(Haplogroup {
            name: node.name.clone(),
            parent,
            snps,
            children,
        })
    }
}

#[derive(Deserialize, Debug)]
pub struct Snp {
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

    progress.set_message("Validating BAM reference genome...");

    // Use IndexedReader instead of Reader
    let mut bam = IndexedReader::from_path(&bam_file)?;
    validate_hg38_reference(&bam)?;

    // Get tree from cache
    let tree_cache = TreeCache::new(tree_type)?;
    let tree_json: HaplogroupTree = tree_cache.get_tree()?;

    let root_count = tree_json.all_nodes.values().filter(|node| node.parent_id == 0).count();
    if root_count > 1 {
        return Err("Multiple root nodes found in tree".into());
    }

    // Find the root node (usually has parent_id = 0)
    let root_node = tree_json
        .all_nodes
        .values()
        .find(|node| node.parent_id == 0)
        .ok_or("No root node found")?;

    // Build the tree structure starting from root
    let tree = tree_json
        .build_tree(root_node.haplogroup_id)
        .ok_or("Failed to build tree")?;

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Snp)>> = HashMap::new();
    collect_snps(&tree, &mut positions);

    // Determine which chromosome we need to process
    let need_y = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "Y"));
    let need_mt = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "MT"));

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    // Create a progress bar for BAM processing
    let progress_style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
        .unwrap()
        .progress_chars("#>-");

    // Process one chromosome based on the tree type
    if need_y {
        let tid = bam.header().tid(b"chrY").ok_or("chrY not found in BAM")?;
        let total_len = bam.header().target_len(tid).unwrap_or(0);

        let progress_bar = ProgressBar::new(total_len);
        progress_bar.set_style(progress_style.clone());
        progress_bar.set_message("Processing chromosome Y...");

        // Create the region string for fetch
        let region = format!("chrY:1-{}", total_len);
        bam.fetch(&region)?;
        process_region(
            &mut bam,
            &positions,
            min_depth,
            min_quality,
            &mut snp_calls,
            &progress_bar,
        )?;

        progress_bar.finish_with_message("Chromosome Y processing complete");
    } else if need_mt {
        let tid = bam.header().tid(b"chrM").ok_or("chrM not found in BAM")?;
        let total_len = bam.header().target_len(tid).unwrap_or(0);

        let progress_bar = ProgressBar::new(total_len);
        progress_bar.set_style(progress_style);
        progress_bar.set_message("Processing mitochondrial DNA...");

        // Create the region string for fetch
        let region = format!("chrM:1-{}", total_len);
        bam.fetch(&region)?;
        process_region(
            &mut bam,
            &positions,
            min_depth,
            min_quality,
            &mut snp_calls,
            &progress_bar,
        )?;

        progress_bar.finish_with_message("Mitochondrial DNA processing complete");
    }

    progress.set_message("Scoring haplogroups...");

    // Score haplogroups
    let mut scores = Vec::new();
    calculate_haplogroup_score(&tree, &snp_calls, &mut scores, false);

    // Sort and write results
    scores.sort_by(|a, b| {
        let score_diff = b.1.partial_cmp(&a.1).unwrap();
        if (b.1 - a.1).abs() < 0.1 {
            b.3.cmp(&a.3)
        } else {
            score_diff
        }
    });

    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(writer, "Haplogroup\tScore\tMatching_SNPs\tTotal_SNPs")?;
    for (name, score, matching, total) in scores {
        writeln!(writer, "{}\t{:.4}\t{}\t{}", name, score, matching, total)?;
    }

    progress.finish_with_message("Analysis complete!");
    Ok(())
}

// Helper function to process a region of the BAM file
fn process_region<R: Read>(
    bam: &mut R,
    positions: &HashMap<u32, Vec<(&str, &Snp)>>,
    min_depth: u32,
    min_quality: u8,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
    progress: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut coverage: HashMap<u32, Vec<char>> = HashMap::new();

    // First pass: collect bases at positions
    for r in bam.records() {
        let record = r?;
        let pos = record.pos() as u32;
        progress.set_position(pos as u64);

        if positions.contains_key(&pos) && record.mapq() >= min_quality {
            if let Some(base) = record.seq().as_bytes().first().map(|&b| b as char) {
                coverage.entry(pos).or_default().push(base);
            }
        }
    }

    // Second pass: analyze coverage
    for (pos, bases) in coverage {
        if bases.len() >= min_depth as usize {
            let mut base_counts: HashMap<char, u32> = HashMap::new();
            for base in bases {
                *base_counts.entry(base).or_insert(0) += 1;
            }

            if let Some((base, count)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                let freq = *count as f64 / base_counts.values().sum::<u32>() as f64;
                if freq >= 0.8 {
                    snp_calls.insert(pos, (*base, base_counts.values().sum(), freq));
                }
            }
        }
    }

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

fn validate_hg38_reference<R: Read>(bam: &R) -> Result<(), Box<dyn std::error::Error>> {
    let header = bam.header();
    let sq_count = header.target_count();

    if sq_count == 0 {
        return Err("BAM file has no reference sequences".into());
    }

    // Check if chrY or chrM exists based on need
    let has_chry = header.tid(b"chrY").is_some();
    let has_chrm = header.tid(b"chrM").is_some();

    if !has_chry && !has_chrm {
        return Err("BAM file missing required chromosome (chrY or chrM)".into());
    }

    Ok(())
}
