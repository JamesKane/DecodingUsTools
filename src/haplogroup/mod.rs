use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(Deserialize, Debug)]
pub struct Variant {
    pub variant: String,
    #[serde(rename = "position", default)]
    pub pos: u32,
    #[serde(default)]
    pub ancestral: String,
    #[serde(default)]
    pub derived: String,
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupNode {
    pub haplogroup_id: u32,
    pub parent_id: u32,
    pub name: String,
    pub is_root: bool,
    pub variants: Vec<Variant>,
    pub children: Vec<u32>,
}

#[derive(Deserialize)]
pub struct HaplogroupTree {
    #[serde(rename = "allNodes")]
    pub all_nodes: HashMap<String, HaplogroupNode>,
}

#[derive(Debug)]
pub struct Haplogroup {
    pub(crate) name: String,
    pub(crate) parent: Option<String>,
    pub(crate) snps: Vec<Snp>,
    pub(crate) children: Vec<Haplogroup>,
}

#[derive(Deserialize, Debug, Clone)]
pub struct Snp {
    pub(crate) position: u32,
    pub(crate) ancestral: String,
    pub(crate) derived: String,
    #[serde(skip)]
    pub(crate) chromosome: String,
    #[serde(skip)]
    pub(crate) build: String,
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
struct HaplogroupScore {
    matches: usize,
    ancestral_matches: usize,
    no_calls: usize,
    total_snps: usize,
    score: f64,
    depth: usize, // Added depth field
}

impl Default for HaplogroupScore {
    fn default() -> Self {
        Self {
            matches: 0,
            ancestral_matches: 0,
            no_calls: 0,
            total_snps: 0,
            score: 0.0,
            depth: 0,
        }
    }
}

#[derive(Debug, Clone)]
struct HaplogroupResult {
    name: String,
    score: f64,
    matching_snps: u32,
    mismatching_snps: u32,
    ancestral_matches: u32,
    no_calls: u32,
    total_snps: u32,
    cumulative_snps: u32, // Total unique SNPs from root to this branch
    depth: u32,
}

pub fn analyze_haplogroup(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    min_quality: u8,
    tree_type: TreeType,
) -> Result<(), Box<dyn std::error::Error>> {
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
    let mut bam = IndexedReader::from_path(&bam_file)?;
    validate_hg38_reference(&bam)?;

    // Get tree from cache
    let tree_cache = TreeCache::new(tree_type)?;
    let tree = tree_cache.get_tree()?;

    let root_count = tree
        .all_nodes
        .values()
        .filter(|node| node.parent_id == 0)
        .count();
    if root_count > 1 {
        return Err("Multiple root nodes found in tree".into());
    }

    let root_node = tree
        .all_nodes
        .values()
        .find(|node| node.parent_id == 0)
        .ok_or("No root node found")?;

    let haplogroup_tree = tree_cache
        .provider
        .build_tree(&tree, root_node.haplogroup_id, tree_type)
        .ok_or("Failed to build tree")?;

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Snp)>> = HashMap::new();
    collect_snps(&haplogroup_tree, &mut positions);

    // Determine which chromosome we need to process
    let need_y = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "chrY"));
    let need_mt = positions
        .iter()
        .any(|(_, snps)| snps.iter().any(|(_, snp)| snp.chromosome == "chrM"));

    let mut snp_calls: HashMap<u32, (char, u32, f64)> = HashMap::new();

    // Create a progress bar for BAM processing
    let progress_style = ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
        .unwrap()
        .progress_chars("#>-");

    let header = bam.header().clone();

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
            &header,
            &mut fasta_reader,
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
            &header,
            &mut fasta_reader,
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
    calculate_haplogroup_score(&haplogroup_tree, &snp_calls, &mut scores, None, 0);

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

fn process_region<R: Read>(
    bam: &mut R,
    header: &rust_htslib::bam::HeaderView,
    fasta: &mut bio::io::fasta::IndexedReader<File>,
    positions: &HashMap<u32, Vec<(&str, &Snp)>>,
    min_depth: u32,
    min_quality: u8,
    snp_calls: &mut HashMap<u32, (char, u32, f64)>,
    progress: &ProgressBar,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut coverage: HashMap<u32, Vec<u8>> = HashMap::new();

    for r in bam.records() {
        let record = r?;
        let start_pos = record.pos() as u32;
        progress.set_position(start_pos as u64);

        if record.mapq() >= min_quality {
            let sequence = record.seq().as_bytes();
            let cigar = record.cigar();
            let mut ref_pos = start_pos;
            let mut read_pos = 0;

            let tid = record.tid();
            let ref_name = String::from_utf8_lossy(header.tid2name(tid as u32));

            for op in cigar.iter() {
                match op.char() {
                    'M' | '=' | 'X' => {
                        let len = op.len() as usize;
                        for i in 0..len {
                            let vcf_pos = ref_pos + 1;
                            if positions.contains_key(&vcf_pos) && read_pos + i < sequence.len() {
                                let base_index = read_pos + i;
                                let base = sequence[base_index].to_ascii_uppercase();

                                if let Ok(_) =
                                    fasta.fetch(&ref_name, ref_pos as u64, (ref_pos + 1) as u64)
                                {
                                    let mut ref_seq = Vec::new();
                                    fasta.read(&mut ref_seq)?;
                                    coverage.entry(vcf_pos).or_default().push(base);
                                }
                            }
                            ref_pos += 1;
                        }
                        read_pos += len;
                    }
                    'D' | 'N' => {
                        ref_pos += op.len() as u32;
                    }
                    'I' | 'S' => {
                        read_pos += op.len() as usize;
                    }
                    _ => {}
                }
            }
        }
    }

    for (pos, bases) in coverage {
        if bases.len() >= min_depth as usize {
            let mut base_counts: HashMap<char, u32> = HashMap::new();

            for &base in &bases {
                *base_counts.entry(base as char).or_insert(0) += 1;
            }

            if let Some((&base, &count)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                let total = bases.len() as u32;
                let freq = count as f64 / total as f64;

                if freq >= 0.7 {
                    snp_calls.insert(pos, (base, total, freq));
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
        "chrM" => true,
        "chrY" => snp.build == "hg38",
        _ => false,
    }
}

fn calculate_haplogroup_score(
    haplogroup: &Haplogroup,
    snp_calls: &HashMap<u32, (char, u32, f64)>,
    scores: &mut Vec<HaplogroupResult>,
    parent_info: Option<&(HaplogroupScore, HashSet<u32>)>,
    depth: u32,
) -> (HaplogroupScore, HashSet<u32>) {
    let mut current_score = HaplogroupScore::default();
    let mut cumulative_snps = parent_info
        .map(|(_, snps)| snps.clone())
        .unwrap_or_default();

    for snp in &haplogroup.snps {
        if is_valid_snp(snp) {
            cumulative_snps.insert(snp.position);
        }
    }

    let defining_snps: Vec<_> = haplogroup
        .snps
        .iter()
        .filter(|snp| is_valid_snp(snp))
        .collect();

    let mut branch_derived = 0;
    let mut branch_ancestral = 0;
    let mut branch_no_calls = 0;
    let mut branch_low_quality = 0;

    // Process defining SNPs
    for snp in &defining_snps {
        if let Some((called_base, depth, freq)) = snp_calls.get(&snp.position) {
            if *depth >= 4 {
                let derived_base = snp.derived.chars().next().unwrap();
                let ancestral_base = snp.ancestral.chars().next().unwrap();

                if *called_base == derived_base {
                    if *freq >= 0.7 {
                        branch_derived += 1;
                    } else if *freq >= 0.5 {
                        branch_derived += 1;
                    } else {
                        branch_low_quality += 1;
                    }
                } else if *called_base == ancestral_base {
                    if *freq >= 0.7 {
                        branch_ancestral += 1;
                    } else {
                        branch_low_quality += 1;
                    }
                } else if *freq >= 0.7 {
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

    // Calculate score with improved logic
    let total_calls = branch_derived + branch_ancestral + branch_low_quality;
    if total_calls > 0 {
        // Revised scoring logic
        let branch_score = if branch_ancestral == 0 {
            // No contradicting evidence
            match branch_derived {
                d if d >= 1 => 3.08, // Any clean derived matches with no ancestral contradictions
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
    current_score.total_snps += defining_snps.len();

    if branch_ancestral > branch_derived * 10 {
        scores.push(HaplogroupResult {
            name: haplogroup.name.clone(),
            score: 0.0,
            matching_snps: branch_derived.try_into().unwrap_or(0),
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: branch_ancestral.try_into().unwrap_or(0),
            no_calls: branch_no_calls.try_into().unwrap_or(0),
            total_snps: defining_snps.len() as u32,
            cumulative_snps: cumulative_snps.len() as u32,
            depth,
        });
        return (current_score, cumulative_snps);
    }

    for child in &haplogroup.children {
        let (child_score, child_cumulative) = calculate_haplogroup_score(
            child,
            snp_calls,
            scores,
            Some(&(current_score.clone(), cumulative_snps.clone())),
            depth + 1,
        );

        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: child_score.score,
            matching_snps: child_score.matches.try_into().unwrap_or(0),
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: child_score.ancestral_matches.try_into().unwrap_or(0),
            no_calls: child_score.no_calls.try_into().unwrap_or(0),
            total_snps: defining_snps.len() as u32,
            cumulative_snps: child_cumulative.len() as u32,
            depth,
        });
    }

    (current_score, cumulative_snps)
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

fn find_path_to_root(
    haplogroup: &Haplogroup,
    target_name: &str,
    scores: &[HaplogroupResult],
) -> Option<Vec<String>> {
    if haplogroup.name == target_name {
        return Some(vec![haplogroup.name.clone()]);
    }

    for child in &haplogroup.children {
        if let Some(mut path) = find_path_to_root(child, target_name, scores) {
            path.push(haplogroup.name.clone());
            return Some(path);
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
        b.cumulative_snps
            .cmp(&a.cumulative_snps)
            .then_with(|| {
                b.score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    // Take the top result and ensure its ancestral path is included first
    if let Some(top_result) = remaining.first().cloned() {
        if let Some(path) = find_path_to_root(tree, &top_result.name, &remaining) {
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
        b.cumulative_snps
            .cmp(&a.cumulative_snps)
            .then_with(|| {
                b.score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    // Add remaining results
    ordered_results.extend(remaining);

    ordered_results
}
