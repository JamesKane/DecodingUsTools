use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

macro_rules! debug_println {
    ($($arg:tt)*) => {
        #[cfg(debug_assertions)]
        println!($($arg)*);
    };
}

#[derive(Deserialize, Debug)]
pub struct Variant {
    variant: String,
    #[serde(rename = "position", default)]
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

#[derive(Deserialize, Debug)]
struct RawVariant {
    #[serde(default)]
    variant: Option<String>,
    #[serde(alias = "pos")]
    position: Option<i32>,
    #[serde(default)]
    ancestral: Option<String>,
    #[serde(default)]
    derived: Option<String>,
    #[serde(default)]
    region: Option<String>,
    #[serde(rename = "snpId", default)]
    snp_id: Option<u32>,
}

// Deal with incoherent variants in FTDNA's tree
impl TryFrom<RawVariant> for Variant {
    type Error = &'static str;

    fn try_from(raw: RawVariant) -> Result<Self, Self::Error> {
        Ok(Variant {
            variant: raw.variant.ok_or("missing variant")?,
            pos: raw.position.ok_or("missing position")?.unsigned_abs(),
            ancestral: raw.ancestral.ok_or("missing ancestral")?,
            derived: raw.derived.ok_or("missing derived")?,
            region: raw.region.unwrap_or_default(),
            snp_id: raw.snp_id.unwrap_or_default(),
        })
    }
}

impl Variant {
    fn to_snp(&self, tree_type: TreeType) -> Option<Snp> {
        // Only convert variants that have all required fields
        if self.variant.is_empty() || self.ancestral.is_empty() || self.derived.is_empty() {
            return None;
        }

        Some(Snp {
            position: self.pos,
            ancestral: self.ancestral.clone(),
            derived: self.derived.clone(),
            chromosome: match tree_type {
                TreeType::YDNA => "chrY".to_string(),
                TreeType::MTDNA => "chrM".to_string(),
            },
            build: match tree_type {
                TreeType::YDNA => "hg38".to_string(),
                TreeType::MTDNA => "rCRS".to_string(),
            },
        })
    }
}

#[derive(Deserialize, Debug)]
pub struct HaplogroupNode {
    #[serde(rename = "haplogroupId")]
    haplogroup_id: u32,
    #[serde(rename = "parentId", default)]
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
    #[serde(deserialize_with = "deserialize_nodes")]
    all_nodes: HashMap<String, HaplogroupNode>,
}

#[derive(Debug)]
pub struct Haplogroup {
    name: String,
    parent: Option<String>,
    snps: Vec<Snp>,
    children: Vec<Haplogroup>,
}

fn deserialize_nodes<'de, D>(deserializer: D) -> Result<HashMap<String, HaplogroupNode>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::Error;

    // First deserialize into a temporary structure with RawVariants
    #[derive(Deserialize)]
    struct TempNode {
        #[serde(rename = "haplogroupId")]
        haplogroup_id: u32,
        #[serde(rename = "parentId", default)]
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
        variants: Vec<RawVariant>,
        #[serde(default)]
        children: Vec<u32>,
    }

    let temp_map: HashMap<String, TempNode> = HashMap::deserialize(deserializer)?;

    // Convert TempNode to HaplogroupNode, processing variants
    temp_map
        .into_iter()
        .map(|(k, v)| {
            let variants = v
                .variants
                .into_iter()
                .filter_map(|raw| Variant::try_from(raw).ok())
                .collect();

            let node = HaplogroupNode {
                haplogroup_id: v.haplogroup_id,
                parent_id: v.parent_id,
                name: v.name,
                is_root: v.is_root,
                root: v.root,
                kits_count: v.kits_count,
                sub_branches: v.sub_branches,
                big_y_count: v.big_y_count,
                variants,
                children: v.children,
            };
            Ok((k, node))
        })
        .collect()
}

impl HaplogroupTree {
    fn build_tree(&self, node_id: u32, tree_type: TreeType) -> Option<Haplogroup> {
        let node = self.all_nodes.get(&node_id.to_string())?;

        let parent = if node.parent_id != 0 {
            self.all_nodes
                .get(&node.parent_id.to_string())
                .map(|p| p.name.clone())
        } else {
            None
        };

        let snps = node
            .variants
            .iter()
            .filter_map(|v| v.to_snp(tree_type))
            .collect();

        let children = node
            .children
            .iter()
            .filter_map(|&child_id| self.build_tree(child_id, tree_type))
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
    #[serde(skip)]
    chromosome: String,
    #[serde(skip)]
    build: String,
}

#[derive(Debug, Clone, Default, PartialEq, PartialOrd)]
struct HaplogroupScore {
    matches: u32,
    mismatches: u32,
    total_snps: u32,
    score: f64,
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

    let root_count = tree_json
        .all_nodes
        .values()
        .filter(|node| node.parent_id == 0)
        .count();
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
        .build_tree(root_node.haplogroup_id, tree_type)
        .ok_or("Failed to build tree")?;

    // Create map of positions to check
    let mut positions: HashMap<u32, Vec<(&str, &Snp)>> = HashMap::new();
    collect_snps(&tree, &mut positions);

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
    calculate_haplogroup_score(&tree, &snp_calls, &mut scores, None);

    // Sort and write results
    scores.sort_by(|a, b| {
        let score_diff = b.1.partial_cmp(&a.1).unwrap();
        if (b.1 - a.1).abs() < 0.1 {
            // If scores are very close, compare by number of matching SNPs
            b.2.cmp(&a.2)
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

    for r in bam.records() {
        let record = r?;
        let start_pos = record.pos() as u32;
        progress.set_position(start_pos as u64);

        if record.mapq() >= min_quality {
            let sequence = record.seq().as_bytes();
            for (read_pos, base) in sequence.iter().enumerate() {
                let ref_pos = start_pos + read_pos as u32;
                if positions.contains_key(&ref_pos) {
                    coverage.entry(ref_pos).or_default().push(*base as char);
                }
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
    debug_println!(
        "Processing haplogroup: {} with {} SNPs",
        haplogroup.name,
        haplogroup.snps.len()
    );
    for snp in &haplogroup.snps {
        debug_println!(
            "  SNP at pos {}: ancestral={}, derived={}, chrom={}, build={}",
            snp.position, snp.ancestral, snp.derived, snp.chromosome, snp.build
        );

        if is_valid_snp(snp) {
            positions
                .entry(snp.position)
                .or_default()
                .push((&haplogroup.name, snp));
        } else {
            debug_println!(
                "  Invalid SNP: pos={}, chrom={}, build={}",
                snp.position, snp.chromosome, snp.build
            );
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
    scores: &mut Vec<(String, f64, u32, u32)>,
    parent_score: Option<HaplogroupScore>,
) -> HaplogroupScore {
    // Start with parent scores if they exist, or create new score
    let mut current_score = parent_score.unwrap_or_default();

    // Add current haplogroup's SNPs to the score
    let valid_snps: Vec<_> = haplogroup.snps.iter().filter(|snp| is_valid_snp(snp)).collect();

    for snp in valid_snps.iter() {
        if let Some((called_base, _depth, _freq)) = snp_calls.get(&snp.position) {
            current_score.total_snps += 1;

            if *called_base == snp.derived.chars().next().unwrap() {
                current_score.matches += 1;
            } else {
                current_score.mismatches += 1;
            }
        }
    }

    // Calculate score
    let score = if current_score.total_snps > 0 {
        current_score.matches as f64 / current_score.total_snps as f64
    } else {
        0.0
    };
    current_score.score = score;

    // Store result in simple format expected by existing code
    scores.push((
        haplogroup.name.clone(),
        score,
        current_score.matches,
        current_score.total_snps,
    ));

    // Process children with current scores as their parent score
    for child in &haplogroup.children {
        calculate_haplogroup_score(child, snp_calls, scores, Some(current_score.clone()));
    }

    current_score
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
