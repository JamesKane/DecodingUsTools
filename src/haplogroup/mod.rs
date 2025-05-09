use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
use serde::Deserialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};

const BAM_CMATCH: u32 = 0;
const BAM_CINS: u32 = 1;
const BAM_CDEL: u32 = 2;
const BAM_CREF_SKIP: u32 = 3;
const BAM_CSOFT_CLIP: u32 = 4;
const BAM_CEQUAL: u32 = 7;
const BAM_CDIFF: u32 = 8;

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

#[derive(Deserialize, Debug, Clone)]
pub struct Snp {
    position: u32,
    ancestral: String,
    derived: String,
    #[serde(skip)]
    chromosome: String,
    #[serde(skip)]
    build: String,
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

#[derive(Debug)]
struct HaplogroupResult {
    name: String,
    score: f64,
    matching_snps: u32,
    mismatching_snps: u32,
    ancestral_matches: u32,
    no_calls: u32,
    total_snps: u32,
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
    calculate_haplogroup_score(&tree, &snp_calls, &mut scores, None, 0);

    // Sort and write results
    scores.sort_by(|a, b| {
        // First compare scores with NaN handling
        match (a.score.is_nan(), b.score.is_nan()) {
            (true, true) => std::cmp::Ordering::Equal,
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            (false, false) => {
                match b
                    .score
                    .partial_cmp(&a.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
                {
                    std::cmp::Ordering::Equal => {
                        // If scores are equal, compare by matching SNPs
                        b.matching_snps.cmp(&a.matching_snps)
                    }
                    ord => ord,
                }
            }
        }
    });

    // Update the report writing section in analyze_haplogroup:
    let mut writer = BufWriter::new(File::create(output_file)?);
    writeln!(
        writer,
        "Haplogroup\tScore\tMatching_SNPs\tMismatching_SNPs\tAncestral_Matches\tNo_Calls\tTotal_SNPs"
    )?;
    for result in scores {
        writeln!(
            writer,
            "{}\t{:.4}\t{}\t{}\t{}\t{}\t{}",
            result.name,
            result.score,
            result.matching_snps,
            result.mismatching_snps,
            result.ancestral_matches,
            result.no_calls,
            result.total_snps
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
    debug_println!(
        "Processing haplogroup: {} with {} SNPs",
        haplogroup.name,
        haplogroup.snps.len()
    );
    for snp in &haplogroup.snps {
        debug_println!(
            "  SNP at pos {}: ancestral={}, derived={}, chrom={}, build={}",
            snp.position,
            snp.ancestral,
            snp.derived,
            snp.chromosome,
            snp.build
        );

        if is_valid_snp(snp) {
            positions
                .entry(snp.position)
                .or_default()
                .push((&haplogroup.name, snp));
        } else {
            debug_println!(
                "  Invalid SNP: pos={}, chrom={}, build={}",
                snp.position,
                snp.chromosome,
                snp.build
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
    scores: &mut Vec<HaplogroupResult>,
    parent_score: Option<HaplogroupScore>,
    consecutive_all_negative: u32,
) -> HaplogroupScore {
    let mut current_score = HaplogroupScore::default();
    current_score.depth = parent_score.map_or(0, |p| p.depth + 1);

    let defining_snps: Vec<_> = haplogroup
        .snps
        .iter()
        .filter(|snp| is_valid_snp(snp))
        .collect();

    let mut branch_derived = 0;
    let mut branch_ancestral = 0;
    let mut branch_no_calls = 0;
    let mut branch_low_quality = 0;

    let is_target = haplogroup.name.contains("FGC29071") || haplogroup.name.contains("FGC29067");

    if is_target {
        println!("\nAnalyzing target branch: {}", haplogroup.name);
    }

    // Process defining SNPs
    for snp in &defining_snps {
        if let Some((called_base, depth, freq)) = snp_calls.get(&snp.position) {
            if *depth >= 4 {
                let derived_base = snp.derived.chars().next().unwrap();
                let ancestral_base = snp.ancestral.chars().next().unwrap();

                if is_target {
                    println!(
                        "SNP at {}: called={} (depth={}, freq={:.2}) derived={} ancestral={}",
                        snp.position, called_base, depth, freq, derived_base, ancestral_base
                    );
                }

                if *called_base == derived_base {
                    if *freq >= 0.7 {
                        branch_derived += 1;
                        if is_target {
                            println!("  → Derived match (high confidence)");
                        }
                    } else if *freq >= 0.5 {
                        branch_derived += 1;
                        if is_target {
                            println!("  → Derived match (medium confidence)");
                        }
                    } else {
                        branch_low_quality += 1;
                        if is_target {
                            println!("  → Low quality call");
                        }
                    }
                } else if *called_base == ancestral_base {
                    if *freq >= 0.7 {
                        branch_ancestral += 1;
                        if is_target {
                            println!("  → Ancestral match (high confidence)");
                        }
                    } else {
                        branch_low_quality += 1;
                        if is_target {
                            println!("  → Low quality ancestral");
                        }
                    }
                } else if *freq >= 0.7 {
                    // Different mutation but high quality
                    branch_derived += 1;
                    if is_target {
                        println!("  → Different mutation (counting as derived)");
                    }
                } else {
                    branch_low_quality += 1;
                    if is_target {
                        println!("  → Low quality different mutation");
                    }
                }
            } else {
                branch_no_calls += 1;
                if is_target {
                    println!("  → Insufficient depth (< 4)");
                }
            }
        } else {
            branch_no_calls += 1;
            if is_target {
                println!("  → No call available");
            }
        }
    }

    // Calculate score with improved logic
    let total_calls = branch_derived + branch_ancestral + branch_low_quality;
    if total_calls > 0 {
        let derived_ratio = branch_derived as f64 / total_calls as f64;

        // Calculate quality metrics
        let coverage_ratio = total_calls as f64 / defining_snps.len() as f64;
        let quality_ratio = (branch_derived + branch_ancestral) as f64
            / (branch_derived + branch_ancestral + branch_low_quality) as f64;

        // Base score calculation with more weight on derived SNPs
        let derived_weight = match branch_derived {
            d if d >= 3 => 1.0, // Full weight for 3+ derived SNPs
            2 => 0.9,           // Strong weight for 2 derived SNPs
            1 => 0.8,           // Moderate weight for 1 derived SNP
            _ => 0.6,           // Low weight for no derived SNPs
        };

        // Penalize less for ancestral SNPs if we have good derived matches
        let ancestral_penalty = if branch_derived >= 3 {
            0.9 // Minor penalty when we have strong derived evidence
        } else if branch_derived >= 2 {
            0.8 // Moderate penalty
        } else {
            0.7 // Stronger penalty when we have weak derived evidence
        };

        // Calculate score for terminal vs non-terminal patterns
        let mut branch_score = 1.0;

        // Terminal branch scoring - reward the pattern of both derived and ancestral SNPs
        if branch_derived >= 2 && branch_ancestral > 0 {
            // This is a likely terminal branch pattern
            branch_score = match (branch_derived, branch_ancestral) {
                (d, a) if d >= 3 && a <= 2 => 2.8, // Ideal terminal pattern
                (d, a) if d >= 2 && a <= 2 => 2.5, // Good terminal pattern
                (d, _) if d >= 2 => 2.2,           // Possible terminal pattern
                _ => 2.0,
            };
        } else if branch_derived > 0 {
            // Non-terminal branch scoring
            branch_score = match branch_derived {
                d if d > 100 => 1.0, // Many SNPs - likely ancient branch
                d if d > 50 => 1.2,  // Too many SNPs for terminal branch
                d if d > 20 => 1.4,  // Getting closer to terminal branch size
                d if d > 10 => 1.6,  // Good number of SNPs
                d if d >= 3 => 1.8,  // Few but confident SNPs
                _ => 1.5,            // Too few SNPs to be confident
            };
        }

        // Adjust for quality and coverage
        let quality_factor = if branch_low_quality == 0 { 1.1 } else { 0.9 };
        let coverage_factor = match total_calls as f64 / defining_snps.len() as f64 {
            r if r >= 0.8 => 1.1,
            r if r >= 0.6 => 1.0,
            _ => 0.9,
        };

        current_score.score = branch_score * quality_factor * coverage_factor;

        if is_target {
            println!("\nScoring details for {}:", haplogroup.name);
            println!("  Derived SNPs: {}", branch_derived);
            println!("  Ancestral SNPs: {}", branch_ancestral);
            println!("  Low quality: {}", branch_low_quality);
            println!("  No calls: {}", branch_no_calls);
            println!("  Quality ratio: {:.2}", quality_ratio);
            println!("  Coverage ratio: {:.2}", coverage_ratio);
            println!("  Derived weight: {:.2}", derived_weight);
            println!(
                "  Ancestral penalty: {:.2}",
                if branch_ancestral > 0 {
                    ancestral_penalty
                } else {
                    1.0
                }
            );
            println!("  Base score: {:.2}", branch_score);
            println!("  Final score: {:.2}", current_score.score);
        }
    }

    // Update running totals
    current_score.matches += branch_derived;
    current_score.ancestral_matches += branch_ancestral;
    current_score.no_calls += branch_no_calls;
    current_score.total_snps += defining_snps.len();

    // Process children
    for child in &haplogroup.children {
        let child_score = calculate_haplogroup_score(
            child,
            snp_calls,
            scores,
            Some(current_score.clone()),
            consecutive_all_negative,
        );

        scores.push(HaplogroupResult {
            name: child.name.clone(),
            score: child_score.score,
            matching_snps: child_score.matches.try_into().unwrap_or(0),
            mismatching_snps: branch_low_quality.try_into().unwrap_or(0),
            ancestral_matches: child_score.ancestral_matches.try_into().unwrap_or(0),
            no_calls: child_score.no_calls.try_into().unwrap_or(0),
            total_snps: child_score.total_snps.try_into().unwrap_or(0),
        });
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
