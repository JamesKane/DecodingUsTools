mod scoring;
pub(crate) mod types;
mod validation;
mod tree;

use crate::haplogroup::types::{Haplogroup, HaplogroupResult, Snp};
use crate::utils::cache::{TreeCache, TreeType};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam::IndexedReader, bam::Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};

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
    validation::validate_hg38_reference(&bam)?;

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
    tree::collect_snps(&haplogroup_tree, &mut positions);

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
                        ref_pos += op.len();
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
