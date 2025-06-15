use crate::utils::histogram_plotter::{self, CoverageRange};
use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read};
use std::collections::HashMap;
use std::fs::File;
use std::io::Result as IoResult;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

// TODO FIX-ME: This is still kinda slow.  About 2 minutes per contig in release mode
#[derive(Clone)]
pub struct CallableOptions {
    pub min_depth: u32,
    pub max_depth: u32,
    pub min_mapping_quality: u8,
    pub min_base_quality: u8,
    pub min_depth_for_low_mapq: u32,
    pub max_low_mapq: u8,
    pub max_low_mapq_fraction: f64,
}

impl CallableOptions {
    pub fn new(
        min_depth: u32,
        max_depth: u32,
        min_mapping_quality: u8,
        min_base_quality: u8,
        min_depth_for_low_mapq: u32,
        max_low_mapq: u8,
        max_low_mapq_fraction: f64,
    ) -> Self {
        Self {
            min_depth,
            max_depth,
            min_mapping_quality,
            min_base_quality,
            min_depth_for_low_mapq,
            max_low_mapq,
            max_low_mapq_fraction,
        }
    }
}

struct ContigStats {
    name: String,
    length: usize,
    n_covered_bases: u64,
    summed_coverage: u64,
    summed_baseq: u64,
    summed_mapq: u64,
    quality_bases: u64,
    n_reads: u32,
    n_selected_reads: u32,
    progress_bar: ProgressBar,
    options: CallableOptions,
}

impl ContigStats {
    fn new(
        name: String,
        length: usize,
        multi_progress: &MultiProgress,
        options: CallableOptions,
    ) -> Self {
        let progress_bar = multi_progress.add(ProgressBar::new(length as u64));
        progress_bar.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-"));
        progress_bar.set_message(format!("Processing {}", name));

        ContigStats {
            name,
            length,
            n_covered_bases: 0,
            summed_coverage: 0,
            summed_baseq: 0,
            summed_mapq: 0,
            quality_bases: 0,
            n_reads: 0,
            n_selected_reads: 0,
            progress_bar,
            options,
        }
    }

    fn process_position(
        &mut self,
        raw_depth: u32,
        alignments: bam::pileup::Alignments,
    ) {
        let mut qc_depth = 0;

        for aln in alignments {
            let record = aln.record();
            let mapq = record.mapq();

            // Count total reads
            self.n_reads += 1;

            // Check quality criteria
            if mapq >= self.options.min_mapping_quality {
                if let Some(qpos) = aln.qpos() {
                    if let Some(&qual) = record.qual().get(qpos) {
                        if qual >= self.options.min_base_quality {
                            qc_depth += 1;
                            self.summed_baseq += qual as u64;
                            self.quality_bases += 1;
                        }
                    }
                }
                self.summed_mapq += mapq as u64;
                self.n_selected_reads += 1;
            }
        }

        if raw_depth > 0 {
            self.n_covered_bases += 1;
            self.summed_coverage += raw_depth as u64;
        }
    }
}

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledState {
    /// The reference base was an N, which is not considered callable
    REF_N,
    /// The base satisfied the minimum depth for calling but had less than maxDepth
    /// to avoid having EXCESSIVE_COVERAGE
    CALLABLE,
    /// Absolutely no reads were seen at this locus, regardless of filtering parameters
    NO_COVERAGE,
    /// There were fewer than minimum depth bases at the locus, after applying filters
    LOW_COVERAGE,
    /// More than maxDepth reads at the locus, indicating some sort of mapping problem
    EXCESSIVE_COVERAGE,
    /// More than maxFractionOfReadsWithLowMAPQ at the locus, indicating poor mapping quality
    POOR_MAPPING_QUALITY,
}

impl CalledState {
    pub fn from_usize(val: usize) -> Option<Self> {
        match val {
            0 => Some(Self::REF_N),
            1 => Some(Self::CALLABLE),
            2 => Some(Self::NO_COVERAGE),
            3 => Some(Self::LOW_COVERAGE),
            4 => Some(Self::EXCESSIVE_COVERAGE),
            5 => Some(Self::POOR_MAPPING_QUALITY),
            _ => None,
        }
    }
}

pub struct CallableLociCounter {
    counts: [u64; 6],
    current_state: Option<(String, u64, u64, CalledState)>,
    bed_writer: BufWriter<File>,
    coverage_ranges: Vec<CoverageRange>,
    current_pos: u32,
    output_dir: PathBuf,
    largest_contig_length: u32,
    contig_counts: HashMap<String, [u64; 6]>,
}

impl CallableLociCounter {
    pub fn new(bed_file: &str, summary_file: &str, largest_contig_length: u32) -> IoResult<Self> {
        let output_dir = PathBuf::from(bed_file)
            .parent()
            .unwrap_or(&PathBuf::from("."))
            .to_path_buf();

        Ok(Self {
            counts: [0; 6],
            current_state: None,
            bed_writer: BufWriter::new(File::create(bed_file)?),
            coverage_ranges: Vec::new(),
            current_pos: 0,
            output_dir,
            largest_contig_length,
            contig_counts: Default::default(),
        })
    }

    fn write_state(&mut self) -> IoResult<()> {
        if let Some((contig, start, end, state)) = &self.current_state {
            writeln!(
                self.bed_writer,
                "{}\t{}\t{}\t{:?}",
                contig, start, end, state
            )?;
        }
        Ok(())
    }

    fn finish_contig(&mut self, contig_name: &str, contig_length: u32) -> IoResult<()> {
        self.write_state()?;

        if !self.coverage_ranges.is_empty() {
            let histogram_path = histogram_plotter::generate_histogram(
                std::mem::take(&mut self.coverage_ranges),
                &format!("{}_coverage", contig_name),
                contig_name,
                contig_length,
                if contig_name == "chrM" { contig_length } else { self.largest_contig_length },
            )?;

            let final_path = self
                .output_dir
                .join(format!("{}_coverage.svg", contig_name));
            std::fs::copy(histogram_path, final_path)?;
        }

        Ok(())
    }

    fn process_position(
        &mut self,
        contig: &str,
        pos: u32,
        depth: u32,
        is_low_mapq: bool,
        state: CalledState,
    ) -> IoResult<()> {
        self.counts[state as usize] += 1;
        // Update per-contig counts
        self.contig_counts
            .entry(contig.to_string())
            .or_insert([0; 6])[state as usize] += 1;


        match self.coverage_ranges.last_mut() {
            Some(range) if range.can_merge(pos, depth, is_low_mapq) => {
                range.extend(pos);
            }
            _ => {
                self.coverage_ranges.push(CoverageRange::new(pos, depth, is_low_mapq));
            }
        }
        
        // Handle state tracking (existing code)
        match &mut self.current_state {
            Some((cur_contig, start, end, cur_state)) => {
                if contig != cur_contig || state != *cur_state || pos as u64 != *end + 1 {
                    self.write_state()?;
                    self.current_state = Some((contig.to_string(), pos as u64, pos as u64, state));
                } else {
                    *end = pos as u64;
                }
            }
            None => {
                self.current_state = Some((contig.to_string(), pos as u64, pos as u64, state));
            }
        }

        Ok(())
    }

    pub fn get_contig_counts(&self, contig: &str) -> [u64; 6] {
        self.contig_counts.get(contig).copied().unwrap_or([0; 6])
    }
}

// Add a helper function for natural sorting of chromosome names
fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    // Special handling for chrM to always put it at the end
    if a == "chrM" {
        return std::cmp::Ordering::Greater;
    }
    if b == "chrM" {
        return std::cmp::Ordering::Less;
    }

    // Special handling for chrX and chrY
    if a == "chrX" && b == "chrY" {
        return std::cmp::Ordering::Less;
    }
    if a == "chrY" && b == "chrX" {
        return std::cmp::Ordering::Greater;
    }
    if a == "chrX" || a == "chrY" {
        return std::cmp::Ordering::Greater;
    }
    if b == "chrX" || b == "chrY" {
        return std::cmp::Ordering::Less;
    }

    // For other chromosomes, split into prefix and number
    let (a_prefix, a_num) = a.split_at(a.chars().take_while(|c| !c.is_ascii_digit()).count());
    let (b_prefix, b_num) = b.split_at(b.chars().take_while(|c| !c.is_ascii_digit()).count());

    // Compare prefixes first
    match a_prefix.cmp(b_prefix) {
        std::cmp::Ordering::Equal => {
            // If prefixes are equal, compare numbers
            let a_num = a_num.parse::<u32>().unwrap_or(0);
            let b_num = b_num.parse::<u32>().unwrap_or(0);
            a_num.cmp(&b_num)
        }
        other => other,
    }
}

/// This function performs a comprehensive analysis of a BAM file to determine callable loci and 
/// generate corresponding statistics. The results are written to both BED and summary output files.
///
/// # Parameters
///
/// - `bam_file`: A `String` representing the path to the input BAM file. The BAM file should be indexed.
/// - `reference_file`: A `String` representing the path to the reference genome file (FASTA format). The FASTA file should be indexed.
/// - `output_bed`: A `String` representing the path to the output BED file where callable loci data will be written.
/// - `output_summary`: A `String` representing the path to the output summary file where contig-level statistics will be stored.
/// - `options`: A `CallableOptions` structure containing various threshold parameters and settings for processing.
///
/// # Returns
///
/// - `Ok(())`: If the function successfully completes the analysis and writes all outputs successfully.
/// - `Err(Box<dyn std::error::Error>)`: If any error occurs during the processing, such as file I/O or data parsing issues.
///
/// # Behavior
///
/// 1. Reads the BAM and reference genome files.
/// 2. Initializes progress bars for tracking coverage analysis and handles up to 25 contigs.
/// 3. Configures parameters for pileup based on the options provided.
/// 4. Iterates over pileup data to compute depth, mapping quality, base quality, and categorizes loci into callable or non-callable states.
/// 5. Updates progress bars in real-time for processing transparency.
/// 6. Writes callable loci results to the BED file.
/// 7. Generates a contig summary in the summary file with the following statistics:
///    - Contig name and length
///    - Count of unique reads
///    - Count of loci marked as reference 'N'
///    - Counts of loci with no coverage, low coverage, excessive coverage, and poor mapping quality
///    - Callable loci count
///    - Coverage percentage, average depth, average mapping quality, and average base quality.
/// 8. Handles file I/O errors gracefully and ensures all progress bars finish properly before exiting.
///
/// # Progress Bars
///
/// - A main spinner tracks overall task progression.
/// - Individual progress bars track per-contig processing progress.
///
/// # Callable Loci State
///
/// For each position in the sequence, the following states are determined:
/// - `REF_N`: Reference base is 'N' or 'n'.
/// - `NO_COVERAGE`: No reads map to the position.
/// - `LOW_COVERAGE`: Insufficient reads with quality mapping to the position.
/// - `POOR_MAPPING_QUALITY`: High fraction of low-quality reads.
/// - `EXCESSIVE_COVERAGE`: More reads than the maximum depth threshold.
/// - `CALLABLE`: Meets all quality and depth criteria for being callable.
///
/// # Example Usage
///
/// ```
/// use my_crate::callable_loci::run;
/// use my_crate::callable_loci::CallableOptions;
///
/// let result = run(
///     "input.bam".to_string(),
///     "
pub fn run(
    bam_file: String,
    reference_file: String,
    output_bed: String,
    output_summary: String,
    options: CallableOptions,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut bam = bam::IndexedReader::from_path(&bam_file)?;
    let header = bam.header().clone();
    let mut fasta = IndexedReader::from_file(&reference_file)?;

    // Create progress bars
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    main_progress.set_message("Analyzing coverage...");

    // Create contig progress bars map using usize as key
    let mut contig_stats = HashMap::new();

    // First pass to get contig lengths and set up progress bars
    for tid in 0..header.target_names().len().min(25) {
        let contig_name = std::str::from_utf8(header.tid2name(tid as u32))?;
        let length = header.target_len(tid as u32).unwrap_or(0) as usize;
        contig_stats.insert(
            tid as usize, // Convert tid to usize
            ContigStats::new(
                contig_name.to_string(),
                length,
                &multi_progress,
                options.clone(),
            ),
        );
    }

    // Find largest non-chrM contig length
    let largest_contig_length = contig_stats.values()
        .filter(|stats| stats.name != "chrM")
        .map(|stats| stats.length)
        .max()
        .unwrap_or(0) as u32;

    let mut counter = CallableLociCounter::new(&output_bed, &output_summary, largest_contig_length)?;

    let mut pileup = bam.pileup();
    pileup.set_max_depth(if options.max_depth > 0 {
        options.max_depth
    } else {
        500
    });

    let mut current_contig = String::new();
    // Process pileup
    for p in pileup {
        let pileup = p?;
        let tid = pileup.tid() as usize; // Convert tid to usize here
        let pos = pileup.pos();
        let contig = std::str::from_utf8(header.tid2name(pileup.tid()))?;
        
        // Detect contig transition
        if current_contig != contig && !current_contig.is_empty() {
            // Use the length from contig_stats for the current contig
            let length = contig_stats.values()
                .find(|stats| stats.name == current_contig)
                .map(|stats| stats.length)
                .unwrap_or(0) as u32;
            counter.finish_contig(&current_contig, length)?;
        }
        current_contig = contig.to_string();


        // Update progress for current contig
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.progress_bar.set_position(pos as u64);
        }

        // Get reference base
        let mut seq = Vec::new();
        fasta.fetch(contig, pos as u64, (pos + 1) as u64)?;
        fasta.read(&mut seq)?;
        let ref_base = seq.first().map(|&b| b).unwrap_or(b'N');

        let mut raw_depth = 0;
        let mut qc_depth = 0;
        let mut low_mapq_count = 0;

        // Process alignments
        for align in pileup.alignments() {
            raw_depth += 1;
            let record = align.record();

            if record.mapq() <= options.max_low_mapq {
                low_mapq_count += 1;
            }

            if record.mapq() >= options.min_mapping_quality {
                if let Some(qpos) = align.qpos() {
                    if let Some(&qual) = record.qual().get(qpos) {
                        if qual >= options.min_base_quality || align.is_del() {
                            qc_depth += 1;
                        }
                    }
                }
            }
        }
        // First determine the state based on the parameters
        let state = if ref_base == b'N' || ref_base == b'n' {
            CalledState::REF_N
        } else if raw_depth == 0 {
            CalledState::NO_COVERAGE
        } else if raw_depth >= options.min_depth_for_low_mapq
            && (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction
        {
            CalledState::POOR_MAPPING_QUALITY
        } else if qc_depth < options.min_depth {
            CalledState::LOW_COVERAGE
        } else if options.max_depth > 0 && qc_depth > options.max_depth {
            CalledState::EXCESSIVE_COVERAGE
        } else {
            CalledState::CALLABLE
        };

        let is_low_mapq = raw_depth >= options.min_depth_for_low_mapq
            && (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction;

        // Now call process_position with the correct parameters
        counter.process_position(
            contig,
            pos,
            qc_depth,
            is_low_mapq,
            state,
        )?;

        // Update contig stats
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.process_position(raw_depth, pileup.alignments());
        }
    }

    // Finish all progress bars
    for (_, stats) in contig_stats.iter() {
        stats.progress_bar.finish_with_message("Complete");
    }
    main_progress.finish_with_message("Coverage analysis complete!");

    // Finish final contig
    let final_length = contig_stats.values()
        .find(|stats| stats.name == current_contig)
        .map(|stats| stats.length)
        .unwrap_or(0) as u32;
    counter.finish_contig(&current_contig, final_length)?;

    let mut summary_writer = BufWriter::new(File::create(output_summary)?);

    // Write header
    writeln!(
        summary_writer,
        "Contig|Start|Stop|UniqueReads|RefN|NoCoverage|LowCoverage|ExcessiveCoverage|PoorMappingQuality|Callable|CoveragePercent|AvgDepth|AvgMapQ|AvgBaseQ"
    )?;

    // Write stats for each contig
    let mut sorted_stats: Vec<_> = contig_stats.iter().collect();
    sorted_stats.sort_by(|a, b| natural_cmp(&a.1.name, &b.1.name));
    for (_, stats) in sorted_stats {
        let counts = counter.get_contig_counts(&stats.name);
        writeln!(
            summary_writer,
            "{}|1|{}|{}|{}|{}|{}|{}|{}|{}|{:.2}|{:.2}|{:.1}|{:.1}",
            stats.name,
            stats.length,
            stats.n_selected_reads,
            counts[CalledState::REF_N as usize],
            counts[CalledState::NO_COVERAGE as usize],
            counts[CalledState::LOW_COVERAGE as usize],
            counts[CalledState::EXCESSIVE_COVERAGE as usize],
            counts[CalledState::POOR_MAPPING_QUALITY as usize],
            counts[CalledState::CALLABLE as usize],
            (stats.n_covered_bases as f64 / stats.length as f64) * 100.0,
            stats.summed_coverage as f64 / stats.length as f64,
            if stats.n_selected_reads > 0 {
                stats.summed_mapq as f64 / stats.n_selected_reads as f64
            } else {
                0.0
            },
            if stats.quality_bases > 0 {
                stats.summed_baseq as f64 / stats.quality_bases as f64
            } else {
                0.0
            }
        )?;
    }

    Ok(())
}
