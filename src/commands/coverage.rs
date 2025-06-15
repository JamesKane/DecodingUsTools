use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::htslib::{BAM_FSECONDARY, BAM_FSUPPLEMENTARY, BAM_FREAD1};
use rust_htslib::{bam, bam::Read};
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::utils::histogram_plotter::{self, CoverageRange};
use std::path::PathBuf;
use std::io::{self, Result as IoResult};


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


#[derive(Default)]
struct CallableStats {
    ref_n: usize,
    no_coverage: usize,
    low_coverage: usize,
    excessive_coverage: usize,
    poor_mapping_quality: usize,
    callable: usize,
}

struct ContigStats {
    name: String,
    length: usize,
    total_depth: u64,
    mapq_sum: u64,
    mapq_count: u64,
    total_base_coverage: u64,
    unique_read_count: u64,
    seen_read_names: HashSet<Vec<u8>>,
    stats: CallableStats,
    progress_bar: ProgressBar,
    options: CallableOptions,
    current_position: u32,
    coverage_ranges: Vec<CoverageRange>,
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
            total_depth: 0,
            mapq_sum: 0,
            mapq_count: 0,
            total_base_coverage: 0,
            unique_read_count: 0,
            seen_read_names: HashSet::new(),
            stats: CallableStats::default(),
            progress_bar,
            options,
            current_position: 0,
            coverage_ranges: Vec::new(),
        }
    }

    fn process_position(&mut self, pileup: &rust_htslib::bam::pileup::Pileup, depth: u32, alignments: bam::pileup::Alignments, ref_base: u8) {
        // Get the position from the pileup alignment
        self.current_position = pileup.pos() as u32;

        self.total_depth += depth as u64;

        // Fast-path checks first
        if ref_base == b'N' || ref_base == b'n' {
            self.stats.ref_n += 1;
            return;
        }
        if depth == 0 {
            self.stats.no_coverage += 1;
            return;
        }

        let mut primary_count = 0;
        let mut qc_depth = 0;
        let mut low_mapq_count = 0;
        let mut mapq_sum: u32 = 0;

        for aln in alignments {
            let record = aln.record();
            let flag = record.flags() as u32;
            if flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0 {
                primary_count += 1;
                let mapq = record.mapq();
                mapq_sum += mapq as u32;

                // Count reads with low MAPQ
                if mapq <= self.options.max_low_mapq {
                    low_mapq_count += 1;
                }

                // Count bases that pass both mapping and base quality thresholds
                if mapq >= self.options.min_mapping_quality {
                    if let Some(qual) = aln.qpos().and_then(|q| record.qual().get(q)) {
                        if *qual >= self.options.min_base_quality {
                            qc_depth += 1;
                        }
                    }
                }

                // Track unique reads
                let mut key = record.qname().to_vec();
                key.push(if flag & BAM_FREAD1 != 0 { b'1' } else { b'2' });
                if !self.seen_read_names.contains(&key) {
                    self.seen_read_names.insert(key);
                    self.unique_read_count += 1;
                }
            }
        }

        // Update coverage statistics
        self.total_base_coverage += primary_count as u64;


        if primary_count > 0 {
            self.mapq_sum += mapq_sum as u64;
            self.mapq_count += primary_count as u64;
            
            // Check for poor mapping quality using the new criteria
            if primary_count >= self.options.min_depth_for_low_mapq &&
                (low_mapq_count as f64 / primary_count as f64) > self.options.max_low_mapq_fraction {
                self.stats.poor_mapping_quality += 1;
                return;
            }

            match qc_depth {
                d if d < self.options.min_depth => self.stats.low_coverage += 1,
                d if d > self.options.max_depth => self.stats.excessive_coverage += 1,
                _ => self.stats.callable += 1,
            }
        } else {
            self.stats.no_coverage += 1;
        }

        let is_low_mapq = primary_count >= self.options.min_depth_for_low_mapq &&
            (low_mapq_count as f64 / primary_count as f64) > self.options.max_low_mapq_fraction;

        self.coverage_ranges.push(CoverageRange {
            start: self.current_position,
            end: self.current_position,
            depth: qc_depth,
            is_low_mapq,
        });
        
    }

    fn format_report_line(&self) -> String {
        let total_bases = (self.length - self.stats.ref_n) as f64;
        let avg_depth = self.total_depth as f64 / total_bases;
        let avg_mapq = if self.mapq_count > 0 {
            self.mapq_sum as f64 / self.mapq_count as f64
        } else {
            0.0
        };
        let coverage_percent =
            ((total_bases - self.stats.no_coverage as f64) / total_bases) * 100.0;

        format!(
            "{}|1|{}|{}|{}|{}|{}|{}|{}|{}|{:.2}|{:.2}|{:.1}",
            self.name,
            self.length,
            self.unique_read_count,
            self.stats.ref_n,
            self.stats.no_coverage,
            self.stats.low_coverage,
            self.stats.excessive_coverage,
            self.stats.poor_mapping_quality,
            self.stats.callable,
            coverage_percent,
            avg_depth,
            avg_mapq
        )
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
    summary_writer: BufWriter<File>,
    // Coverage tracking for histograms
    coverage_ranges: Vec<CoverageRange>,
    current_pos: u32,
    output_dir: PathBuf,
}

impl CallableLociCounter {
    pub fn new(bed_file: &str, summary_file: &str) -> IoResult<Self> {
        let output_dir = PathBuf::from(bed_file).parent()
            .unwrap_or(&PathBuf::from("."))
            .to_path_buf();

        Ok(Self {
            counts: [0; 6],
            current_state: None,
            bed_writer: BufWriter::new(File::create(bed_file)?),
            summary_writer: BufWriter::new(File::create(summary_file)?),
            coverage_ranges: Vec::new(),
            current_pos: 0,
            output_dir,
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

    fn finish_contig(&mut self, contig_name: &str) -> IoResult<()> {
        self.write_state()?;

        if !self.coverage_ranges.is_empty() {
            let histogram_path = histogram_plotter::generate_histogram(
                std::mem::take(&mut self.coverage_ranges),
                &format!("{}_coverage", contig_name),
                contig_name,
            )?;

            let final_path = self.output_dir
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
        // Update position tracking
        self.current_pos = pos;

        // Add coverage range
        self.coverage_ranges.push(CoverageRange {
            start: pos,
            end: pos,
            depth,
            is_low_mapq,
        });

        // Update state tracking
        self.counts[state as usize] += 1;
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
}

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
    let mut counter = CallableLociCounter::new(&output_bed, &output_summary)?;

    // Create progress bars
    let multi_progress = MultiProgress::new();
    let main_progress = multi_progress.add(ProgressBar::new_spinner());
    main_progress.set_style(ProgressStyle::default_spinner()
        .template("{spinner:.green} [{elapsed_precise}] {msg}")
        .unwrap());
    main_progress.set_message("Analyzing coverage...");

    // Create contig progress bars map using usize as key
    let mut contig_stats = HashMap::new();

    // First pass to get contig lengths and set up progress bars
    for tid in 0..header.target_names().len() {
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
        } else if raw_depth >= options.min_depth_for_low_mapq &&
            (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction {
            CalledState::POOR_MAPPING_QUALITY
        } else if qc_depth < options.min_depth {
            CalledState::LOW_COVERAGE
        } else if options.max_depth > 0 && qc_depth > options.max_depth {
            CalledState::EXCESSIVE_COVERAGE
        } else {
            CalledState::CALLABLE
        };

        let is_low_mapq = raw_depth >= options.min_depth_for_low_mapq &&
            (low_mapq_count as f64 / raw_depth as f64) > options.max_low_mapq_fraction;

        // Now call process_position with the correct parameters
        counter.process_position(
            contig,
            pos as u32,  // Convert to u32 as required
            qc_depth,
            is_low_mapq,
            state,
        )?;

        // Update contig stats
        if let Some(stats) = contig_stats.get_mut(&tid) {
            stats.process_position(&pileup, raw_depth, pileup.alignments(), ref_base);
        }
    }

    // Finish all progress bars
    for (_, stats) in contig_stats.iter() {
        stats.progress_bar.finish_with_message("Complete");
    }
    main_progress.finish_with_message("Coverage analysis complete!");

    // Finish final contig
    counter.finish_contig(&*current_contig)?;

    let mut summary_writer = BufWriter::new(File::create(output_summary)?);

    // Write header
    writeln!(
        summary_writer,
        "Contig|Version|Length|UniqueReads|RefN|NoCoverage|LowCoverage|ExcessiveCoverage|PoorMappingQuality|Callable|CoveragePercent|AvgDepth|AvgMapQ"
    )?;

    // Write stats for each contig
    for (_, stats) in &contig_stats {
        writeln!(summary_writer, "{}", stats.format_report_line())?;
    }

    Ok(())
}
