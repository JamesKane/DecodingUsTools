use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::htslib::{BAM_FSECONDARY, BAM_FSUPPLEMENTARY};
use rust_htslib::{bam, bam::Read};
use std::collections::HashSet;
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufWriter, Write};

// TODO FIX-ME: This is still kinda slow.  About 2 minutes per contig in release mode

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
    min_depth: u32,
    max_depth: u32,
    min_mapping_quality: u8,
}

impl ContigStats {
    fn new(
        name: String,
        length: usize,
        multi_progress: &MultiProgress,
        min_depth: u32,
        max_depth: u32,
        min_mapping_quality: u8,
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
            min_depth,
            max_depth,
            min_mapping_quality,
        }
    }

    fn process_position(&mut self, depth: u32, alignments: bam::pileup::Alignments, ref_base: u8) {
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
        let mut mapq_sum: u32 = 0;

        for aln in alignments {
            let record = aln.record();
            let flag = record.flags() as u32;
            if flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY) == 0 {
                primary_count += 1;
                mapq_sum += record.mapq() as u32;

                // Track unique reads
                if !self.seen_read_names.contains(record.qname()) {
                    self.seen_read_names.insert(record.qname().to_vec());
                    self.unique_read_count += 1;
                }
            }
        }

        // Update coverage statistics
        self.total_base_coverage += primary_count as u64;

        if primary_count > 0 {
            let avg_mapq = (mapq_sum / primary_count) as u8;
            self.mapq_sum += mapq_sum as u64;
            self.mapq_count += primary_count as u64;

            match depth {
                d if d <= self.min_depth => self.stats.low_coverage += 1,
                d if d > self.max_depth => self.stats.excessive_coverage += 1,
                _ => {
                    if avg_mapq < self.min_mapping_quality {
                        self.stats.poor_mapping_quality += 1;
                    } else {
                        self.stats.callable += 1;
                    }
                }
            }
        } else {
            self.stats.no_coverage += 1;
        }
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

pub fn run(
    bam_file: String,
    reference_file: String,
    output_file: String,
    min_depth: u32,
    max_depth: u32,
    min_mapping_quality: u8,
) -> Result<(), Box<dyn std::error::Error>> {
    // Open the output file
    let output_file = File::create(&output_file)?;
    let mut writer = BufWriter::new(output_file);

    // Write the header
    writeln!(writer, "{}",
             "CONTIG|START_POS|END_POS|NUM_READS|REF_N|NO_COV|LOW_COV|EXCESSIVE_COV|POOR_MQ|CALLABLE|COV_PERCENT|MEAN_DEPTH|MEAN_MQ"
    )?;

    let mut bam = bam::Reader::from_path(&bam_file)?;
    bam.set_threads(4)?;
    let header = bam::Header::from_template(bam.header());
    let bam_header = bam.header().clone();

    // Get contig information
    let mut contig_lengths = Vec::new();
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for record in records {
                if let (Some(sn), Some(ln)) = (record.get("SN"), record.get("LN")) {
                    if let Ok(length) = ln.parse::<usize>() {
                        contig_lengths.push((sn.to_string(), length));
                    }
                }
            }
        }
    }

    let mut pileup = bam.pileup();
    pileup.set_max_depth(1000000);

    let multi_progress = MultiProgress::new();

    let mut fasta_reader = IndexedReader::from_file(&reference_file)?;
    let mut current_contig: Option<ContigStats> = None;
    let mut current_seq = Vec::new();

    // Process pileups
    for p in pileup {
        let pileup = p.unwrap();
        let tid: i32 = pileup.tid().try_into().unwrap();
        let tid_u32: u32 = tid.try_into().unwrap();
        let ref_name = String::from_utf8(bam_header.tid2name(tid_u32).to_owned()).unwrap();
        let pos = pileup.pos() as usize;

        // If we've moved to a new contig
        if current_contig.as_ref().map_or(true, |c| c.name != ref_name) {
            if let Some(stats) = current_contig.take() {
                stats
                    .progress_bar
                    .finish_with_message(format!("Completed {}", stats.name));
                writeln!(writer, "{}", stats.format_report_line()).unwrap();
            }

            if let Some(&(_, length)) = contig_lengths.iter().find(|(name, _)| name == &ref_name) {
                // Read the sequence for the new contig
                current_seq.clear();
                fasta_reader.fetch(&ref_name, 0, length as u64)?;
                fasta_reader.read(&mut current_seq)?;

                current_contig = Some(ContigStats::new(
                    ref_name,
                    length,
                    &multi_progress,
                    min_depth,
                    max_depth,
                    min_mapping_quality,
                ));
            }
        }

        if let Some(ref mut contig_stats) = current_contig {
            if pos < contig_stats.length {
                let depth = pileup.depth();
                let ref_base = if pos < current_seq.len() {
                    current_seq[pos]
                } else {
                    b'N'
                };
                contig_stats.process_position(
                    depth,
                    pileup.alignments(), // Pass iterator directly instead of collecting
                    ref_base,
                );
                contig_stats.progress_bar.set_position(pos as u64);
            }
        }
    }

    // Print stats for the last contig
    if let Some(stats) = current_contig {
        stats
            .progress_bar
            .finish_with_message(format!("Completed {}", stats.name));
        writeln!(writer, "{}", stats.format_report_line()).unwrap();
    }

    // Wait for all progress bars to finish
    multi_progress.clear().unwrap();

    Ok(())
}
