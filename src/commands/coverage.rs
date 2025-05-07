use bio::io::fasta::IndexedReader;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::htslib::{BAM_FSECONDARY, BAM_FSUPPLEMENTARY};
use rust_htslib::{bam, bam::Read};
use std::convert::TryInto;
use std::fs::File;
use std::io::{BufWriter, Write};

// TODO FIX-ME: This is SLOW!  Around 6-8 minutes per contig on hs1 aligned samples.

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
    primary_alignments: u64, // Add this field to track primary alignments
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
            primary_alignments: 0,
            stats: CallableStats::default(),
            progress_bar,
            min_depth,
            max_depth,
            min_mapping_quality,
        }
    }

    fn process_position(
        &mut self,
        depth: u32,
        alignments: &[bam::pileup::Alignment],
        ref_base: u8,
    ) {
        self.total_depth += depth as u64;

        // Check if reference base is N
        if ref_base == b'N' || ref_base == b'n' {
            self.stats.ref_n += 1;
            return;
        }

        // Count primary alignments and their mapping qualities
        let primary_alns: Vec<_> = alignments
            .iter()
            .filter(|aln| {
                let flag = aln.record().flags() as u32;
                !(flag & BAM_FSECONDARY != 0 || flag & BAM_FSUPPLEMENTARY != 0)
            })
            .collect();

        // Update primary alignment count
        self.primary_alignments += primary_alns.len() as u64;

        // Calculate average mapping quality using only primary alignments
        let mapping_qualities: Vec<u8> =
            primary_alns.iter().map(|aln| aln.record().mapq()).collect();

        let mapq_sum: u32 = mapping_qualities.iter().map(|&q| q as u32).sum();
        let avg_mapq = if !mapping_qualities.is_empty() {
            (mapq_sum / mapping_qualities.len() as u32) as u8
        } else {
            0
        };

        // Update mapping quality statistics for primary alignments only
        for &mapq in &mapping_qualities {
            if mapq > 0 {
                self.mapq_sum += mapq as u64;
                self.mapq_count += 1;
            }
        }

        // Update callable statistics
        match depth {
            0 => self.stats.no_coverage += 1,
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
    }

    fn format_report_line(&self) -> String {
        let total_bases = self.length as f64;
        let avg_depth = self.total_depth as f64 / total_bases;
        let avg_mapq = if self.mapq_count > 0 {
            self.mapq_sum as f64 / self.mapq_count as f64
        } else {
            0.0
        };
        let coverage_percent = (self.stats.callable as f64 / total_bases) * 100.0;

        format!(
            "{}|1|{}|{}|{}|{}|{}|{}|{}|{}|{:.2}|{:.2}|{:.1}",
            self.name,
            self.length,
            self.primary_alignments,
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
                    &pileup.alignments().collect::<Vec<_>>(),
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
