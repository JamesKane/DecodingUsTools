use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use indicatif::{ProgressBar, ProgressStyle};

pub struct BamStats {
    read_count: usize,
    total_read_length: usize,
    paired_reads: usize,
    total_insert_size: i64,
    paired_count: usize,
    pub(crate) max_samples: usize,
    length_distribution: HashMap<usize, usize>,
    insert_size_distribution: HashMap<i64, usize>,
}

impl BamStats {
    pub fn new(max_samples: usize) -> Self {
        BamStats {
            read_count: 0,
            total_read_length: 0,
            paired_reads: 0,
            total_insert_size: 0,
            paired_count: 0,
            max_samples,
            length_distribution: HashMap::new(),
            insert_size_distribution: HashMap::new(),
        }
    }

    pub fn collect_stats(&mut self, bam_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut bam = bam::Reader::from_path(bam_path)?;

        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] Collecting BAM statistics...")
                .unwrap(),
        );

        for (i, record_result) in bam.records().enumerate() {
            if i >= self.max_samples {
                break;
            }

            let record = record_result?;

            // Only count primary alignments
            if !record.is_secondary() && !record.is_supplementary() {
                let seq = record.seq();
                let seq_len = seq.len();

                *self.length_distribution.entry(seq_len).or_insert(0) += 1;

                self.read_count += 1;
                self.total_read_length += seq_len;

                if record.is_paired() {
                    self.paired_reads += 1;

                    if record.is_proper_pair() && record.is_first_in_template() {
                        let insert_size = record.insert_size().abs();
                        if insert_size > 0 {
                            *self.insert_size_distribution.entry(insert_size).or_insert(0) += 1;
                            self.total_insert_size += insert_size;
                            self.paired_count += 1;
                        }
                    }
                }
            }

            if i % 10000 == 0 {
                progress.set_message(format!("Processed {} reads...", i));
            }
        }

        progress.finish_with_message("BAM Statistics collected");
        Ok(())
    }

    pub fn get_stats(&self) -> HashMap<String, f64> {
        let mut stats = HashMap::new();

        if self.read_count > 0 {
            let modal_length = self.length_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(length, _)| *length)
                .unwrap_or(0);

            stats.insert(
                "average_read_length".to_string(),
                modal_length as f64,
            );
            stats.insert(
                "paired_percentage".to_string(),
                (self.paired_reads as f64 / self.read_count as f64) * 100.0,
            );
        }

        if self.paired_count > 0 {
            let modal_insert_size = self.insert_size_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(size, _)| *size)
                .unwrap_or(0);

            stats.insert(
                "average_insert_size".to_string(),
                modal_insert_size as f64,
            );
            stats.insert(
                "proper_pair_percentage".to_string(),
                (self.paired_count as f64 * 2.0 / self.paired_reads as f64) * 100.0,
            );
        }

        stats
    }
}