use crate::cli::Region;
use crate::sequence_processor::core::{ProcessingStats, Sequence, SequenceProcessor};
use crate::sequence_processor::collectors::base::StatsCollector;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use super::utils::{hash_kmer, generate_hexdigest};

#[derive(Debug, Clone)]
pub struct FastFingerprint {
    pub(crate) ksize: usize,
    pub(crate) scaled: usize,
    pub(crate) hashes: HashMap<u64, u32>,
    pub(crate) max_hash: u64,
    pub(crate) max_frequency: Option<u32>,
    pub(crate) progress: Arc<ProgressBar>,
    pub(crate) region: Region,
    current_kmer: Vec<u8>,
}

impl FastFingerprint {
    pub fn new(ksize: usize, scaled: usize, max_frequency: Option<u32>, region: Region) -> Self {
        let max_hash = ((u64::MAX as f64) / scaled as f64) as u64;
        let progress = ProgressBarBuilder::new("")
            .with_template(
                "{spinner:.green} [{elapsed_precise}] {msg} [{wide_bar}] {pos}/{len} ({per_sec})",
            )
            .with_progress_bar()
            .build()
            .unwrap();

        FastFingerprint {
            ksize,
            scaled,
            hashes: HashMap::new(),
            max_hash,
            max_frequency,
            progress: Arc::new(progress),
            region,
            current_kmer: Vec::with_capacity(ksize),
        }
    }

    pub fn get_hexdigest(&self) -> String {
        generate_hexdigest(self)
    }

    pub fn save_hashes(&self, output_path: &Path) -> Result<()> {
        let progress = ProgressBarBuilder::new("Saving k-mer hashes...").build()?;
        let sorted_entries = self.sorted_entries();

        super::utils::save_hashes_to_file(
            self,
            &sorted_entries,
            output_path,
            progress,
        )
    }

    pub(crate) fn sorted_entries(&self) -> Vec<(&u64, &u32)> {
        let mut sorted_entries: Vec<_> = self
            .hashes
            .iter()
            .filter(|(_, &count)| {
                self.max_frequency
                    .map_or(true, |max_freq| count <= max_freq)
            })
            .collect();
        sorted_entries.sort_by_key(|&(hash, _)| hash);
        sorted_entries
    }
}

// Add these implementations to the existing processor.rs file

impl StatsCollector for FastFingerprint {
    fn process_base(&mut self, position: u64, base: u8, _quality: u8, _mapping_quality: u8) -> Result<()> {
        self.current_kmer.push(base);

        if self.current_kmer.len() == self.ksize {
            if !self.current_kmer.iter().any(|&b| b == b'N') {
                let hash = super::utils::hash_kmer(self.ksize, &self.current_kmer);
                if hash <= self.max_hash {
                    *self.hashes.entry(hash).or_insert(0) += 1;
                }
            }
            self.current_kmer.remove(0);
        }
        Ok(())
    }

    fn skip_base(&mut self, _position: u64) -> Result<()> {
        self.current_kmer.clear();
        Ok(())
    }

    fn get_min_base_quality(&self) -> u8 {
        0 // We don't filter by base quality
    }

    fn get_min_mapping_quality(&self) -> u8 {
        0 // We don't filter by mapping quality
    }

    fn merge_with(&mut self, other: &Self) -> Result<()> {
        for (&hash, &count) in &other.hashes {
            *self.hashes.entry(hash).or_insert(0) += count;
        }
        Ok(())
    }
}

impl SequenceProcessor for FastFingerprint {
    fn process_sequence(&mut self, sequence: &Sequence) -> Result<()> {
        if sequence.data.len() < self.ksize {
            return Ok(());
        }

        let num_kmers = sequence.data.len() - self.ksize + 1;
        self.progress.inc(1);
        self.progress
            .set_message(format!("Processing {} bp sequence", sequence.data.len()));

        for i in 0..num_kmers {
            let kmer = &sequence.data[i..i + self.ksize];
            if kmer.iter().any(|&b| b == b'N') {
                continue;
            }

            let hash = super::utils::hash_kmer(self.ksize, kmer);
            if hash <= self.max_hash {
                *self.hashes.entry(hash).or_insert(0) += 1;
            }
        }

        Ok(())
    }

    fn get_min_length(&self) -> usize {
        self.ksize
    }

    fn update_progress(&mut self, stats: &ProcessingStats) {
        self.progress.set_position(stats.processed);
    }

    fn supports_parallel(&self) -> bool {
        true
    }

    fn merge_processor(&mut self, other: &Self) -> Result<()> {
        self.merge_with(other)
    }
}