use anyhow::{Context, Result};
use bio::io::fastq;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read as BamRead};
use seahash::SeaHasher;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::hash::Hasher;
use std::path::{Path, PathBuf};


#[derive(Debug)]
struct FastFingerprint {
    ksize: usize,
    scaled: usize,
    hashes: HashSet<u64>,
    max_hash: u64,
    progress: ProgressBar,
}

impl FastFingerprint {
    fn new(ksize: usize, scaled: usize) -> Self {
        let max_hash = ((u64::MAX as f64) / scaled as f64) as u64;
        let progress = ProgressBar::new(0);
        progress.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] {msg} [{wide_bar}] {pos}/{len} ({per_sec})")
                .unwrap()
                .progress_chars("=>-")
        );
        FastFingerprint {
            ksize,
            scaled,
            hashes: HashSet::new(),
            max_hash,
            progress,
        }
    }

    fn add_sequence(&mut self, sequence: &[u8]) {
        if sequence.len() < self.ksize {
            return;
        }

        let num_kmers = sequence.len() - self.ksize + 1;
        self.progress.inc(1);
        self.progress.set_message(format!("Processing {} bp sequence", sequence.len()));

        for i in 0..num_kmers {
            let kmer = &sequence[i..i + self.ksize];
            if kmer.iter().any(|&b| b == b'N') {
                continue;
            }

            let hash = self.hash_kmer(kmer);
            if hash <= self.max_hash {
                self.hashes.insert(hash);
            }
        }
    }

    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
        let canonical = {
            let rev_comp: Vec<u8> = kmer.iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'N',
                })
                .collect();

            if kmer < rev_comp.as_slice() {
                kmer.to_vec()
            } else {
                rev_comp
            }
        };

        let mut hasher = SeaHasher::new();
        hasher.write(&canonical);
        hasher.finish()
    }

    fn hexdigest(&self) -> String {
        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap()
        );
        progress.set_message("Finalizing fingerprint...");

        let mut hasher = Sha256::new();

        let mut sorted_hashes: Vec<_> = self.hashes.iter().collect();
        sorted_hashes.sort();

        for hash in sorted_hashes {
            hasher.update(hash.to_le_bytes());
        }

        progress.finish_with_message(format!("Generated fingerprint from {} distinct k-mers", self.hashes.len()));
        format!("{:x}", hasher.finalize())
    }
}

use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

pub fn run(
    input_file: String,
    reference_file: Option<String>,
    ksize: usize,
    scaled: usize,
) -> Result<()> {
    let input_path = PathBuf::from(&input_file);
    let mut fp = FastFingerprint::new(ksize, scaled);

    match input_path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") | Some("cram") => process_bam(&input_path, &mut fp, reference_file)?,
        Some(ext) if ext == "fastq" || ext == "fq" || ext == "gz" => process_fastq(&input_path, &mut fp)?,
        _ => anyhow::bail!("Unsupported file format. Must be .fastq, .fq, .fastq.gz, .fq.gz, .bam, or .cram")
    }

    if fp.hashes.is_empty() {
        anyhow::bail!("No valid sequences found in file");
    }

    println!("{}", fp.hexdigest());
    Ok(())
}

fn process_bam(input_path: &Path, fp: &mut FastFingerprint, reference_file: Option<String>) -> Result<()> {
    let mut reader = bam::Reader::from_path(input_path)
        .context("Failed to open alignment file")?;

    if input_path.extension().map_or(false, |ext| ext == "cram") {
        if let Some(ref_path) = &reference_file {
            reader.set_reference(&ref_path)
                .context("Failed to set CRAM reference")?;
        }
    }

    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} [{elapsed_precise}] Processed {human_count} records ({per_sec})")?
    );

    progress.enable_steady_tick(std::time::Duration::from_secs(3));
    progress.set_message("Processing BAM/CRAM file");

    // Create channels for streaming records with larger buffer
    let (tx, rx) = crossbeam_channel::bounded(100_000);
    let (result_tx, result_rx) = crossbeam_channel::bounded(100_000);
    let counter = Arc::new(AtomicU64::new(0));
    let counter_clone = counter.clone();

    // Spawn thread to read records
    std::thread::spawn(move || {
        for record in reader.records() {
            if let Ok(record) = record {
                if tx.send(record).is_err() {
                    break;
                }
            }
        }
    });

    // Spawn worker threads for processing
    let num_workers = rayon::current_num_threads();
    let rx = Arc::new(rx);
    let mut workers = Vec::with_capacity(num_workers);

    for _ in 0..num_workers {
        let rx = rx.clone();
        let result_tx = result_tx.clone();
        let ksize = fp.ksize;

        workers.push(std::thread::spawn(move || {
            while let Ok(record) = rx.recv() {
                if !record.is_unmapped() {
                    let seq = record.seq().as_bytes();
                    if seq.len() >= ksize {
                        // Send processed sequences to result channel
                        if result_tx.send(seq.to_vec()).is_err() {
                            break;
                        }
                    }
                }
            }
        }));
    }

    // Spawn thread to close result sender when all workers are done
    std::thread::spawn(move || {
        for worker in workers {
            let _ = worker.join();
        }
        // result_tx is dropped here, closing the channel
    });

    // Process results as they arrive
    while let Ok(seq) = result_rx.recv() {
        fp.add_sequence(&seq);
        counter_clone.fetch_add(1, Ordering::Relaxed);
        progress.set_position(counter_clone.load(Ordering::Relaxed));
    }

    progress.finish_and_clear();
    Ok(())
}

fn process_fastq(input_path: &Path, fp: &mut FastFingerprint) -> Result<()> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")?
    );

    progress.set_message("Processing FASTQ file...");

    let reader = fastq::Reader::new(std::fs::File::open(input_path)
        .context("Failed to open FASTQ file")?);

    // Create channels for streaming records
    let (tx, rx) = crossbeam_channel::bounded(10_000);
    let counter = Arc::new(AtomicU64::new(0));
    let counter_clone = counter.clone();

    // Spawn thread to read records
    std::thread::spawn(move || {
        for record in reader.records() {
            if let Ok(record) = record {
                if tx.send(record).is_err() {
                    break;
                }
            }
        }
    });

    // Process records as they arrive
    for record in rx {
        let seq = record.seq();
        if seq.len() >= fp.ksize {
            fp.add_sequence(seq);
        }
        counter_clone.fetch_add(1, Ordering::Relaxed);
        progress.set_position(counter_clone.load(Ordering::Relaxed));
    }

    progress.finish_and_clear();
    Ok(())
}
