use anyhow::{Context, Result};
use bio::io::fastq;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::{bam, bam::Read as BamRead};
use seahash::SeaHasher;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::hash::Hasher;
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct FastFingerprint {
    ksize: usize,
    scaled: usize,
    hashes: HashSet<u64>,
    max_hash: u64,
    progress: ProgressBar,
    region: Region,
}

impl FastFingerprint {
    fn new(ksize: usize, scaled: usize, region: Region) -> Self {
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
            region,
        }
    }

    fn add_sequence(&mut self, sequence: &[u8]) {
        if sequence.len() < self.ksize {
            return;
        }

        let num_kmers = sequence.len() - self.ksize + 1;
        self.progress.inc(1);
        self.progress
            .set_message(format!("Processing {} bp sequence", sequence.len()));

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
        // Ensure aligned memory access by copying to a new Vec
        let mut canonical = Vec::with_capacity(self.ksize);
        let rev_comp: Vec<u8> = kmer
            .iter()
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
            canonical.extend_from_slice(kmer);
        } else {
            canonical = rev_comp;
        }

        let mut hasher = SeaHasher::new();
        hasher.write(&canonical);
        hasher.finish()
    }

    fn hexdigest(&self) -> String {
        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        progress.set_message("Finalizing fingerprint...");

        let mut hasher = Sha256::new();

        let mut sorted_hashes: Vec<_> = self.hashes.iter().collect();
        sorted_hashes.sort();

        for hash in sorted_hashes {
            hasher.update(hash.to_le_bytes());
        }

        progress.finish_with_message(format!(
            "Generated fingerprint from {} distinct k-mers",
            self.hashes.len()
        ));
        format!("{:x}", hasher.finalize())
    }

    fn save_hashes(&self, output_path: &Path) -> Result<()> {
        let progress = ProgressBar::new_spinner();
        progress.set_style(ProgressStyle::default_spinner().template("{spinner:.green} {msg}")?);
        progress.set_message("Saving k-mer hashes...");

        let mut sorted_hashes: Vec<_> = self.hashes.iter().collect();
        sorted_hashes.sort();

        let file = std::fs::File::create(output_path).context("Failed to create output file")?;
        let mut writer = std::io::BufWriter::new(file);

        // Write header with parameters
        writeln!(writer, "#ksize={}", self.ksize).context("Failed to write ksize")?;
        writeln!(writer, "#scaled={}", self.scaled).context("Failed to write scaled factor")?;
        writeln!(writer, "#region={:?}", self.region).context("Failed to write region")?;

        // Write hashes one per line
        for &hash in &sorted_hashes {
            writeln!(writer, "{}", hash).context("Failed to write hash")?;
        }

        progress.finish_with_message(format!(
            "Saved {} k-mer hashes to {}",
            self.hashes.len(),
            output_path.display()
        ));
        Ok(())
    }
}

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use crate::cli::Region;

pub fn run(
    input_file: String,
    reference_file: Option<String>,
    ksize: usize,
    scaled: usize,
    output_file: Option<String>,
    region: Region,
) -> Result<()> {
    let input_path = PathBuf::from(&input_file);
    let mut fp = FastFingerprint::new(ksize, scaled, region);

    match input_path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") | Some("cram") => process_bam(&input_path, &mut fp, reference_file)?,
        Some(ext) if ext == "fastq" || ext == "fq" || ext == "gz" => {
            process_fastq(&input_path, &mut fp)?
        }
        _ => anyhow::bail!(
            "Unsupported file format. Must be .fastq, .fq, .fastq.gz, .fq.gz, .bam, or .cram"
        ),
    }

    if fp.hashes.is_empty() {
        anyhow::bail!("No valid sequences found in file");
    }

    println!("{}", fp.hexdigest());

    // Save hashes if output file is specified
    if let Some(output_path) = output_file {
        fp.save_hashes(Path::new(&output_path))?;
    }

    Ok(())
}

fn process_bam(
    input_path: &Path,
    fp: &mut FastFingerprint,
    reference_file: Option<String>,
) -> Result<()> {
    let mut reader = bam::IndexedReader::from_path(input_path).context("Failed to open alignment file")?;

    if input_path.extension().map_or(false, |ext| ext == "cram") {
        if let Some(ref_path) = &reference_file {
            reader
                .set_reference(&ref_path)
                .context("Failed to set CRAM reference")?;
        }
    }

    let target_chromosomes = fp.region.to_chromosome_names();

    // If no specific region is requested, process everything
    if target_chromosomes.is_empty() {
        return process_all_regions(&mut reader, fp);
    }

    // Process each requested chromosome
    for chr_name in target_chromosomes {
        // Get chromosome ID and length from a fresh header reference
        let tid = match reader.header().tid(chr_name.as_bytes()) {
            Some(tid) => tid,
            None => continue, // Skip if chromosome not found
        };

        let chr_len = match reader.header().target_len(tid) {
            Some(len) => len,
            None => continue, // Skip if no length information
        };

        // Set up region for fetching
        match reader.fetch((tid, 0, chr_len)) {
            Ok(_) => {
                println!("Processing chromosome: {} (length: {})", chr_name, chr_len);
                process_region(&mut reader, fp)?;
            },
            Err(e) => {
                eprintln!(
                    "Warning: Failed to fetch chromosome {}. Error: {}. Skipping.",
                    chr_name, e
                );
                continue;
            }
        }
    }

    if fp.hashes.is_empty() {
        return Err(anyhow::anyhow!(
            "No valid sequences found in specified regions"
        ));
    }

    Ok(())
}

fn process_region<R: bam::Read>(
    reader: &mut R,
    fp: &mut FastFingerprint,
) -> Result<()> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(ProgressStyle::default_spinner().template(
        "{spinner:.green} [{elapsed_precise}] Processed {human_count} records ({per_sec})",
    )?);

    progress.enable_steady_tick(std::time::Duration::from_secs(3));

    for r in reader.records() {
        let record = r?;
        if !record.is_unmapped() {
            let seq = record.seq().as_bytes();
            if seq.len() >= fp.ksize {
                fp.add_sequence(&seq);
                progress.inc(1);
            }
        }
    }

    progress.finish_and_clear();
    Ok(())
}

fn process_all_regions<R: bam::Read>(
    reader: &mut R,
    fp: &mut FastFingerprint,
) -> Result<()> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(ProgressStyle::default_spinner().template(
        "{spinner:.green} [{elapsed_precise}] Processed {human_count} records ({per_sec})",
    )?);

    progress.enable_steady_tick(std::time::Duration::from_secs(3));
    progress.set_message("Processing entire BAM file");

    for r in reader.records() {
        let record = r?;
        if !record.is_unmapped() {
            let seq = record.seq().as_bytes();
            if seq.len() >= fp.ksize {
                fp.add_sequence(&seq);
                progress.inc(1);
            }
        }
    }

    progress.finish_and_clear();
    Ok(())
}

fn process_fastq(input_path: &Path, fp: &mut FastFingerprint) -> Result<()> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(ProgressStyle::default_spinner().template("{spinner:.green} {msg}")?);

    progress.set_message("Processing FASTQ file...");

    let reader =
        fastq::Reader::new(std::fs::File::open(input_path).context("Failed to open FASTQ file")?);

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
