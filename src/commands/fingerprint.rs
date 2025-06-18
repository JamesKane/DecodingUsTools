use crate::cli::Region;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use crate::utils::sequence_processor::core::SequenceReader;
use crate::utils::sequence_processor::core::{ProcessingStats, Sequence, SequenceProcessor};
use crate::utils::sequence_processor::readers::{BamReader, FastqReader};
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use seahash::SeaHasher;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::hash::Hasher;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
struct FastFingerprint {
    ksize: usize,
    scaled: usize,
    hashes: HashSet<u64>,
    max_hash: u64,
    progress: Arc<ProgressBar>,
    region: Region,
}

impl FastFingerprint {
    fn new(ksize: usize, scaled: usize, region: Region) -> Self {
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
            hashes: HashSet::new(),
            max_hash,
            progress: Arc::new(progress),
            region,
        }
    }

    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
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
        let progress = ProgressBarBuilder::new("Finalizing fingerprint...")
            .build()
            .unwrap();

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
        let progress = ProgressBarBuilder::new("Saving k-mer hashes...").build()?;

        let mut sorted_hashes: Vec<_> = self.hashes.iter().collect();
        sorted_hashes.sort();

        let file = std::fs::File::create(output_path).context("Failed to create output file")?;
        let mut writer = std::io::BufWriter::new(file);

        writeln!(writer, "#ksize={}", self.ksize)?;
        writeln!(writer, "#scaled={}", self.scaled)?;
        writeln!(writer, "#region={}", self.region.to_output_name())?;

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

            let hash = self.hash_kmer(kmer);
            if hash <= self.max_hash {
                self.hashes.insert(hash);
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
        self.hashes.extend(other.hashes.iter().copied());
        Ok(())
    }
}

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
    let progress = ProgressBarBuilder::new("Processing").build()?;

    let stats = match input_path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") | Some("cram") => {
            let mut reader = BamReader::new(&input_path, reference_file.as_deref())?;
            reader.read_sequences(&mut fp, &progress)?
        }
        Some(ext) if ext == "fastq" || ext == "fq" || ext == "gz" => {
            let mut reader = FastqReader::new(&input_path)?;
            reader.read_sequences(&mut fp, &progress)?
        }
        _ => anyhow::bail!(
            "Unsupported file format. Must be .fastq, .fq, .fastq.gz, .fq.gz, .bam, or .cram"
        ),
    };

    if fp.hashes.is_empty() {
        anyhow::bail!("No valid sequences found in file");
    }

    println!("Processed {} sequences", stats.processed);
    println!("{}", fp.hexdigest());

    if let Some(output_path) = output_file {
        fp.save_hashes(Path::new(&output_path))?;
    }

    Ok(())
}
