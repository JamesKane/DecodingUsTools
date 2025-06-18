use crate::cli::Region;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::{Context, Result};
use bio::io::fastq;
use indicatif::{ProgressBar, ProgressStyle};
use niffler::get_reader;
use rust_htslib::{bam, bam::Read as BamRead};
use seahash::SeaHasher;
use sha2::{Digest, Sha256};
use std::collections::HashSet;
use std::fs::File;
use std::hash::Hasher;
use std::io::{BufReader, Write};
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
        let progress = ProgressBarBuilder::new("Saving k-mer hashes...")
            .build()?;

        let mut sorted_hashes: Vec<_> = self.hashes.iter().collect();
        sorted_hashes.sort();

        let file = std::fs::File::create(output_path).context("Failed to create output file")?;
        let mut writer = std::io::BufWriter::new(file);

        // Write header with parameters
        writeln!(writer, "#ksize={}", self.ksize)?;
        writeln!(writer, "#scaled={}", self.scaled)?;
        writeln!(writer, "#region={}", self.region.to_output_name())?;

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
            let count = process_fastq(&input_path, &mut fp)?;
            println!("Processed {} valid sequences", count);
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
    let mut reader =
        bam::IndexedReader::from_path(input_path).context("Failed to open alignment file")?;

    if input_path.extension().map_or(false, |ext| ext == "cram") {
        if let Some(ref_path) = &reference_file {
            reader
                .set_reference(&ref_path)
                .context("Failed to set CRAM reference")?;
        }
    }

    let header = reader.header().clone();
    let target_chromosomes = fp.region.to_chromosome_names();

    let progress = ProgressBarBuilder::new("Processing")
        .with_template("{spinner:.green} [{elapsed_precise}] Processing {msg}")
        .with_tick()
        .build()?;
    progress.enable_steady_tick(std::time::Duration::from_secs(5));

    // If no specific region is requested (full genome), process all records
    println!("Processing Sequences file... {:?}", target_chromosomes);
    if target_chromosomes.is_empty() {
        progress.set_message("entire genome");

        process_all_records(input_path, fp, &progress)?;

        return Ok(());
    }

    // Process specific chromosomes
    let mut found_any = false;
    for chr_name in target_chromosomes {
        if let Some(tid) = header.tid(chr_name.as_bytes()) {
            if let Some(chr_len) = header.target_len(tid) {
                progress.set_message(format!("chromosome {}", chr_name));
                match reader.fetch((tid, 0, chr_len)) {
                    Ok(_) => {
                        process_records(&mut reader, fp, &progress)?;
                        found_any = true;
                    }
                    Err(e) => {
                        eprintln!(
                            "Warning: Failed to fetch {}. Error: {}. Skipping.",
                            chr_name, e
                        );
                    }
                }
            }
        }
    }

    if !found_any {
        return Err(anyhow::anyhow!(
            "No valid sequences found in specified regions"
        ));
    }

    progress.finish_and_clear();
    Ok(())
}

fn process_records<R: bam::Read>(
    reader: &mut R,
    fp: &mut FastFingerprint,
    progress: &ProgressBar,
) -> Result<bool> {
    // Changed return type to Result<bool>
    let mut processed = 0;
    let chunk_size = 10000; // Process records in chunks

    let mut records = Vec::with_capacity(chunk_size);

    loop {
        records.clear();

        // Collect a chunk of records
        for _ in 0..chunk_size {
            match reader.records().next() {
                Some(Ok(record)) => {
                    let seq = record.seq().as_bytes();
                    if seq.len() >= fp.ksize {
                        records.push(seq.to_vec());
                    }
                }
                Some(Err(e)) => eprintln!("Warning: Error reading record: {}", e),
                None => break,
            }
        }

        if records.is_empty() && processed == 0 {
            break;
        }

        // Process the chunk
        for seq in &records {
            fp.add_sequence(seq);
            processed += 1;
            if processed % 1000 == 0 {
                progress.set_position(processed);
            }
        }

        if records.len() < chunk_size {
            break;
        }
    }

    Ok(processed > 0) // Return true if we processed any records
}

fn process_all_records(
    input_path: &Path,
    fp: &mut FastFingerprint,
    progress: &ProgressBar,
) -> Result<()> {
    // Use regular Reader for full genome processing
    let mut reader = bam::Reader::from_path(input_path).context("Failed to open BAM file")?;

    let mut processed = 0;
    let mut unmapped = 0;
    let mut too_short = 0;

    println!("Starting full genome processing...");

    let mut record = bam::Record::new();

    while let Some(result) = reader.read(&mut record) {
        match result {
            Ok(()) => {
                if record.is_unmapped() {
                    unmapped += 1;
                }

                let seq = record.seq().as_bytes();
                if seq.len() >= fp.ksize {
                    fp.add_sequence(&seq);
                    processed += 1;

                    if processed % 10000 == 0 {
                        progress.set_message(format!("Processed {} sequences", processed));
                        progress.set_position(processed as u64);
                    }
                } else {
                    too_short += 1;
                }
            }
            Err(e) => {
                eprintln!("Warning: Error reading record: {}", e);
            }
        }
    }

    println!("\nFinal statistics:");
    println!("  Total sequences processed: {}", processed);
    println!("  Unmapped records: {}", unmapped);
    println!("  Too short sequences: {}", too_short);

    if processed == 0 {
        return Err(anyhow::anyhow!("No valid sequences were processed"));
    }

    Ok(())
}

fn process_fastq(input_path: &Path, fp: &mut FastFingerprint) -> Result<u64> {
    let progress = ProgressBarBuilder::new("Processing FASTQ file...").build()?;

    println!("Opening FASTQ file: {}", input_path.display());

    let mut processed_sequences = 0u64;
    println!("Starting FASTQ processing");

    // Create boxed reader for niffler
    let file = File::open(input_path)?;
    let (reader, _compression) = get_reader(Box::new(file))?;
    let reader = BufReader::new(reader);
    let fastq_reader = fastq::Reader::new(reader);

    for (i, record) in fastq_reader.records().enumerate() {
        match record {
            Ok(record) => {
                fp.add_sequence(&record.seq());
                processed_sequences += 1;

                if processed_sequences % 1_000_000 == 0 {
                    progress.set_message(format!(
                        "Processed {} million sequences",
                        processed_sequences / 1_000_000
                    ));
                }
            }
            Err(e) => {
                eprintln!("Warning: Error reading record at position {}: {}", i, e);
            }
        }
    }

    progress.finish_with_message(format!(
        "Finished processing {} sequences",
        processed_sequences
    ));
    Ok(processed_sequences)
}
