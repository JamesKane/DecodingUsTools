use super::processor::FastFingerprint;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::{Context, Result};
use indicatif::ProgressBar;
use seahash::SeaHasher;
use sha2::{Digest, Sha256};
use std::hash::Hasher;
use std::io::Write;
use std::path::Path;

pub(crate) fn hash_kmer(ksize: usize, kmer: &[u8]) -> u64 {
    let mut canonical = Vec::with_capacity(ksize);
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

pub(crate) fn generate_hexdigest(fp: &FastFingerprint) -> String {
    let progress = ProgressBarBuilder::new("Finalizing fingerprint...")
        .build()
        .unwrap();

    let mut hasher = Sha256::new();
    let sorted_entries = fp.sorted_entries();

    for (hash, count) in sorted_entries {
        hasher.update(hash.to_le_bytes());
        hasher.update(count.to_le_bytes());
    }

    progress.finish_with_message(format!(
        "Generated fingerprint from {} distinct k-mers",
        fp.hashes.len()
    ));
    format!("{:x}", hasher.finalize())
}

pub(crate) fn save_hashes_to_file(
    fp: &FastFingerprint,
    sorted_entries: &[(&u64, &u32)],
    output_path: &Path,
    progress: ProgressBar,
) -> Result<()> {
    let file = std::fs::File::create(output_path).context("Failed to create output file")?;
    let mut writer = std::io::BufWriter::new(file);

    writeln!(writer, "#ksize={}", fp.ksize)?;
    writeln!(writer, "#scaled={}", fp.scaled)?;
    writeln!(writer, "#region={}", fp.region.to_output_name())?;
    if let Some(max_freq) = fp.max_frequency {
        writeln!(writer, "#max_frequency={}", max_freq)?;
    }

    for (&hash, &count) in sorted_entries {
        writeln!(writer, "{}\t{}", hash, count).context("Failed to write hash and count")?;
    }

    progress.finish_with_message(format!(
        "Saved {} k-mer hashes with abundances to {}",
        sorted_entries.len(),
        output_path.display()
    ));
    Ok(())
}