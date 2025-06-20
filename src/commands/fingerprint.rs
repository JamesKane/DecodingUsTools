use crate::cli::Region;
use crate::sequence_processor::core::SequenceReader;
use crate::sequence_processor::collectors::fingerprint::processor::FastFingerprint;
use crate::sequence_processor::readers::{BamReader, FastqReader, GamReader};
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::Result;
use std::path::{Path, PathBuf};

pub fn run(
    input_file: String,
    reference_file: Option<String>,
    ksize: usize,
    scaled: usize,
    max_frequency: Option<u32>,
    output_file: Option<String>,
    region: Region,
) -> Result<()> {
    let input_path = PathBuf::from(&input_file);
    let mut fp = FastFingerprint::new(ksize, scaled, max_frequency, region);
    let progress = ProgressBarBuilder::new("Processing").build()?;

    let num_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);

    let stats = match input_path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") | Some("cram") => {
            let mut reader = BamReader::new(&input_path, reference_file.as_deref())?;
            reader.read_sequences_with_threads(&mut fp, &progress, num_threads)?
        }
        Some("gam") => {
            let mut reader = GamReader::new(&input_path)?;
            reader.read_sequences_with_threads(&mut fp, &progress, num_threads)?
        }
        Some(ext) if ext == "fastq" || ext == "fq" || ext == "gz" => {
            let mut reader = FastqReader::new(&input_path)?;
            reader.read_sequences_with_threads(&mut fp, &progress, num_threads)?
        }
        _ => anyhow::bail!(
            "Unsupported file format. Must be .fastq, .fq, .fastq.gz, .fq.gz, .bam, or .cram"
        ),
    };

    println!("Processed {} sequences", stats.processed);
    println!("{}", fp.get_hexdigest());

    if let Some(output_path) = output_file {
        fp.save_hashes(Path::new(&output_path))?;
    }

    Ok(())
}
