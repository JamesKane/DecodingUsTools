use anyhow::{anyhow, Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::{self, Header, Read, Reader, Writer};
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::process::Command;

#[derive(Debug, Clone, Copy)]
pub enum ReferenceGenome {
    GRCh38,
    CHM13,
}

impl ReferenceGenome {
    fn prefix(&self) -> &str {
        match self {
            ReferenceGenome::GRCh38 => "GRCh38#0#",
            ReferenceGenome::CHM13 => "CHM13#0#",
        }
    }
}

pub struct BamFixer {
    reference_file: String,
    input_bam: String,
    output_bam: String,
    keep_temp: bool,
    temp_original_non_sq_header: String,
    temp_reference_sq_header: String,
    temp_combined_header: String,
    reference_genome: ReferenceGenome,
}

impl BamFixer {
    pub fn new(
        reference_file: String,
        input_bam: String,
        output_bam: String,
        keep_temp: bool,
    ) -> Self {
        let reference_genome = Self::detect_reference_genome(&input_bam)
            .expect("Failed to detect reference genome type");

        let temp_dir = std::env::temp_dir();
        let base_name = Path::new(&input_bam)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();

        BamFixer {
            reference_file,
            input_bam,
            output_bam,
            keep_temp,
            temp_original_non_sq_header: temp_dir
                .join(format!("{}.non_sq.sam", base_name))
                .to_string_lossy()
                .to_string(),
            temp_reference_sq_header: temp_dir
                .join(format!("{}.ref_sq.sam", base_name))
                .to_string_lossy()
                .to_string(),
            temp_combined_header: temp_dir
                .join(format!("{}.combined.sam", base_name))
                .to_string_lossy()
                .to_string(),
            reference_genome,
        }
    }

    fn detect_reference_genome(bam_path: &str) -> Result<ReferenceGenome> {
        let reader = Reader::from_path(bam_path)?;
        let header = reader.header();

        // Check first reference sequence name
        if let Some(name) = header.target_names().first() {
            let name = String::from_utf8_lossy(name);
            if name.starts_with("GRCh38#0#") {
                Ok(ReferenceGenome::GRCh38)
            } else if name.starts_with("CHM13#0#") {
                Ok(ReferenceGenome::CHM13)
            } else {
                Err(anyhow!("Unknown reference genome format"))
            }
        } else {
            Err(anyhow!("No reference sequences found in BAM header"))
        }
    }

    pub fn run(&self) -> Result<()> {
        // Generate reference header
        self.generate_reference_header()?;

        // Extract non-SQ headers from original BAM
        self.extract_non_sq_headers()?;

        // Combine headers
        self.combine_headers()?;

        // Process and fix the BAM file
        self.fix_bam()?;

        // Clean up unless keep_temp is true
        if !self.keep_temp {
            self.cleanup()?;
        }

        Ok(())
    }

    fn generate_reference_header(&self) -> Result<()> {
        // Check if the FASTA index exists
        let fai_path = format!("{}.fai", self.reference_file);
        if !Path::new(&fai_path).exists() {
            anyhow::bail!(
                "Reference FASTA index not found at {}. \
            Please run 'samtools faidx {}' first",
                fai_path,
                self.reference_file
            );
        }

        // Use samtools dict to generate header
        let output = Command::new("samtools")
            .args(&[
                "dict",
                "-o",
                &self.temp_reference_sq_header,
                &self.reference_file,
            ])
            .output()
            .context("Failed to run samtools dict")?;

        if !output.status.success() {
            anyhow::bail!(
                "samtools dict failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }

        Ok(())
    }

    fn extract_non_sq_headers(&self) -> Result<()> {
        let mut reader = Reader::from_path(&self.input_bam).context("Failed to open input BAM")?;

        let header = reader.header().clone();
        let mut output = File::create(&self.temp_original_non_sq_header)
            .context("Failed to create temporary non-SQ header file")?;

        // Convert header to text and filter @SQ lines
        let header_text = String::from_utf8_lossy(header.as_bytes());
        for line in header_text.lines() {
            if !line.starts_with("@SQ") {
                writeln!(output, "{}", line)?;
            }
        }

        Ok(())
    }

    fn combine_headers(&self) -> Result<()> {
        let mut combined = File::create(&self.temp_combined_header)
            .context("Failed to create combined header file")?;

        // Copy non-SQ headers
        io::copy(
            &mut File::open(&self.temp_original_non_sq_header)?,
            &mut combined,
        )?;

        // Copy SQ headers
        io::copy(
            &mut File::open(&self.temp_reference_sq_header)?,
            &mut combined,
        )?;

        Ok(())
    }

    fn fix_bam(&self) -> Result<()> {
        // Read the new header
        let header_reader = Reader::from_path(&self.temp_combined_header)?;
        let header = Header::from_template(header_reader.header());
        let header_view = header_reader.header();

        // Create valid reference name set
        let valid_refs: HashSet<_> = (0..header_view.target_names().len())
            .filter_map(|tid| {
                let name = header_view.tid2name(tid as u32);
                if name.is_empty() {
                    None
                } else {
                    Some(String::from_utf8_lossy(name).to_string())
                }
            })
            .collect();

        // Create mapping of reference names to their new TIDs
        let ref_name_to_tid: HashMap<String, i32> = (0..header_view.target_names().len())
            .filter_map(|tid| {
                let name = header_view.tid2name(tid as u32);
                if name.is_empty() {
                    None
                } else {
                    Some((String::from_utf8_lossy(name).to_string(), tid as i32))
                }
            })
            .collect();

        // Setup BAM reader and writer
        let mut reader = Reader::from_path(&self.input_bam)?;
        let mut writer = Writer::from_path(&self.output_bam, &header, bam::Format::Bam)?;

        // Clone the header before starting the loop
        let input_header = reader.header().clone();

        // Setup progress bar
        let total_reads = reader.records().count() as u64;
        let progress = ProgressBar::new(total_reads);
        progress.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
                )?
                .progress_chars("#>-"),
        );

        // Reset reader for the main loop
        let mut reader = Reader::from_path(&self.input_bam)?;

        let mut renamed_count = 0;
        let mut unmapped_count = 0;

        for (i, record_result) in reader.records().enumerate() {
            let mut record = record_result?;

            if !record.is_unmapped() {
                let ref_name = input_header.tid2name(record.tid() as u32);
                if !ref_name.is_empty() {
                    let ref_name = String::from_utf8_lossy(ref_name).to_string();
                    if ref_name.starts_with(self.reference_genome.prefix()) {
                        let stripped_name = ref_name
                            .trim_start_matches(self.reference_genome.prefix())
                            .to_string();

                        if valid_refs.contains(&stripped_name) {
                            if let Some(&new_tid) = ref_name_to_tid.get(&stripped_name) {
                                record.set_tid(new_tid);
                                renamed_count += 1;
                            }
                        } else {
                            record.set_unmapped();
                            unmapped_count += 1;
                        }
                    }
                }
            }

            writer.write(&record)?;

            if i % 1000 == 0 {
                progress.set_position(i as u64);
            }
        }

        progress.finish_with_message("BAM processing complete");

        // Sort and index the output BAM
        Command::new("samtools")
            .args(&["sort", "-o", &self.output_bam, &self.output_bam])
            .status()
            .context("Failed to sort output BAM")?;

        Command::new("samtools")
            .args(&["index", &self.output_bam])
            .status()
            .context("Failed to index output BAM")?;

        println!("Processing complete:");
        println!("  - Renamed {} contigs", renamed_count);
        println!("  - Marked {} reads as unmapped", unmapped_count);

        Ok(())
    }

    fn cleanup(&self) -> Result<()> {
        for temp_file in [
            &self.temp_reference_sq_header,
            &self.temp_original_non_sq_header,
            &self.temp_combined_header,
        ] {
            if Path::new(temp_file).exists() {
                fs::remove_file(temp_file)
                    .with_context(|| format!("Failed to remove temporary file: {}", temp_file))?;
            }
        }
        Ok(())
    }
}
