use anyhow::{anyhow, Context, Result};
use crossbeam_channel::{bounded, unbounded};
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::{self, Header, Read, Reader, Writer};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::Path;
use std::process::Command;
use std::thread;
use std::time::Duration;


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
    temp_fixed_bam: String,
    reference_genome: ReferenceGenome,
}

impl BamFixer {
    pub fn new(
        reference_file: String,
        input_bam: String,
        output_bam: String,
        keep_temp: bool,
    ) -> Result<Self> {
        // Check for samtools availability first
        crate::utils::external_tools::check_samtools()?;

        let reference_genome = Self::detect_reference_genome(&input_bam)
            .context("Failed to detect reference genome type")?;

        let temp_dir = std::env::temp_dir();
        let base_name = Path::new(&input_bam)
            .file_name()
            .unwrap()
            .to_string_lossy()
            .to_string();

        Ok(BamFixer {
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
            temp_fixed_bam: temp_dir
                .join(format!("{}.fixed.bam", base_name))
                .to_string_lossy()
                .to_string(),
            reference_genome,
        })
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

    fn extract_reference_sq_headers(&self) -> Result<()> {
        self.generate_reference_header()?;
        Ok(())
    }

    pub fn run(&self) -> Result<()> {
        let progress = ProgressBar::new_spinner();
        progress.enable_steady_tick(Duration::from_millis(100));
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );

        // Process headers
        progress.set_message("Processing BAM headers...");
        self.extract_non_sq_headers()?;
        self.extract_reference_sq_headers()?;
        self.combine_headers()?;

        // Fix BAM
        progress.set_message("Fixing BAM file...");
        self.fix_bam()?;

        // Sort BAM
        self.run_samtools_with_progress(
            &["sort", "-o", &self.output_bam, &self.temp_fixed_bam],
            "Sorting BAM file (this may take several minutes for large files)...".to_string(),
        )?;

        // Index BAM
        self.run_samtools_with_progress(
            &["index", &self.output_bam],
            "Indexing BAM file...".to_string(),
        )?;

        // Cleanup
        if !self.keep_temp {
            progress.set_message("Cleaning up temporary files...");
            self.cleanup()?;
            progress.finish_with_message("Cleanup complete");
        }

        progress.finish_with_message("BAM file processing complete!");
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
        let reader = Reader::from_path(&self.input_bam).context("Failed to open input BAM")?;

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

        // Create valid reference name mappings
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

        // Setup BAM reader
        let reader = Reader::from_path(&self.input_bam)?;
        let input_header = reader.header().clone();

        // Setup BAM writer
        let mut writer = Writer::from_path(
            &self.temp_fixed_bam,
            &header,
            bam::Format::Bam
        )?;

        // Get reference genome prefix
        let prefix = self.reference_genome.prefix().to_string();

        // Setup progress bar
        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] Processing BAM records...")
                .unwrap(),
        );

        // Create channels for parallel processing
        let (record_tx, record_rx) = bounded::<bam::Record>(10000);
        let (stat_tx, stat_rx) = unbounded();

        // Clone necessary data for the reader thread
        let input_bam = self.input_bam.clone();

        // Spawn reader thread
        let reader_thread = thread::spawn(move || {
            let mut reader = Reader::from_path(&input_bam).unwrap();
            for record_result in reader.records() {
                match record_result {
                    Ok(record) => {
                        if record_tx.send(record).is_err() {
                            break;
                        }
                    }
                    Err(_) => break,
                }
            }
        });

        // Extract reference names before spawning the thread
        let ref_names: Vec<(u32, String)> = (0..input_header.target_names().len() as u32)
            .filter_map(|tid| {
                let name = input_header.tid2name(tid);
                if name.is_empty() {
                    None
                } else {
                    Some((tid, String::from_utf8_lossy(name).to_string()))
                }
            })
            .collect();

        // Clone necessary data for the processing thread
        let ref_name_to_tid = ref_name_to_tid.clone();
        let prefix = prefix.clone();

        // Process records in the main thread
        let processing_thread = thread::spawn(move || {
            let mut renamed_count = 0;
            let mut unmapped_count = 0;

            while let Ok(mut record) = record_rx.recv() {
                if !record.is_unmapped() {
                    let tid = record.tid() as u32;
                    if let Some(ref_name) = ref_names.iter()
                        .find(|(t, _)| *t == tid)
                        .map(|(_, name)| name.as_str()) {
                        if ref_name.starts_with(&prefix) {
                            let stripped_name = ref_name
                                .trim_start_matches(&prefix)
                                .to_string();

                            if let Some(&new_tid) = ref_name_to_tid.get(&stripped_name) {
                                record.set_tid(new_tid);
                                renamed_count += 1;
                            } else {
                                record.set_unmapped();
                                unmapped_count += 1;
                            }
                        }
                    }
                }

                writer.write(&record).unwrap();
            }

            stat_tx.send((renamed_count, unmapped_count)).unwrap();
        });

        // Wait for threads to complete
        reader_thread.join().unwrap();
        processing_thread.join().unwrap();

        // Get statistics
        let (renamed_count, unmapped_count) = stat_rx.recv()?;

        progress.finish_with_message("BAM processing complete");

        println!("Processing complete:");
        println!("  - Renamed {} contigs", renamed_count);
        println!("  - Marked {} reads as unmapped", unmapped_count);

        Ok(())
    }


    fn run_samtools_with_progress(&self, args: &[&str], message: String) -> Result<()> {
        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        progress.set_message(message.clone());

        let output = Command::new("samtools")
            .args(args)
            .output()
            .context(format!("Failed to execute samtools {}", args[0]))?;

        if !output.status.success() {
            progress.finish_with_message(format!("Error in samtools {}", args[0]));
            return Err(anyhow!(
                "samtools {} failed: {}",
                args[0],
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        progress.finish_with_message(format!("Completed: {}", message));
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
