use super::platform_inference::{PlatformInference, SequencingPlatform};
use crate::callable_loci::detect_aligner;
use crate::types::ReferenceGenome;
use crate::utils::bam_reader::BamReaderFactory;
use crate::utils::progress_manager::ProgressManager;
use rust_htslib::bam::Read;
use std::collections::HashMap;

pub struct BamStats {
    read_count: usize,
    total_read_length: usize,
    paired_reads: usize,
    total_insert_size: i64,
    paired_count: usize,
    pub(crate) max_samples: usize,
    length_distribution: HashMap<usize, usize>,
    insert_size_distribution: HashMap<i64, usize>,
    aligner: String,
    reference_build: String,
    flow_cells: HashMap<String, usize>,
    instruments: HashMap<String, usize>,
    platform_counts: HashMap<SequencingPlatform, usize>,
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
            aligner: String::new(),
            reference_build: String::new(),
            flow_cells: HashMap::new(),
            instruments: HashMap::new(),
            platform_counts: HashMap::new(),
        }
    }

    pub fn collect_stats(
        &mut self,
        bam_path: &str,
        reference_path: Option<&str>,
        progress_mgr: &ProgressManager,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut bam = BamReaderFactory::open(bam_path, reference_path)?;
        let header = bam.header().clone();

        self.aligner = detect_aligner(&header);
        self.reference_build = ReferenceGenome::from_header(&header)
            .map(|g| g.name().to_string())
            .unwrap_or_else(|| "Unknown".to_string());

        let progress = progress_mgr.add_spinner("Collecting BAM statistics...");

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

                // Extract flow cell and instrument from read name
                if let Ok(qname) = std::str::from_utf8(record.qname()) {
                    let platform = PlatformInference::detect_platform_from_qname(qname);
                    *self.platform_counts.entry(platform.clone()).or_insert(0) += 1;

                    match platform {
                        SequencingPlatform::Illumina => {
                            if let Some((instrument, flow_cell)) =
                                PlatformInference::parse_illumina_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                                *self.flow_cells.entry(flow_cell.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::PacBio => {
                            if let Some(instrument) =
                                PlatformInference::parse_pacbio_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::Nanopore => {
                            if let Some(instrument) =
                                PlatformInference::parse_nanopore_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::MGI => {
                            if let Some((instrument, flow_cell)) =
                                PlatformInference::parse_mgi_read_name(qname)
                            {
                                *self.instruments.entry(instrument.to_string()).or_insert(0) += 1;
                                *self.flow_cells.entry(flow_cell.to_string()).or_insert(0) += 1;
                            }
                        }
                        SequencingPlatform::Unknown => {}
                    }
                }

                if record.is_paired() {
                    self.paired_reads += 1;

                    if record.is_proper_pair() && record.is_first_in_template() {
                        let insert_size = record.insert_size().abs();
                        if insert_size > 0 {
                            *self
                                .insert_size_distribution
                                .entry(insert_size)
                                .or_insert(0) += 1;
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
            let modal_length = self
                .length_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(length, _)| *length)
                .unwrap_or(0);

            stats.insert("average_read_length".to_string(), modal_length as f64);
            stats.insert(
                "paired_percentage".to_string(),
                (self.paired_reads as f64 / self.read_count as f64) * 100.0,
            );
        }

        if self.paired_count > 0 {
            let modal_insert_size = self
                .insert_size_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(size, _)| *size)
                .unwrap_or(0);

            stats.insert("average_insert_size".to_string(), modal_insert_size as f64);
            stats.insert(
                "proper_pair_percentage".to_string(),
                (self.paired_count as f64 * 2.0 / self.paired_reads as f64) * 100.0,
            );
        }

        stats
    }

    pub fn aligner(&self) -> &str {
        &self.aligner
    }

    pub fn reference_build(&self) -> &str {
        &self.reference_build
    }

    pub fn average_read_length(&self) -> usize {
        if self.read_count > 0 {
            self.total_read_length / self.read_count
        } else {
            0
        }
    }

    pub fn modal_read_length(&self) -> usize {
        if self.read_count > 0 {
            self.length_distribution
                .iter()
                .max_by_key(|(_, &count)| count)
                .map(|(length, _)| *length)
                .unwrap_or(0)
        } else {
            0
        }
    }

    pub fn flow_cells(&self) -> &HashMap<String, usize> {
        &self.flow_cells
    }

    pub fn instruments(&self) -> &HashMap<String, usize> {
        &self.instruments
    }

    pub fn get_primary_platform(&self) -> SequencingPlatform {
        self.platform_counts
            .iter()
            .max_by_key(|(_, &count)| count)
            .map(|(platform, _)| platform.clone())
            .unwrap_or(SequencingPlatform::Unknown)
    }

    pub fn infer_platform(&self) -> String {
        let primary_platform = self.get_primary_platform();
        PlatformInference::infer_specific_platform(&primary_platform, &self.instruments)
    }

    pub fn get_instrument_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .instruments
            .iter()
            .map(|(id, &count)| {
                let percentage = (count as f64 / total_reads) * 100.0;
                (id.clone(), count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }

    pub fn get_flow_cell_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .flow_cells
            .iter()
            .map(|(id, &count)| {
                let percentage = (count as f64 / total_reads) * 100.0;
                (id.clone(), count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }

    pub fn get_platform_summary(&self) -> Vec<(String, usize, f64)> {
        let total_reads = self.read_count as f64;
        let mut summary: Vec<_> = self
            .platform_counts
            .iter()
            .map(|(platform, &count)| {
                let platform_name = match platform {
                    SequencingPlatform::Illumina => "Illumina",
                    SequencingPlatform::PacBio => "PacBio",
                    SequencingPlatform::Nanopore => "Oxford Nanopore",
                    SequencingPlatform::MGI => "MGI DNBseq",
                    SequencingPlatform::Unknown => "Unknown",
                }
                .to_string();
                let percentage = (count as f64 / total_reads) * 100.0;
                (platform_name, count, percentage)
            })
            .collect();

        summary.sort_by(|a, b| b.1.cmp(&a.1));
        summary
    }
}
