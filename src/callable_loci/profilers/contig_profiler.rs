use std::collections::HashSet;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::bam;
use crate::callable_loci::options::CallableOptions;
use crate::callable_loci::types::ContigStateCounts;
use crate::utils::progress_manager::ProgressManager;

pub struct ContigProfiler {
    pub name: String,
    pub length: usize,
    pub n_covered_bases: u64,
    pub summed_coverage: u64,
    pub summed_baseq: u64,
    pub summed_mapq: u64,
    pub quality_bases: u64,
    pub n_reads: u32,
    pub n_selected_reads: u32,
    pub progress_bar: ProgressBar,
    options: CallableOptions,
    seen_read_names: HashSet<Vec<u8>>,
}

impl ContigProfiler {
    pub fn new(
        name: String,
        length: usize,
        progress_mgr: &ProgressManager,
        options: CallableOptions,
    ) -> Self {
        let progress_bar = progress_mgr.add_contig_bar(&name, length);

        Self {
            name,
            length,
            n_covered_bases: 0,
            summed_coverage: 0,
            summed_baseq: 0,
            summed_mapq: 0,
            quality_bases: 0,
            n_reads: 0,
            n_selected_reads: 0,
            progress_bar,
            options,
            seen_read_names: HashSet::new(),
        }
    }

    pub fn process_position(
        &mut self,
        raw_depth: u32,
        alignments: bam::pileup::Alignments,
    ) {
        let mut qc_depth = 0;

        for aln in alignments {
            let record = aln.record();
            let mapq = record.mapq();

            // Count unique reads by tracking read names
            let read_name = record.qname().to_vec();
            if self.seen_read_names.insert(read_name) {
                self.n_reads += 1;
            }

            if mapq >= self.options.min_mapping_quality {
                if let Some(qpos) = aln.qpos() {
                    if let Some(&qual) = record.qual().get(qpos) {
                        if qual >= self.options.min_base_quality {
                            qc_depth += 1;
                            self.summed_baseq += qual as u64;
                            self.quality_bases += 1;
                        }
                    }
                }
                self.summed_mapq += mapq as u64;
                self.n_selected_reads += 1;
            }
        }

        if raw_depth > 0 {
            self.n_covered_bases += 1;
            self.summed_coverage += raw_depth as u64;
        }
    }

    pub fn finish(&self) {
        self.progress_bar.finish_with_message("Complete");
    }

    /// Calculate coverage statistics for this contig
    pub fn get_coverage_stats(&self) -> crate::callable_loci::types::ContigCoverageStats {
        use crate::callable_loci::types::{ContigCoverageStats, ContigStateCounts};

        let average_depth = if self.n_covered_bases > 0 {
            self.summed_coverage as f64 / self.n_covered_bases as f64
        } else {
            0.0
        };

        let coverage_percent = if self.length > 0 {
            (self.n_covered_bases as f64 / self.length as f64) * 100.0
        } else {
            0.0
        };

        ContigCoverageStats {
            unique_reads: self.n_reads as u64,
            state_counts: ContigStateCounts {
                reference_n: 0, // Will be populated from CallableProfiler
                no_coverage: 0,
                low_coverage: 0,
                excessive_coverage: 0,
                poor_mapping_quality: 0,
                callable: 0,
            },
            coverage_percent,
            average_depth,
            covered_bases: self.n_covered_bases,
            total_bases: self.length as u64,
        }
    }

    /// Calculate quality statistics for this contig
    pub fn get_quality_stats(&self) -> crate::callable_loci::types::ContigQualityStats {
        use crate::callable_loci::types::ContigQualityStats;

        let average_mapq = if self.quality_bases > 0 {
            self.summed_mapq as f64 / self.quality_bases as f64
        } else {
            0.0
        };

        let average_baseq = if self.quality_bases > 0 {
            self.summed_baseq as f64 / self.quality_bases as f64
        } else {
            0.0
        };

        // Estimate Q30 percentage based on average base quality
        // This is an approximation - for exact calculation you'd need to track Q30 bases separately
        let q30_percentage = if self.quality_bases > 0 {
            if average_baseq >= 30.0 {
                100.0
            } else if average_baseq < 20.0 {
                0.0
            } else {
                // Linear interpolation between Q20 and Q30
                ((average_baseq - 20.0) / 10.0) * 100.0
            }
        } else {
            0.0
        };

        ContigQualityStats {
            average_mapq,
            average_baseq,
            q30_percentage,
        }
    }
}