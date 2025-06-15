use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use rust_htslib::bam;
use crate::callable_loci::options::CallableOptions;

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
}

impl ContigProfiler {
    pub fn new(
        name: String,
        length: usize,
        multi_progress: &MultiProgress,
        options: CallableOptions,
    ) -> Self {
        let progress_bar = multi_progress.add(ProgressBar::new(length as u64));
        progress_bar.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-"));
        progress_bar.set_message(format!("Processing {}", name));

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

            self.n_reads += 1;

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
}