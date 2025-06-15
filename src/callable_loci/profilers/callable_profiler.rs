use crate::callable_loci::utils::histogram_plotter;
use crate::callable_loci::utils::histogram_plotter::CoverageRange;
use std::collections::HashMap;
use std::fs::File;
use std::io::Result as IoResult;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledState {
    REF_N,
    CALLABLE,
    NO_COVERAGE,
    LOW_COVERAGE,
    EXCESSIVE_COVERAGE,
    POOR_MAPPING_QUALITY,
}

pub struct CallableProfiler {
    counts: [u64; 6],
    current_state: Option<(String, u64, u64, CalledState)>,
    bed_writer: BufWriter<File>,
    coverage_ranges: Vec<CoverageRange>,
    output_dir: PathBuf,
    largest_contig_length: u32,
    contig_counts: HashMap<String, [u64; 6]>,
}

impl CallableProfiler {
    pub fn new(bed_file: &str, largest_contig_length: u32) -> IoResult<Self> {
        let output_dir = PathBuf::from(bed_file)
            .parent()
            .unwrap_or(&PathBuf::from("."))
            .to_path_buf();

        Ok(Self {
            counts: [0; 6],
            current_state: None,
            bed_writer: BufWriter::new(File::create(bed_file)?),
            coverage_ranges: Vec::new(),
            output_dir,
            largest_contig_length,
            contig_counts: Default::default(),
        })
    }

    fn write_state(&mut self) -> IoResult<()> {
        if let Some((contig, start, end, state)) = &self.current_state {
            writeln!(
                self.bed_writer,
                "{}\t{}\t{}\t{:?}",
                contig, start, end, state
            )?;
        }
        Ok(())
    }

    pub(crate) fn finish_contig(&mut self, contig_name: &str, contig_length: u32) -> IoResult<()> {
        self.write_state()?;

        if !self.coverage_ranges.is_empty() {
            let histogram_path = histogram_plotter::generate_histogram(
                std::mem::take(&mut self.coverage_ranges),
                &format!("{}_coverage", contig_name),
                contig_name,
                contig_length,
                if contig_name == "chrM" {
                    contig_length
                } else {
                    self.largest_contig_length
                },
            )?;

            let final_path = self
                .output_dir
                .join(format!("{}_coverage.svg", contig_name));
            std::fs::copy(histogram_path, final_path)?;
        }

        Ok(())
    }

    pub(crate) fn process_position(
        &mut self,
        contig: &str,
        pos: u32,
        depth: u32,
        is_low_mapq: bool,
        state: CalledState,
    ) -> IoResult<()> {
        self.counts[state as usize] += 1;

        self.contig_counts
            .entry(contig.to_string())
            .or_insert([0; 6])[state as usize] += 1;

        match self.coverage_ranges.last_mut() {
            Some(range) if range.can_merge(pos, depth, is_low_mapq) => {
                range.extend(pos);
            }
            _ => {
                self.coverage_ranges
                    .push(CoverageRange::new(pos, depth, is_low_mapq));
            }
        }

        // Handle state tracking (existing code)
        match &mut self.current_state {
            Some((cur_contig, _start, end, cur_state)) => {
                if contig != cur_contig || state != *cur_state || pos as u64 != *end + 1 {
                    self.write_state()?;
                    self.current_state = Some((contig.to_string(), pos as u64, pos as u64, state));
                } else {
                    *end = pos as u64;
                }
            }
            None => {
                self.current_state = Some((contig.to_string(), pos as u64, pos as u64, state));
            }
        }

        Ok(())
    }

    pub fn get_contig_counts(&self, contig: &str) -> [u64; 6] {
        self.contig_counts.get(contig).copied().unwrap_or([0; 6])
    }
}
