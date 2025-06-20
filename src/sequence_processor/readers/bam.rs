use crate::sequence_processor::{core::*, threading::*};
use anyhow::Result;
use indicatif::ProgressBar;
use rust_htslib::bam::{self, Read};
use std::path::Path;

pub struct BamReader {
    reader: bam::Reader,
}

impl BamReader {
    pub fn new(path: &Path, reference: Option<&str>) -> Result<Self> {
        let mut reader = bam::Reader::from_path(path)?;

        if path.extension().map_or(false, |ext| ext == "cram") {
            if let Some(ref_path) = reference {
                reader.set_reference(ref_path)?;
            }
        }

        Ok(Self {
            reader
        })
    }

    fn create_sequence_from_record(&self, record: &bam::Record) -> Sequence {
        let seq = record.seq().as_bytes();
        let qual = record.qual().to_vec();
        let read_name = String::from_utf8_lossy(record.qname()).into_owned();

        Sequence {
            data: seq.to_vec(),
            id: Some(read_name),
            quality: Some(qual),
            metadata: SequenceMetadata {
                is_mapped: !record.is_unmapped(),
                chromosome: if record.tid() >= 0 {
                    let name = self.reader.header().tid2name(record.tid() as u32);
                    Some(String::from_utf8_lossy(name).into_owned())
                } else {
                    None
                },
                position: Some(record.pos() as u64),
                mapping_quality: Some(record.mapq()),
            },
        }
    }
}

impl SequenceReader for BamReader {
    fn read_sequences_single_thread<P: SequenceProcessor>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
    ) -> Result<ProcessingStats> {
        let mut stats = ProcessingStats::default();
        let mut record = bam::Record::new();

        while let Some(result) = self.reader.read(&mut record) {
            match result {
                Ok(()) => {
                    let seq = record.seq().as_bytes();
                    if seq.len() >= processor.get_min_length() {
                        let sequence = self.create_sequence_from_record(&record);
                        if let Err(e) = processor.process_sequence(&sequence) {
                            eprintln!("Error processing sequence: {}", e);
                            stats.errors += 1;
                        } else {
                            stats.processed += 1;
                        }
                    } else {
                        stats.too_short += 1;
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Error reading record: {}", e);
                    stats.errors += 1;
                }
            }

            if stats.processed % 1000 == 0 {
                processor.update_progress(&stats);
            }
        }

        Ok(stats)
    }

    fn read_sequences_with_threads<P: SequenceProcessor + Clone + 'static>(
        &mut self,
        processor: &mut P,
        _progress: &ProgressBar,
        num_threads: usize,
    ) -> Result<ProcessingStats> {
        if num_threads <= 1 || !processor.supports_parallel() {
            return self.read_sequences_single_thread(processor, _progress);
        }

        let pool = ThreadPool::new(processor.clone(), num_threads)?;
        let mut stats = ProcessingStats::default();
        let mut record = bam::Record::new();

        while let Some(result) = self.reader.read(&mut record) {
            match result {
                Ok(()) => {
                    let seq = record.seq().as_bytes();
                    if seq.len() >= processor.get_min_length() {
                        let sequence = self.create_sequence_from_record(&record);
                        pool.send(sequence)?;
                    } else {
                        stats.too_short += 1;
                    }
                }
                Err(e) => {
                    eprintln!("Warning: Error reading record: {}", e);
                    stats.errors += 1;
                }
            }
        }

        let (thread_stats, processors) = pool.finish()?;
        stats.processed += thread_stats.processed;
        stats.errors += thread_stats.errors;

        merge_processors(processors, processor)?;

        Ok(stats)
    }
}