use crate::sequence_processor::{core::*, threading::*};
use anyhow::Result;
use bio::io::fastq::{self, FastqRead};
use indicatif::ProgressBar;
use niffler::get_reader;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub struct FastqReader {
    reader: fastq::Reader<BufReader<Box<dyn std::io::Read>>>,
}

impl FastqReader {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let (inner_reader, _compression) = get_reader(Box::new(file))?;
        Ok(Self {
            reader: fastq::Reader::new(Box::new(BufReader::with_capacity(
                16 * 1024 * 1024,
                inner_reader,
            ))),
        })
    }

    fn create_sequence_from_record(&self, record: &fastq::Record) -> Sequence {
        Sequence {
            data: record.seq().to_vec(),
            id: Some(record.id().to_string()),
            quality: Some(record.qual().to_vec()),
            metadata: SequenceMetadata::default(),
        }
    }
}

impl SequenceReader for FastqReader {
    fn read_sequences_single_thread<P: SequenceProcessor>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
    ) -> Result<ProcessingStats> {
        let mut stats = ProcessingStats::default();
        let mut record = fastq::Record::new();

        while let Ok(()) = self.reader.read(&mut record) {
            if record.seq().len() >= processor.get_min_length() {
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
        let mut batch = Vec::with_capacity(1000);
        let mut record = fastq::Record::new();

        while self.reader.read(&mut record).is_ok() {
            if record.id().is_empty() {
                break; // End of file or invalid record
            }

            if record.seq().len() >= processor.get_min_length() {
                let sequence = self.create_sequence_from_record(&record);
                batch.push(sequence);

                if batch.len() >= 1000 {
                    for seq in batch.drain(..) {
                        if let Err(e) = pool.send(seq) {
                            eprintln!("Error sending sequence to worker: {}", e);
                            stats.errors += 1;
                        }
                    }
                }
            } else {
                stats.too_short += 1;
            }
        }

        // Process remaining sequences in batch
        for seq in batch.drain(..) {
            if let Err(e) = pool.send(seq) {
                eprintln!("Error sending sequence to worker: {}", e);
                stats.errors += 1;
            }
        }

        let (thread_stats, processors) = pool.finish()?;
        stats.processed += thread_stats.processed;
        stats.errors += thread_stats.errors;

        merge_processors(processors, processor)?;

        Ok(stats)
    }
}