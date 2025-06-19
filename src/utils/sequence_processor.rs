pub mod core {
    use anyhow::Result;
    use indicatif::ProgressBar;

    #[derive(Debug, Clone)]
    pub struct Sequence {
        pub data: Vec<u8>,
        pub id: Option<String>,
        pub quality: Option<Vec<u8>>,
        pub metadata: SequenceMetadata,
    }

    #[derive(Debug, Clone, Default)]
    pub struct SequenceMetadata {
        pub is_mapped: bool,
        pub chromosome: Option<String>,
        pub position: Option<u64>,
    }

    #[derive(Debug, Default)]
    pub struct ProcessingStats {
        pub processed: u64,
        pub skipped: u64,
        pub errors: u64,
        pub too_short: u64,
        // Add other common statistics
    }

    pub trait SequenceProcessor: Send + Clone + 'static {
        fn process_sequence(&mut self, sequence: &Sequence) -> Result<()>;
        fn get_min_length(&self) -> usize;
        fn update_progress(&mut self, stats: &ProcessingStats);
        fn finalize(&mut self) -> Result<()> {
            Ok(())
        }
        fn supports_parallel(&self) -> bool {
            false
        }
        fn merge_processor(&mut self, other: &Self) -> Result<()> {
            Ok(())
        }
    }

    pub trait SequenceReader {
        fn read_sequences_single_thread<P: SequenceProcessor>(
            &mut self,
            processor: &mut P,
            progress: &ProgressBar,
        ) -> Result<ProcessingStats>;

        fn read_sequences<P: SequenceProcessor + Clone + 'static>(
            &mut self,
            processor: &mut P,
            progress: &ProgressBar,
        ) -> Result<ProcessingStats> {
            self.read_sequences_with_threads(processor, progress, 1)
        }

        fn read_sequences_with_threads<P: SequenceProcessor + Clone + 'static>(
            &mut self,
            processor: &mut P,
            progress: &ProgressBar,
            num_threads: usize,
        ) -> Result<ProcessingStats>;
    }
}

pub mod readers {
    use super::core::*;
    use crate::utils::progress_bar_builder::ProgressBarBuilder;
    use anyhow::Result;
    use bio::io::fastq;
    use bio::io::fastq::FastqRead;
    use crossbeam_channel::bounded;
    use indicatif::ProgressBar;
    use niffler::get_reader;
    use rust_htslib::bam;
    use rust_htslib::bam::Read;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;
    use std::thread;

    pub struct BamReader {
        reader: bam::Reader,
        chunk_size: usize,
    }

    pub struct FastqReader {
        reader: fastq::Reader<BufReader<Box<dyn std::io::Read>>>,
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
                reader,
                chunk_size: 10000,
            })
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
                            let sequence = Sequence {
                                data: seq.to_vec(),
                                id: None,
                                quality: None,
                                metadata: SequenceMetadata {
                                    is_mapped: !record.is_unmapped(),
                                    chromosome: if record.tid() >= 0 {
                                        let name =
                                            self.reader.header().tid2name(record.tid() as u32);
                                        Some(String::from_utf8_lossy(name).into_owned())
                                    } else {
                                        None
                                    },
                                    position: Some(record.pos() as u64),
                                },
                            };
                            processor.process_sequence(&sequence)?;
                            stats.processed += 1;
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

        fn read_sequences_with_threads<P: SequenceProcessor + Clone>(
            &mut self,
            processor: &mut P,
            progress: &ProgressBar,
            num_threads: usize,
        ) -> Result<ProcessingStats> {
            if num_threads <= 1 || !processor.supports_parallel() {
                return self.read_sequences_single_thread(processor, progress);
            }

            let (tx, rx) = bounded(num_threads * 2);
            let mut handles = vec![];
            let mut processors = vec![];

            // Spawn worker threads
            for _ in 0..num_threads {
                let rx = rx.clone();
                let mut worker_processor = processor.clone();
                let handle = thread::spawn(move || {
                    let mut local_stats = ProcessingStats::default();
                    while let Ok(sequence) = rx.recv() {
                        if let Err(e) = worker_processor.process_sequence(&sequence) {
                            eprintln!("Error processing sequence: {}", e);
                            local_stats.errors += 1;
                        } else {
                            local_stats.processed += 1;
                        }
                    }
                    (worker_processor, local_stats)
                });
                handles.push(handle);
            }
            // Read and send sequences
            let mut record = bam::Record::new();
            let mut stats = ProcessingStats::default();

            while let Some(result) = self.reader.read(&mut record) {
                match result {
                    Ok(()) => {
                        let seq = record.seq().as_bytes();
                        if seq.len() >= processor.get_min_length() {
                            let sequence = Sequence {
                                data: seq.to_vec(),
                                id: None,
                                quality: None,
                                metadata: SequenceMetadata {
                                    is_mapped: !record.is_unmapped(),
                                    chromosome: if record.tid() >= 0 {
                                        let name =
                                            self.reader.header().tid2name(record.tid() as u32);
                                        Some(String::from_utf8_lossy(name).into_owned())
                                    } else {
                                        None
                                    },
                                    position: Some(record.pos() as u64),
                                },
                            };
                            tx.send(sequence)?;
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

            drop(tx); // Close channel to signal workers to finish

            // Progress bar for collecting worker results
            let collect_progress = ProgressBarBuilder::new("Collecting worker results")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                .with_tick()
                .build()?;

            let num_handles = handles.len();
            collect_progress.set_length(num_handles as u64);

            // Collect results
            for (idx, handle) in handles.into_iter().enumerate() {
                collect_progress.set_message(format!(
                    "Collecting worker {} of {}",
                    idx + 1,
                    num_handles
                ));
                let (worker_processor, worker_stats) = handle.join().unwrap();
                processors.push(worker_processor);
                stats.processed += worker_stats.processed;
                stats.errors += worker_stats.errors;

                if stats.processed % 1000 == 0 {
                    processor.update_progress(&stats);
                }
                collect_progress.inc(1);
            }
            collect_progress.finish_with_message(format!("Collected {} workers", processors.len()));

            // Progress bar for merging results
            let merge_progress = ProgressBarBuilder::new("Merging processors")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                .with_tick()
                .build()?;

            merge_progress.set_length(processors.len() as u64);

            // Merge results from all processors
            for (idx, worker_processor) in processors.iter().enumerate() {
                merge_progress.set_message(format!(
                    "Merging processor {} of {}",
                    idx + 1,
                    processors.len()
                ));
                processor.merge_processor(worker_processor)?;
                merge_progress.inc(1);
            }
            merge_progress.finish_with_message(format!("Merged {} processors", processors.len()));

            Ok(stats)
        }
    }

    impl FastqReader {
        pub fn new(path: &Path) -> Result<Self> {
            let file = File::open(path)?;
            let (inner_reader, _compression) = get_reader(Box::new(file))?;
            // Increase buffer size to 16MB
            Ok(Self {
                reader: fastq::Reader::new(Box::new(BufReader::with_capacity(
                    16 * 1024 * 1024,
                    inner_reader,
                ))),
            })
        }

        // Private helper methods that will support the main implementation
        fn spawn_worker_threads<P: SequenceProcessor>(
            num_threads: usize,
            rx: &crossbeam_channel::Receiver<Sequence>,
            processor: &P,
            progress: &ProgressBar,
        ) -> Vec<thread::JoinHandle<(P, ProcessingStats)>> {
            let mut handles = vec![];

            for thread_id in 0..num_threads {
                let rx = rx.clone();  // Now we're cloning the reference which is fine
                let mut worker_processor = processor.clone();
                let worker_progress = progress.clone();

                let handle = thread::spawn(move || {
                    let mut local_stats = ProcessingStats::default();
                    while let Ok(sequence) = rx.recv() {
                        if let Err(e) = worker_processor.process_sequence(&sequence) {
                            eprintln!("Thread {}: Error processing sequence: {}", thread_id, e);
                            local_stats.errors += 1;
                        } else {
                            local_stats.processed += 1;
                            if local_stats.processed % 10000 == 0 {
                                worker_progress.inc(10000);
                            }
                        }
                    }
                    eprintln!("Worker {} finished", thread_id);
                    (worker_processor, local_stats)
                });
                handles.push(handle);
            }
            handles
        }

        fn process_batch(
            batch: &mut Vec<Sequence>,
            tx: &crossbeam_channel::Sender<Sequence>,
            record_count: u64,
            read_progress: &ProgressBar,
        ) -> Result<bool> {
            for seq in batch.drain(..) {
                if tx.send(seq).is_err() {
                    read_progress.finish_and_clear();
                    eprintln!("Channel send error after {} records", record_count);
                    return Ok(false);
                }
            }
            Ok(true)
        }

        fn collect_worker_results<P: SequenceProcessor>(
            handles: Vec<thread::JoinHandle<(P, ProcessingStats)>>,
            stats: &mut ProcessingStats,
            collect_progress: &ProgressBar,
        ) -> Result<Vec<P>> {
            let mut results = Vec::new();
            let num_handles = handles.len();
            collect_progress.set_length(num_handles as u64);

            for (idx, handle) in handles.into_iter().enumerate() {
                eprintln!("Waiting for worker thread {} of {}...", idx + 1, num_handles);
                collect_progress.set_message(format!(
                    "Collecting worker {} of {} ({} processed so far)",
                    idx + 1, num_handles, stats.processed
                ));

                match handle.join() {
                    Ok((processor, local_stats)) => {
                        eprintln!("Worker {} completed with {} sequences processed",
                                  idx + 1, local_stats.processed);
                        stats.processed += local_stats.processed;
                        stats.errors += local_stats.errors;
                        results.push(processor);
                    }
                    Err(e) => {
                        eprintln!("Worker thread {} panicked: {:?}", idx + 1, e);
                        stats.errors += 1;
                    }
                }
                collect_progress.inc(1);
            }
            collect_progress.finish_with_message(format!(
                "Collection complete - processed {} sequences",
                stats.processed
            ));
            Ok(results)
        }

        fn merge_worker_results<P: SequenceProcessor>(
            results: &[P],
            processor: &mut P,
            merge_progress: &ProgressBar,
        ) -> Result<()> {
            merge_progress.set_length(results.len() as u64);

            for (idx, worker_processor) in results.iter().enumerate() {
                let start = std::time::Instant::now();
                merge_progress.set_message(format!("Merging worker {} of {}", idx + 1, results.len()));
                eprintln!("Starting merge of worker {} of {}", idx + 1, results.len());
                processor.merge_processor(worker_processor)?;
                let duration = start.elapsed();
                eprintln!("Completed merge of worker {} in {:?}", idx + 1, duration);
                merge_progress.set_message(format!("Merged worker {} of {} (took {:?})",
                                                   idx + 1, results.len(), duration));
                merge_progress.inc(1);
            }

            merge_progress.finish_with_message(format!("Merged {} workers", results.len()));
            Ok(())
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
                    let sequence = Sequence {
                        data: record.seq().to_vec(),
                        id: Some(record.id().to_string()),
                        quality: Some(record.qual().to_vec()),
                        metadata: SequenceMetadata::default(),
                    };
                    processor.process_sequence(&sequence)?;
                    stats.processed += 1;
                } else {
                    stats.too_short += 1;
                }

                if stats.processed % 1000 == 0 {
                    processor.update_progress(&stats);
                }
            }

            Ok(stats)
        }

        fn read_sequences_with_threads<P: SequenceProcessor>(
            &mut self,
            processor: &mut P,
            progress: &ProgressBar,
            num_threads: usize,
        ) -> Result<ProcessingStats> {
            if num_threads <= 1 || !processor.supports_parallel() {
                return self.read_sequences_single_thread(processor, progress);
            }

            let (tx, rx) = bounded(num_threads * 10000);
            let handles = Self::spawn_worker_threads(num_threads, &rx, processor, progress); // Pass reference to rx
            drop(rx);

            let mut stats = ProcessingStats::default();
            let mut batch = Vec::with_capacity(1000);
            let mut record = fastq::Record::new();
            let mut record_count = 0;

            let read_progress = ProgressBarBuilder::new("Reading records")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg} [{bar:40.cyan/blue}] {pos}/{len} ({per_sec}/s)\r")
                .with_tick()
                .build()?;
            read_progress.enable_steady_tick(std::time::Duration::from_secs(10));

            while let Ok(()) = self.reader.read(&mut record) {
                record_count += 1;
                if record_count % 100000 == 0 {
                    read_progress.set_message(format!("Processing {} bp sequence", record.seq().len()));
                    read_progress.set_position(record_count);
                }

                if record.seq().len() >= processor.get_min_length() {
                    batch.push(Sequence {
                        data: record.seq().to_vec(),
                        id: Some(record.id().to_string()),
                        quality: Some(record.qual().to_vec()),
                        metadata: SequenceMetadata::default(),
                    });

                    if batch.len() >= 1000 && !Self::process_batch(&mut batch, &tx, record_count, &read_progress)? {
                        return Ok(stats);
                    }
                } else {
                    stats.too_short += 1;
                }
            }

            read_progress.finish_and_clear();

            // Process remaining sequences
            if !batch.is_empty() && !Self::process_batch(&mut batch, &tx, record_count, &read_progress)? {
                return Ok(stats);
            }

            drop(tx);

            let collect_progress = ProgressBarBuilder::new("Collecting worker results")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                .with_tick()
                .build()?;

            let results = Self::collect_worker_results(handles, &mut stats, &collect_progress)?;

            let merge_progress = ProgressBarBuilder::new("Merging worker results")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                .with_tick()
                .build()?;

            Self::merge_worker_results(&results, processor, &merge_progress)?;

            Ok(stats)
        }
    }
}

