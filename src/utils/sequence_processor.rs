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
            let progress_update_interval = 10000; // Update progress every 10k sequences

            for thread_id in 0..num_threads {
                let rx = rx.clone();
                let mut worker_processor = processor.clone();
                let worker_progress = progress.clone();

                let handle = thread::spawn(move || {
                    let mut local_stats = ProcessingStats::default();

                    // Track processing rate
                    let start_time = std::time::Instant::now();
                    let mut last_update = start_time;

                    while let Ok(sequence) = rx.recv() {
                        if let Err(e) = worker_processor.process_sequence(&sequence) {
                            eprintln!("Thread {}: Error processing sequence: {}", thread_id, e);
                            local_stats.errors += 1;
                        } else {
                            local_stats.processed += 1;

                            // Update progress periodically
                            if local_stats.processed % progress_update_interval == 0 {
                                let now = std::time::Instant::now();
                                let elapsed = now.duration_since(last_update);
                                let rate = progress_update_interval as f64 / elapsed.as_secs_f64();

                                worker_progress.set_message(format!(
                                    "Worker {} | {} sequences ({:.0} seq/s)",
                                    thread_id, local_stats.processed, rate
                                ));

                                last_update = now;
                            }
                        }
                    }

                    // Final progress update for this worker
                    let total_time = start_time.elapsed();
                    let overall_rate = local_stats.processed as f64 / total_time.as_secs_f64();
                    eprintln!(
                        "Worker {} finished: {} sequences processed ({:.0} seq/s average)",
                        thread_id, local_stats.processed, overall_rate
                    );

                    (worker_processor, local_stats)
                });
                handles.push(handle);
            }
            handles
        }

        fn thousands_separator(n: u64) -> String {
            n.to_string()
                .as_bytes()
                .rchunks(3)
                .rev()
                .map(std::str::from_utf8)
                .collect::<Result<Vec<&str>, _>>()
                .unwrap()
                .join(",")
        }

        fn process_batch(
            batch: &mut Vec<Sequence>,
            tx: &crossbeam_channel::Sender<Sequence>,
            record_count: u64,
            read_progress: &ProgressBar,
        ) -> Result<bool> {
            // Update progress before starting batch processing
            read_progress.set_message(format!("Processing batch at record {}", record_count));

            // Send sequences in chunks to avoid excessive blocking
            const CHUNK_SIZE: usize = 100;
            for chunk in batch.drain(..).collect::<Vec<_>>().chunks(CHUNK_SIZE) {
                for seq in chunk {
                    match tx.try_send(seq.clone()) {
                        Ok(_) => {}
                        Err(crossbeam_channel::TrySendError::Full(_)) => {
                            eprintln!("Channel full at record {}. Active batch size: {}. Attempting retry...",
                                      record_count, batch.len());

                            // Try a few times with backoff before giving up
                            for retry in 1..=3 {
                                thread::sleep(std::time::Duration::from_millis(100 * retry));
                                match tx.try_send(seq.clone()) {
                                    Ok(_) => break,
                                    Err(crossbeam_channel::TrySendError::Full(_)) if retry == 3 => {
                                        read_progress.finish_and_clear();
                                        return Err(anyhow::anyhow!(
                        "Channel remained full after retries at record {}. Worker threads may be stuck.",
                        record_count
                    ));
                                    }
                                    _ => continue,
                                }
                            }
                        }
                        Err(e) => {
                            read_progress.finish_and_clear();
                            return Err(anyhow::anyhow!(
                                "Channel send error at record {}: {:?}",
                                record_count,
                                e
                            ));
                        }
                    }
                }
                // Update progress after each chunk
                read_progress.set_message(format!("Processed {} sequences", record_count));
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
                eprintln!(
                    "Waiting for worker thread {} of {}...",
                    idx + 1,
                    num_handles
                );
                collect_progress.set_message(format!(
                    "Collecting worker {} of {} ({} processed so far)",
                    idx + 1,
                    num_handles,
                    stats.processed
                ));

                match handle.join() {
                    Ok((processor, local_stats)) => {
                        eprintln!(
                            "Worker {} completed with {} sequences processed",
                            idx + 1,
                            local_stats.processed
                        );
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
                merge_progress.set_message(format!(
                    "Merging worker {} of {}",
                    idx + 1,
                    results.len()
                ));
                eprintln!("Starting merge of worker {} of {}", idx + 1, results.len());
                processor.merge_processor(worker_processor)?;
                let duration = start.elapsed();
                eprintln!("Completed merge of worker {} in {:?}", idx + 1, duration);
                merge_progress.set_message(format!(
                    "Merged worker {} of {} (took {:?})",
                    idx + 1,
                    results.len(),
                    duration
                ));
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

            // Use a smaller channel capacity to prevent memory issues
            let (tx, rx) = bounded(num_threads * 1000);
            let handles = Self::spawn_worker_threads(num_threads, &rx, processor, progress);

            let mut stats = ProcessingStats::default();
            let mut batch = Vec::with_capacity(1000);
            let mut record_count: u64 = 0;

            // Create a progress bar for reading phase
            let read_progress = ProgressBarBuilder::new("Reading FASTQ")
                .with_template(
                    "{spinner:.green} [{elapsed_precise}] {msg} ({records} records processed)",
                )
                .with_tick()
                .build()?;

            // Read and process sequences
            eprintln!("Starting sequence reading phase...");
            let mut record = fastq::Record::new();
            while self.reader.read(&mut record).is_ok() {
                if record.id().is_empty() {
                    break; // End of file or invalid record
                }

                record_count += 1;

                if record_count % 10000 == 0 {
                    read_progress.set_message(format!(
                        "Reading sequences ({} records processed)",
                        record_count
                    ));
                }

                if record.seq().len() >= processor.get_min_length() {
                    batch.push(Sequence {
                        data: record.seq().to_vec(),
                        id: Some(record.id().to_string()),
                        quality: Some(record.qual().to_vec()),
                        metadata: SequenceMetadata::default(),
                    });

                    if batch.len() >= 1000 {
                        if !Self::process_batch(&mut batch, &tx, record_count, &read_progress)? {
                            return Ok(stats);
                        }
                    }
                } else {
                    stats.too_short += 1;
                }
            }

            read_progress.finish_with_message(format!(
                "Read phase complete: processed {} sequences",
                record_count
            ));

            // Drop sender to signal workers no more data is coming
            drop(tx);

            // Create a new progress bar for worker completion phase
            let worker_progress = ProgressBarBuilder::new("Processing sequences")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg} | Workers: {prefix}")
                .with_tick()
                .build()?;

            worker_progress.enable_steady_tick(std::time::Duration::from_millis(100));

            // Track active workers
            let mut active_workers = handles.len();
            worker_progress.set_prefix(format!("{}/{} active", active_workers, num_threads));

            // Collect results with progress monitoring
            let mut results = Vec::with_capacity(handles.len());
            for (idx, handle) in handles.into_iter().enumerate() {
                match handle.join() {
                    Ok((processor, local_stats)) => {
                        active_workers -= 1;
                        worker_progress
                            .set_prefix(format!("{}/{} active", active_workers, num_threads));
                        worker_progress.set_message(format!(
                            "Worker {} finished: {} sequences processed",
                            idx + 1,
                            local_stats.processed
                        ));

                        stats.processed += local_stats.processed;
                        stats.errors += local_stats.errors;
                        results.push(processor);
                    }
                    Err(e) => {
                        eprintln!("Worker thread {} panicked: {:?}", idx + 1, e);
                        stats.errors += 1;
                    }
                }
            }

            worker_progress.finish_with_message(format!(
                "All workers complete. Total processed: {} sequences",
                stats.processed
            ));

            // Create merge progress bar
            let merge_progress = ProgressBarBuilder::new("Merging results")
                .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
                .with_tick()
                .build()?;

            // Merge results with progress updates
            Self::merge_worker_results(&results, processor, &merge_progress)?;

            merge_progress.finish_with_message(format!(
                "Processing complete. Total sequences: {}",
                stats.processed
            ));

            Ok(stats)
        }
    }
}
