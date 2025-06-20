use super::core::*;
use crate::utils::progress_bar_builder::ProgressBarBuilder;
use anyhow::Result;
use crossbeam_channel::{bounded, Sender};
use std::thread;

pub struct ThreadPool<P: SequenceProcessor> {
    handles: Vec<thread::JoinHandle<(P, ProcessingStats)>>,
    tx: Sender<Sequence>,
    processor: P,
    num_threads: usize,
}

impl<P: SequenceProcessor + Clone + 'static> ThreadPool<P> {
    pub fn new(processor: P, num_threads: usize) -> Result<Self> {
        let (tx, rx) = bounded(num_threads * 2);
        let mut handles = Vec::with_capacity(num_threads);

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

        Ok(ThreadPool {
            handles,
            tx,
            processor,
            num_threads,
        })
    }

    pub fn send(&self, sequence: Sequence) -> Result<()> {
        Ok(self.tx.send(sequence)?)
    }

    pub fn finish(self) -> Result<(ProcessingStats, Vec<P>)> {
        drop(self.tx);

        let collect_progress = ProgressBarBuilder::new("Collecting results")
            .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
            .with_tick()
            .build()?;

        let mut stats = ProcessingStats::default();
        let mut processors = Vec::with_capacity(self.handles.len());

        for (idx, handle) in self.handles.into_iter().enumerate() {
            collect_progress.set_message(format!(
                "Collecting worker {} of {}",
                idx + 1,
                self.num_threads
            ));

            let (worker_processor, worker_stats) = handle.join().unwrap();
            processors.push(worker_processor);
            stats.processed += worker_stats.processed;
            stats.errors += worker_stats.errors;
            collect_progress.inc(1);
        }

        Ok((stats, processors))
    }
}

pub fn merge_processors<P: SequenceProcessor>(
    processors: Vec<P>,
    main_processor: &mut P,
) -> Result<()> {
    let merge_progress = ProgressBarBuilder::new("Merging results")
        .with_template("{spinner:.green} [{elapsed_precise}] {msg}")
        .with_tick()
        .build()?;

    for (idx, worker_processor) in processors.iter().enumerate() {
        merge_progress.set_message(format!(
            "Merging processor {} of {}",
            idx + 1,
            processors.len()
        ));
        main_processor.merge_processor(worker_processor)?;
        merge_progress.inc(1);
    }

    Ok(())
}