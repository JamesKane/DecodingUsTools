use crate::sequence_processor::collectors::base::StatsCollector;
use crate::sequence_processor::core::{ProcessingStats, Sequence, SequenceProcessor};
use indicatif::ProgressBar;
use std::sync::Arc;
use anyhow::Result;

/// A processor that handles sequences base-by-base using a StatsCollector
#[derive(Clone)]
pub struct BaseProcessor<C: StatsCollector> {
    collector: C,
    progress: Arc<ProgressBar>,
}

impl<C: StatsCollector> BaseProcessor<C> {
    /// Create a new BaseProcessor with the given collector and progress bar
    pub fn new(collector: C, progress: ProgressBar) -> Self {
        Self {
            collector,
            progress: Arc::new(progress),
        }
    }

    /// Get a reference to the underlying collector
    pub fn collector(&self) -> &C {
        &self.collector
    }

    /// Get a mutable reference to the underlying collector
    pub fn collector_mut(&mut self) -> &mut C {
        &mut self.collector
    }
}

impl<C: StatsCollector> SequenceProcessor for BaseProcessor<C> {
    fn process_sequence(&mut self, sequence: &Sequence) -> Result<()> {
        let position = sequence.metadata.position.unwrap_or(0);
        let mapping_quality = sequence.metadata.mapping_quality.unwrap_or(0);

        for (i, &base) in sequence.data.iter().enumerate() {
            if let Some(quality) = sequence.quality.as_ref().map(|q| q[i]) {
                if quality >= self.collector.get_min_base_quality()
                    && mapping_quality >= self.collector.get_min_mapping_quality() {
                    self.collector.process_base(
                        position + i as u64,
                        base,
                        quality,
                        mapping_quality,
                    )?;
                } else {
                    self.collector.skip_base(position + i as u64)?;
                }
            }
        }

        Ok(())
    }

    fn get_min_length(&self) -> usize {
        1 // Process all sequences
    }

    fn update_progress(&mut self, stats: &ProcessingStats) {
        self.progress.set_position(stats.processed);
    }

    fn supports_parallel(&self) -> bool {
        true
    }

    fn merge_processor(&mut self, other: &Self) -> Result<()> {
        // Delegate merging to the collector
        self.collector.merge_with(&other.collector)
    }
}