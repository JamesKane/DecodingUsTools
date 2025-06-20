use super::types::*;
use crate::sequence_processor::core::{Sequence, SequenceProcessor};
use crate::sequence_processor::collectors::base::StatsCollector;
use crate::sequence_processor::ProcessingStats;
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::sync::Arc;
use indicatif::ProgressBar;

#[derive(Default, Debug, Serialize, Deserialize, Clone)]
pub struct GraphCoverageProcessor {
    stats: GraphCoverageStats,
    min_mapping_quality: u8,
    min_coverage: u32,
    #[serde(skip)]
    progress: Option<Arc<ProgressBar>>,
}

impl GraphCoverageProcessor {
    pub fn new(min_mapping_quality: u8, min_coverage: u32) -> Self {
        Self {
            stats: GraphCoverageStats::default(),
            min_mapping_quality,
            min_coverage,
            progress: None,
        }
    }

    pub fn with_progress(mut self, progress: ProgressBar) -> Self {
        self.progress = Some(Arc::new(progress));
        self
    }

    pub fn get_stats(&self) -> &GraphCoverageStats {
        &self.stats
    }
}

impl StatsCollector for GraphCoverageProcessor {
    fn process_base(&mut self, _position: u64, _base: u8, _quality: u8, mapping_quality: u8) -> Result<()> {
        if mapping_quality >= self.min_mapping_quality {
            self.stats.mapped_reads += 1;
        } else {
            self.stats.low_quality_reads += 1;
        }
        Ok(())
    }

    fn skip_base(&mut self, _position: u64) -> Result<()> {
        self.stats.unmapped_reads += 1;
        Ok(())
    }

    fn get_min_base_quality(&self) -> u8 {
        0 // We don't filter by base quality
    }

    fn get_min_mapping_quality(&self) -> u8 {
        self.min_mapping_quality
    }

    fn merge_with(&mut self, other: &Self) -> Result<()> {
        // Merge node coverage
        for (node_id, other_coverage) in &other.stats.nodes {
            let coverage = self.stats.nodes
                .entry(node_id.clone())
                .or_default();
            coverage.coverage += other_coverage.coverage;
            coverage.forward_coverage += other_coverage.forward_coverage;
            coverage.reverse_coverage += other_coverage.reverse_coverage;
            coverage.total_quality += other_coverage.total_quality;
            coverage.unique_reads.extend(other_coverage.unique_reads.clone());
        }

        // Merge edges
        self.stats.edges.extend(other.stats.edges.clone());

        // Merge counters
        self.stats.total_reads += other.stats.total_reads;
        self.stats.mapped_reads += other.stats.mapped_reads;
        self.stats.unmapped_reads += other.stats.unmapped_reads;
        self.stats.low_quality_reads += other.stats.low_quality_reads;

        Ok(())
    }
}

impl SequenceProcessor for GraphCoverageProcessor {
    fn process_sequence(&mut self, sequence: &Sequence) -> Result<()> {
        self.stats.total_reads += 1;
        if let Some(progress) = &self.progress {
            progress.inc(1);
        }

        if !sequence.metadata.is_mapped {
            self.stats.unmapped_reads += 1;
            return Ok(());
        }

        if let Some(mapping_quality) = sequence.metadata.mapping_quality {
            if mapping_quality < self.min_mapping_quality {
                self.stats.low_quality_reads += 1;
                return Ok(());
            }
        }

        // Process graph path information if available
        if let Some(chromosome) = &sequence.metadata.chromosome {
            if !chromosome.is_empty() {
                let node_coverage = self.stats.nodes
                    .entry(chromosome.clone())
                    .or_default();

                node_coverage.coverage += 1;
                if let Some(mapping_quality) = sequence.metadata.mapping_quality {
                    node_coverage.total_quality += mapping_quality as f64;
                }
                if let Some(id) = &sequence.id {
                    node_coverage.unique_reads.insert(id.clone());
                }
            }
        }

        Ok(())
    }

    fn get_min_length(&self) -> usize {
        0 // Process all sequences
    }

    fn update_progress(&mut self, stats: &ProcessingStats) {
        if let Some(progress) = &self.progress {
            progress.set_position(stats.processed);
        }
    }

    fn supports_parallel(&self) -> bool {
        true
    }

    fn merge_processor(&mut self, other: &Self) -> Result<()> {
        self.merge_with(other)
    }
}