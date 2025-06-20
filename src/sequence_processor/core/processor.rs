use super::stats::ProcessingStats;
use super::sequence::Sequence;
use anyhow::Result;

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