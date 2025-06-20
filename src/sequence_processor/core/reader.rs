use super::processor::SequenceProcessor;
use super::stats::ProcessingStats;
use anyhow::Result;
use indicatif::ProgressBar;

pub trait SequenceReader {
    fn read_sequences_single_thread<P: SequenceProcessor>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
    ) -> Result<ProcessingStats>;

    fn read_sequences_with_threads<P: SequenceProcessor + Clone + 'static>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
        num_threads: usize,
    ) -> Result<ProcessingStats>;
}