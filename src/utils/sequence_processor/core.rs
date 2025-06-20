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

    fn read_sequences_with_threads<P: SequenceProcessor + Clone + 'static>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
        num_threads: usize,
    ) -> Result<ProcessingStats>;
}