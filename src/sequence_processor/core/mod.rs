pub(crate) mod processor;
pub(crate) mod sequence;
pub(crate) mod stats;
mod reader;

pub use processor::SequenceProcessor;
pub use sequence::{Sequence, SequenceMetadata};
pub use stats::ProcessingStats;
pub use reader::SequenceReader;
