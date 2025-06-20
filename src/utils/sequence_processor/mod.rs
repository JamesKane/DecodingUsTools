pub mod core;
pub mod threading;
pub mod readers;

// Re-export commonly used items
pub use core::{Sequence, SequenceProcessor, SequenceReader, ProcessingStats, SequenceMetadata};
pub use readers::{BamReader, FastqReader, GamReader};