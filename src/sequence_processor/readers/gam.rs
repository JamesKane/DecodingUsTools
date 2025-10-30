use crate::sequence_processor::core::*;
use crate::sequence_processor::threading::{merge_processors, ThreadPool};
use crate::vg::framing::GroupIterator;
use anyhow::Result;
use indicatif::ProgressBar;
use niffler::get_reader;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

pub struct GamReader {
    reader: BufReader<Box<dyn std::io::Read>>,
}

impl GamReader {
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let (inner_reader, _compression) = get_reader(Box::new(file))?;
        let reader = BufReader::with_capacity(16 * 1024 * 1024, inner_reader);

        Ok(Self { reader })
    }
}

impl SequenceReader for GamReader {
    fn read_sequences_single_thread<P: SequenceProcessor>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
    ) -> Result<ProcessingStats> {
        let mut stats = ProcessingStats::default();
        let mut group_iter = GroupIterator::new(&mut self.reader);
        let mut last_error_was_buffer = false;

        while let Some(group_result) = group_iter.next() {
            let group = match group_result {
                Ok(g) => {
                    last_error_was_buffer = false;
                    g
                }
                Err(e) => {
                    // If we get a buffer fill error and we've already seen one,
                    // assume we've hit EOF and break
                    if e.to_string().contains("failed to fill whole buffer") {
                        if last_error_was_buffer {
                            break;
                        }
                        last_error_was_buffer = true;
                    } else {
                        eprintln!("Warning: Error reading GAM group: {}", e);
                        stats.errors += 1;
                        last_error_was_buffer = false;
                    }
                    continue;
                }
            };

            // Skip non-GAM groups
            if group.type_tag.as_deref() != Some("GAM") {
                continue;
            }

            for msg_bytes in group.messages {
                match protobuf::Message::parse_from_bytes(&msg_bytes) {
                    Ok(alignment) => {
                        let alignment: crate::vg::gam::Alignment = alignment;
                        let seq_data = alignment.sequence.as_bytes();

                        if seq_data.len() >= processor.get_min_length() {
                            let sequence = Sequence {
                                data: seq_data.to_vec(),
                                id: Some(alignment.name.clone()),
                                quality: Some(alignment.quality.clone()),
                                metadata: SequenceMetadata {
                                    is_mapped: alignment.read_mapped,
                                    chromosome: Option::from(alignment.path.as_ref()
                                        .map(|p| p.name.clone())
                                        .unwrap_or_default()),
                                    position: Some(alignment.query_position as u64),
                                    mapping_quality: Some(alignment.mapping_quality as u8),
                                },
                            };

                            if let Err(e) = processor.process_sequence(&sequence) {
                                eprintln!("Error processing sequence: {}", e);
                                stats.errors += 1;
                            } else {
                                stats.processed += 1;
                            }
                        } else {
                            stats.too_short += 1;
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: Error parsing GAM record: {}", e);
                        stats.errors += 1;
                    }
                }

                if stats.processed % 1000 == 0 {
                    processor.update_progress(&stats);
                }
            }
        }

        // Ensure we update the progress one final time at the end
        processor.update_progress(&stats);

        Ok(stats)
    }

    fn read_sequences_with_threads<P: SequenceProcessor + Clone + 'static>(
        &mut self,
        processor: &mut P,
        progress: &ProgressBar,
        num_threads: usize,
    ) -> Result<ProcessingStats> {
        if num_threads <= 1 || !processor.supports_parallel() {
            return self.read_sequences_single_thread(processor, progress);
        }

        let pool = ThreadPool::new(processor.clone(), num_threads)?;
        let mut stats = ProcessingStats::default();
        let mut group_iter = GroupIterator::new(&mut self.reader);
        let mut last_error_was_buffer = false;

        while let Some(group_result) = group_iter.next() {
            let group = match group_result {
                Ok(g) => {
                    last_error_was_buffer = false;
                    g
                }
                Err(e) => {
                    if e.to_string().contains("failed to fill whole buffer") {
                        if last_error_was_buffer {
                            break;
                        }
                        last_error_was_buffer = true;
                    } else {
                        eprintln!("Warning: Error reading GAM group: {}", e);
                        stats.errors += 1;
                        last_error_was_buffer = false;
                    }
                    continue;
                }
            };

            // Skip non-GAM groups
            if group.type_tag.as_deref() != Some("GAM") {
                continue;
            }

            for msg_bytes in group.messages {
                match protobuf::Message::parse_from_bytes(&msg_bytes) {
                    Ok(alignment) => {
                        let alignment: crate::vg::gam::Alignment = alignment;
                        let seq_data = alignment.sequence.as_bytes();

                        if seq_data.len() >= processor.get_min_length() {
                            let sequence = Sequence {
                                data: seq_data.to_vec(),
                                id: Some(alignment.name.clone()),
                                quality: Some(alignment.quality.clone()),
                                metadata: SequenceMetadata {
                                    is_mapped: alignment.read_mapped,
                                    chromosome: Option::from(alignment.path.as_ref()
                                        .map(|p| p.name.clone())
                                        .unwrap_or_default()),
                                    position: Some(alignment.query_position as u64),
                                    mapping_quality: Some(alignment.mapping_quality as u8),
                                },
                            };

                            pool.send(sequence)?;
                        } else {
                            stats.too_short += 1;
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: Error parsing GAM record: {}", e);
                        stats.errors += 1;
                    }
                }
            }
        }

        let (thread_stats, processors) = pool.finish()?;
        stats.processed += thread_stats.processed;
        stats.errors += thread_stats.errors;

        merge_processors(processors, processor)?;

        Ok(stats)
    }
}
