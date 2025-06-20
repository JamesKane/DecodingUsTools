use anyhow::Result;

/// A trait for processing individual bases in a sequence while collecting statistics
pub trait StatsCollector: Send + Clone + 'static {
    /// Process a single base at the given position
    fn process_base(
        &mut self,
        position: u64,
        base: u8,
        quality: u8,
        mapping_quality: u8,
    ) -> Result<()>;

    /// Skip processing a base at the given position
    fn skip_base(&mut self, position: u64) -> Result<()>;

    /// Get the minimum quality threshold for bases to be processed
    fn get_min_base_quality(&self) -> u8;

    /// Get the minimum mapping quality threshold for reads to be processed
    fn get_min_mapping_quality(&self) -> u8;

    /// Merge statistics from another collector of the same type
    fn merge_with(&mut self, other: &Self) -> Result<()>;
}
