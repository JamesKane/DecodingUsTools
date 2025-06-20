#[derive(Debug, Default)]
pub struct ProcessingStats {
    pub processed: u64,
    pub skipped: u64,
    pub errors: u64,
    pub too_short: u64,
}