#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledState {
    REF_N,
    CALLABLE,
    NO_COVERAGE,
    LOW_COVERAGE,
    EXCESSIVE_COVERAGE,
    POOR_MAPPING_QUALITY,
}

#[derive(Debug)]
pub(crate) struct CoverageRange {
    pub(crate) start: u32,
    pub(crate) end: u32,
    pub(crate) state: CalledState,
}

impl CoverageRange {
    pub(crate) fn new(pos: u32, state: CalledState) -> Self {
        Self {
            start: pos,
            end: pos,
            state,
        }
    }

    pub(crate) fn can_merge(&self, pos: u32, state: CalledState) -> bool {
        self.end + 1 == pos && self.state == state
    }

    pub(crate) fn extend(&mut self, pos: u32) {
        self.end = pos;
    }
}