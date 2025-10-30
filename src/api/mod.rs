pub mod coverage;

use serde::Serialize;
use std::sync::Arc;
use thiserror::Error;

/// Progress callback for GUI updates
pub type ProgressCallback = Arc<dyn Fn(ProgressEvent) + Send + Sync>;

/// Progress events that can be sent to GUI or CLI
#[derive(Clone, Debug)]
pub enum ProgressEvent {
    Started { task: String },
    Progress { task: String, current: u64, total: u64 },
    Message { task: String, message: String },
    Completed { task: String },
    Error { task: String, error: String },
}

/// Result types that can be serialized to JSON
pub type ApiResult<T> = Result<T, ApiError>;

/// API-level errors
#[derive(Debug, Error, Serialize)]
pub enum ApiError {
    #[error("IO error: {0}")]
    Io(String),

    #[error("Analysis error: {0}")]
    Analysis(String),

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("BAM/CRAM error: {0}")]
    BamError(String),

    #[error("Reference genome error: {0}")]
    ReferenceError(String),
}

// Implement conversions from common error types
impl From<String> for ApiError {
    fn from(s: String) -> Self {
        ApiError::Analysis(s)
    }
}

impl From<std::io::Error> for ApiError {
    fn from(err: std::io::Error) -> Self {
        ApiError::Io(err.to_string())
    }
}

impl From<Box<dyn std::error::Error>> for ApiError {
    fn from(err: Box<dyn std::error::Error>) -> Self {
        ApiError::Analysis(err.to_string())
    }
}

impl From<anyhow::Error> for ApiError {
    fn from(err: anyhow::Error) -> Self {
        ApiError::Analysis(err.to_string())
    }
}