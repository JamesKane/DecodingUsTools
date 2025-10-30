pub mod client;
pub mod auth;
pub mod records;

use crate::export::AnalysisExport;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone)]
pub struct PdsConfig {
    pub endpoint: String,
    pub identifier: String,  // DID or handle
    pub password: Option<String>,
    pub app_password: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AtProtoRecord {
    pub repo: String,      // DID
    pub collection: String, // e.g., "app.bsky.genomics.analysis"
    pub rkey: Option<String>,
    pub record: AnalysisExport,
}

pub struct PdsClient {
    config: PdsConfig,
    http_client: reqwest::Client,
    access_token: Option<String>,
}

impl PdsClient {
    pub fn new(config: PdsConfig) -> Self {
        Self {
            config,
            http_client: reqwest::Client::new(),
            access_token: None,
        }
    }

    pub async fn authenticate(&mut self) -> Result<(), PdsError> {
        // Implement AT Protocol authentication
        // POST to /xrpc/com.atproto.server.createSession
        todo!()
    }

    pub async fn create_record(&self, record: AtProtoRecord) -> Result<String, PdsError> {
        // POST to /xrpc/com.atproto.repo.createRecord
        todo!()
    }

    pub async fn upload_blob(&self, data: Vec<u8>, mime_type: &str) -> Result<BlobRef, PdsError> {
        // POST to /xrpc/com.atproto.repo.uploadBlob
        // Useful for uploading coverage plots, etc.
        todo!()
    }
}

#[derive(Debug, thiserror::Error)]
pub enum PdsError {
    #[error("Authentication failed: {0}")]
    AuthError(String),
    #[error("Network error: {0}")]
    NetworkError(#[from] reqwest::Error),
    #[error("Invalid record: {0}")]
    InvalidRecord(String),
}

#[derive(Debug, Serialize, Deserialize)]
pub struct BlobRef {
    #[serde(rename = "$type")]
    pub ref_type: String,
    pub r#ref: String,
    pub mime_type: String,
    pub size: usize,
}