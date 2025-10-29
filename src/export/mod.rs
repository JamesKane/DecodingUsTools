pub mod json;
pub mod pds;
pub mod formats;

use crate::export::formats::coverage::CoverageExport;
use crate::export::formats::fingerprint::FingerprintExport;
use crate::export::formats::haplogroup::HaplogroupExport;
use chrono::{DateTime, Utc};
use serde::de::{Deserialize as DeserializeTrait, Deserializer, Error};
use serde::ser::{Serialize as SerializeTrait, Serializer};
use serde::{Deserialize, Serialize};

/// Root structure for all exports - AT Protocol compatible
#[derive(Debug, Serialize, Deserialize)]
pub struct AnalysisExport {
    #[serde(rename = "$type")]
    pub record_type: String,  // e.g., "app.bsky.genomics.analysis"

    pub version: String,
    #[serde(serialize_with = "serialize_datetime", deserialize_with = "deserialize_datetime")]
    pub created_at: DateTime<Utc>,
    pub tool_version: String,

    #[serde(flatten)]
    pub data: AnalysisData,

    pub metadata: ExportMetadata,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AnalysisData {
    Coverage(CoverageExport),
    Haplogroup(HaplogroupExport),
    Fingerprint(FingerprintExport),
}

fn serialize_datetime<S>(date: &DateTime<Utc>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&date.to_rfc3339())
}

fn deserialize_datetime<'de, D>(deserializer: D) -> Result<DateTime<Utc>, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    DateTime::parse_from_rfc3339(&s)
        .map(|dt| dt.with_timezone(&Utc))
        .map_err(D::Error::custom)
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ExportMetadata {
    pub sample_id: Option<String>,
    pub lab_id: Option<String>,
    pub sequencing_platform: Option<String>,
    pub reference_genome: String,
    pub tags: Vec<String>,
}