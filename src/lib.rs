pub mod api;
pub mod bam_fixer;
pub mod callable_loci;
pub mod cli;  // Add this
pub mod commands;
pub mod config;
pub mod haplogroup;
pub mod sequence_processor;
pub mod types;
pub mod utils;
pub mod vendor;
pub mod vg;
pub mod generated;
pub mod export;
mod error;

// Re-export main API
pub use api::*;