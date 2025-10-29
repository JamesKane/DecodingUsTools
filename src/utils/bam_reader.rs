// src/utils/bam_reader.rs
use rust_htslib::bam::{self, Read};
use std::error::Error;

pub struct BamReaderFactory;

impl BamReaderFactory {
    pub fn open_indexed(bam_path: &str, reference_path: Option<&str>) -> Result<bam::IndexedReader, Box<dyn Error>> {
        if let Some(ref_path) = reference_path {
            if bam_path.ends_with(".cram") {
                std::env::set_var("REF_PATH", ref_path);
            }
        }
        Ok(bam::IndexedReader::from_path(bam_path)?)
    }

    pub fn open(bam_path: &str, reference_path: Option<&str>) -> Result<bam::Reader, Box<dyn Error>> {
        if let Some(ref_path) = reference_path {
            if bam_path.ends_with(".cram") {
                std::env::set_var("REF_PATH", ref_path);
            }
        }
        Ok(bam::Reader::from_path(bam_path)?)
    }
}