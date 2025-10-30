use anyhow::{Context, Result};
use std::process::Command;

pub fn check_samtools() -> Result<()> {
    Command::new("samtools")
        .arg("--version")
        .output()
        .context("samtools not found. Please install samtools (http://www.htslib.org/) and ensure it's in your PATH")
        .map(|_| ())
}