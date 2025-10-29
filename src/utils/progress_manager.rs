use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

pub struct ProgressManager {
    multi: MultiProgress,
}

impl ProgressManager {
    pub fn new() -> Self {
        Self { multi: MultiProgress::new() }
    }

    pub fn add_contig_bar(&self, name: &str, length: usize) -> ProgressBar {
        let pb = self.multi.add(ProgressBar::new(length as u64));
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} {msg}")
                .unwrap()
                .progress_chars("#>-")
        );
        pb.set_message(format!("Processing {}", name));
        pb
    }

    pub fn add_spinner(&self, message: &str) -> ProgressBar {
        let pb = self.multi.add(ProgressBar::new_spinner());
        pb.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap()
        );
        pb.set_message(message.to_string());
        pb
    }
}