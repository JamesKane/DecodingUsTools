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
        //pb.set_style(/* your style */);
        pb.set_message(format!("Processing {}", name));
        pb
    }

    pub fn add_spinner(&self, message: &str) -> ProgressBar {
        let pb = self.multi.add(ProgressBar::new_spinner());
        //pb.set_style(/* your spinner style */);
        pb.set_message(message.to_string());
        pb
    }
}