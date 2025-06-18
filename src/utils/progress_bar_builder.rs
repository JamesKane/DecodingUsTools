use indicatif::{ProgressBar, ProgressStyle};
use anyhow::Result;
use std::time::Duration;

pub(crate) struct ProgressBarBuilder {
    style_template: &'static str,
    message: String,
    spinner: bool,
    enable_tick: bool,
}

impl ProgressBarBuilder {
    pub(crate) fn new(message: impl Into<String>) -> Self {
        Self {
            style_template: "{spinner:.green} {msg}",
            message: message.into(),
            spinner: true,
            enable_tick: false,
        }
    }

    pub(crate) fn with_template(mut self, template: &'static str) -> Self {
        self.style_template = template;
        self
    }

    pub(crate) fn with_progress_bar(mut self) -> Self {
        self.spinner = false;
        self
    }

    pub(crate) fn with_tick(mut self) -> Self {
        self.enable_tick = true;
        self
    }

    pub(crate) fn build(self) -> Result<ProgressBar> {
        let pb = if self.spinner {
            ProgressBar::new_spinner()
        } else {
            ProgressBar::new(0)
        };

        pb.set_style(ProgressStyle::default_spinner().template(self.style_template)?);
        pb.set_message(self.message);

        if self.enable_tick {
            pb.enable_steady_tick(Duration::from_secs(5));
        }

        Ok(pb)
    }
}