use crate::config::config::Config;
use crate::haplogroup::types::HaplogroupTree;
pub(crate) use crate::vendor::TreeProvider;
use chrono::{Datelike, Local};
use directories::ProjectDirs;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::{self, File};
use std::io::Read;
use std::path::PathBuf;


#[derive(Clone, Copy, PartialEq, Eq)]
pub enum TreeType {
    YDNA,
    MTDNA,
}

pub struct TreeCache {
    cache_dir: PathBuf,
    tree_type: TreeType,
    pub(crate) provider: Box<dyn TreeProvider>,
    config: Config,
}


impl TreeCache {
    pub fn new(tree_type: TreeType, provider: crate::cli::TreeProvider) -> Result<Self, Box<dyn std::error::Error>> {
        let proj_dirs = ProjectDirs::from("com", "decodingus", "decodingus-tools")
            .ok_or("Failed to determine project directories")?;

        let provider = crate::vendor::get_provider(provider);
        let cache_dir = proj_dirs.cache_dir().join(provider.cache_prefix(tree_type));
        fs::create_dir_all(&cache_dir)?;
        let config = Config::load();

        Ok(TreeCache {
            cache_dir,
            tree_type,
            provider,
            config,
        })
    }


    fn get_cache_path(&self) -> PathBuf {
        let now = Local::now();
        let year = now.year();
        let week = now.iso_week().week();
        self.cache_dir.join(format!(
            "{}_{year}_w{week:02}.json",
            self.provider.cache_prefix(self.tree_type)
        ))
    }

    fn is_cache_valid(&self, path: &PathBuf) -> bool {
        if !path.exists() {
            return false;
        }

        match fs::metadata(path) {
            Ok(metadata) => {
                let modified = metadata.modified().ok();

                if let Some(modified) = modified {
                    if let Ok(modified) = modified.elapsed() {
                        return modified.as_secs() < 7 * 24 * 60 * 60;
                    }
                }
            }
            Err(_) => return false,
        }
        false
    }

    pub fn get_tree(&self) -> Result<HaplogroupTree, Box<dyn std::error::Error>> {
        let cache_path = self.get_cache_path();

        if self.is_cache_valid(&cache_path) {
            println!("Using cached tree from {}", cache_path.display());
            let mut file = File::open(&cache_path)?;
            let mut contents = String::new();
            file.read_to_string(&mut contents)?;

            if let Ok(tree) = self.provider.parse_tree(&contents) {
                return Ok(tree);
            }
            println!("Cached tree invalid or outdated");
        }

        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        progress.set_message(self.provider.progress_message(self.tree_type));

        println!("Downloading tree from {}", self.provider.url(self.tree_type));
        let client = reqwest::blocking::Client::builder()
            .timeout(std::time::Duration::from_secs(self.config.download_timeout))
            .build()
            .map_err(|e| format!("Failed to build HTTP client: {}", e))?;

        let response = client.get(self.provider.url(self.tree_type))
            .send()
            .map_err(|e| format!("Failed to download tree: {}", e))?;

        let tree_json = response.text()
            .map_err(|e| format!("Failed to read response body: {}", e))?;

        println!("Caching tree to {}", cache_path.display());
        fs::write(&cache_path, &tree_json)?;

        progress.finish_with_message(format!(
            "{} tree downloaded and cached",
            match self.tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        ));

        self.provider.parse_tree(&tree_json)
            .map_err(|e| format!("Failed to parse tree: {}", e).into())
    }
}