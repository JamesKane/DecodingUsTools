use chrono::{Datelike, Local};
use directories::ProjectDirs;
use indicatif::{ProgressBar, ProgressStyle};
use serde::de::DeserializeOwned;
use std::fs::{self, File};
use std::io::Read;
use std::path::PathBuf;

#[derive(Clone, Copy)]
pub enum TreeType {
    YDNA,
    MTDNA,
}

impl TreeType {
    fn url(&self) -> &'static str {
        match self {
            TreeType::YDNA => "https://www.familytreedna.com/public/y-dna-haplotree/get",
            TreeType::MTDNA => "https://www.familytreedna.com/public/mt-dna-haplotree/get",
        }
    }

    fn cache_prefix(&self) -> &'static str {
        match self {
            TreeType::YDNA => "ytree",
            TreeType::MTDNA => "mttree",
        }
    }
}

pub struct TreeCache {
    cache_dir: PathBuf,
    tree_type: TreeType,
}

impl TreeCache {
    pub fn new(tree_type: TreeType) -> Result<Self, Box<dyn std::error::Error>> {
        let proj_dirs = ProjectDirs::from("com", "decodingus", "decodingus-tools")
            .ok_or("Failed to determine project directories")?;

        let cache_dir = proj_dirs.cache_dir().join(tree_type.cache_prefix());
        fs::create_dir_all(&cache_dir)?;

        Ok(TreeCache {
            cache_dir,
            tree_type,
        })
    }

    fn get_cache_path(&self) -> PathBuf {
        let now = Local::now();
        let year = now.year();
        let week = now.iso_week().week();
        self.cache_dir.join(format!(
            "{}_{year}_w{week:02}.json",
            self.tree_type.cache_prefix()
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

    pub fn get_tree<T>(&self) -> Result<T, Box<dyn std::error::Error>>
    where
        T: DeserializeOwned,
    {
        let cache_path = self.get_cache_path();

        if self.is_cache_valid(&cache_path) {
            let mut file = File::open(&cache_path)?;
            let mut contents = String::new();
            file.read_to_string(&mut contents)?;

            match serde_json::from_str(&contents) {
                Ok(tree) => return Ok(tree),
                Err(_) => {}
            }
        }

        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        progress.set_message(format!(
            "Downloading FTDNA {} tree...",
            match self.tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        ));

        let tree_json = reqwest::blocking::get(self.tree_type.url())?.text()?;
        fs::write(&cache_path, &tree_json)?;

        progress.finish_with_message(format!(
            "{} tree downloaded and cached",
            match self.tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        ));

        Ok(serde_json::from_str(&tree_json)?)
    }
}
