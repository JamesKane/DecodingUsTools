use crate::haplogroup::{Haplogroup, HaplogroupTree};
pub(crate) use crate::vendor::TreeProvider;
use chrono::{Datelike, Local};
use directories::ProjectDirs;
use indicatif::{ProgressBar, ProgressStyle};
use std::error::Error;
use std::fs::{self, File};
use std::io::Read;
use std::path::PathBuf;

#[derive(Clone, Copy)]
pub enum TreeType {
    YDNA,
    MTDNA,
}

struct FtdnaTreeProvider;

impl TreeProvider for FtdnaTreeProvider {
    fn url(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "https://www.familytreedna.com/public/y-dna-haplotree/get",
            TreeType::MTDNA => "https://www.familytreedna.com/public/mt-dna-haplotree/get",
        }
    }

    fn cache_prefix(&self, tree_type: TreeType) -> &str {
        match tree_type {
            TreeType::YDNA => "ytree",
            TreeType::MTDNA => "mttree",
        }
    }

    fn progress_message(&self, tree_type: TreeType) -> String {
        format!(
            "Downloading FTDNA {} tree...",
            match tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        )
    }

    fn parse_tree(&self, data: &str) -> Result<HaplogroupTree, Box<dyn Error>> {
        crate::vendor::ftdna::FtdnaTreeProvider::new().parse_tree(data)
    }

    fn build_tree(&self, tree: &HaplogroupTree, node_id: u32, tree_type: TreeType) -> Option<Haplogroup> {
        crate::vendor::ftdna::FtdnaTreeProvider::new().build_tree(tree, node_id, tree_type)
    }
}

pub struct TreeCache {
    cache_dir: PathBuf,
    tree_type: TreeType,
    pub(crate) provider: Box<dyn TreeProvider>,
}

impl TreeCache {
    pub fn new(tree_type: TreeType) -> Result<Self, Box<dyn std::error::Error>> {
        let proj_dirs = ProjectDirs::from("com", "decodingus", "decodingus-tools")
            .ok_or("Failed to determine project directories")?;

        // Default to FTDNA provider for now
        let provider = Box::new(FtdnaTreeProvider);
        let cache_dir = proj_dirs.cache_dir().join(provider.cache_prefix(tree_type));
        fs::create_dir_all(&cache_dir)?;

        Ok(TreeCache {
            cache_dir,
            tree_type,
            provider,
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
            let mut file = File::open(&cache_path)?;
            let mut contents = String::new();
            file.read_to_string(&mut contents)?;

            if let Ok(tree) = self.provider.parse_tree(&contents) {
                return Ok(tree);
            }
        }

        let progress = ProgressBar::new_spinner();
        progress.set_style(
            ProgressStyle::default_spinner()
                .template("{spinner:.green} {msg}")
                .unwrap(),
        );
        progress.set_message(self.provider.progress_message(self.tree_type));

        let client = reqwest::blocking::Client::builder()
            .timeout(std::time::Duration::from_secs(300))
            .build()?;
        let tree_json = client.get(self.provider.url(self.tree_type)).send()?.text()?;

        fs::write(&cache_path, &tree_json)?;

        progress.finish_with_message(format!(
            "{} tree downloaded and cached",
            match self.tree_type {
                TreeType::YDNA => "Y-DNA",
                TreeType::MTDNA => "MT-DNA",
            }
        ));

        self.provider.parse_tree(&tree_json)
    }
}