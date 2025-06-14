use directories::ProjectDirs;
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    #[serde(default = "default_download_timeout")]
    pub download_timeout: u64,
}

fn default_download_timeout() -> u64 {
    300
}

impl Default for Config {
    fn default() -> Self {
        Self {
            download_timeout: default_download_timeout(),
        }
    }
}

impl Config {
    pub fn load() -> Self {
        if let Some(proj_dirs) = ProjectDirs::from("com", "decodingus", "decodingus-tools") {
            let config_dir = proj_dirs.config_dir();
            let config_path = config_dir.join("config.toml");

            if config_path.exists() {
                if let Ok(content) = fs::read_to_string(config_path) {
                    if let Ok(config) = toml::from_str(&content) {
                        return config;
                    }
                }
            }
        }
        Config::default()
    }

    pub fn save(&self) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(proj_dirs) = ProjectDirs::from("com", "decodingus", "decodingus-tools") {
            let config_dir = proj_dirs.config_dir();
            fs::create_dir_all(config_dir)?;

            let config_path = config_dir.join("config.toml");
            let content = toml::to_string_pretty(self)?;
            fs::write(config_path, content)?;
        }
        Ok(())
    }
}