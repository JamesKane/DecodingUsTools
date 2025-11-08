use directories::ProjectDirs;
use indicatif::{ProgressBar, ProgressStyle};
use niffler;
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use crate::types::ReferenceGenome;

// External liftover machinery
use chainfile as chain;
use chain::liftover;
use omics::coordinate::interval::interbase::Interval;

#[derive(Debug)]
pub struct Liftover {
    // Path to cached chain file; None means identity mapping
    chain_path: Option<PathBuf>,
    // Retain legacy fields for compatibility with existing code paths/logging; kept empty
    chains_by_query: HashMap<String, Vec<()>>,
    intervals_by_query: HashMap<String, Vec<()>>,
}

impl Liftover {
    pub fn load_or_fetch(source: ReferenceGenome, target: ReferenceGenome) -> Result<Self, Box<dyn std::error::Error>> {
        if source == target {
            return Ok(Liftover { chain_path: None, chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() });
        }

        let verbose = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().filter(|v| !v.is_empty() && v != "0").is_some();
        let refresh = std::env::var("DECODINGUS_LIFTOVER_REFRESH").ok().filter(|v| !v.is_empty() && v != "0").is_some();

        // Overrides: explicit local file or URL
        if let Ok(path) = std::env::var("DECODINGUS_LIFTOVER_FILE") {
            let p = PathBuf::from(&path);
            if verbose { eprintln!("[liftover] override: using local chain file {}", p.display()); }
            return Ok(Liftover { chain_path: Some(p), chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() });
        }
        if let Ok(url) = std::env::var("DECODINGUS_LIFTOVER_URL") {
            // Sanitize URL into a cache file name
            let mut tag = url.replace(|c: char| !c.is_ascii_alphanumeric(), "_");
            if tag.len() > 80 { tag.truncate(80); }
            const CACHE_VERSION_SALT: &str = "v2";
            let cache_name = format!("override__{}__{}.chain.gz", CACHE_VERSION_SALT, tag);
            let path = get_cache_path(cache_name)?;
            if refresh { let _ = fs::remove_file(&path); }
            if !is_cache_valid(&path) {
                if verbose { eprintln!("[liftover] override: downloading {} -> {:?}", url, path); }
                download_chain(&url, &path)?;
            } else if verbose { eprintln!("[liftover] override: using cached {:?}", path); }
            return Ok(Liftover { chain_path: Some(path), chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() });
        }

        let candidates = chain_url_candidates(&source, &target);
        if candidates.is_empty() {
            return Err(format!("No chain mapping known from {} to {}", source.name(), target.name()).into());
        }

        let mut last_err: Option<Box<dyn std::error::Error>> = None;
        for (url, cache_tag, filename) in candidates {
            // add cache version salt to avoid stale semantics reuse
            const CACHE_VERSION_SALT: &str = "v2";
            // create a unique cache file name per candidate
            let cache_name = format!("{}__{}__{}", CACHE_VERSION_SALT, cache_tag, filename);
            let path = get_cache_path(cache_name)?;
            if refresh {
                let _ = fs::remove_file(&path);
            }
            if !is_cache_valid(&path) {
                if verbose { eprintln!("[liftover] downloading chain: {} -> {:?}", url, path); }
                if let Err(e) = download_chain(&url, &path) { last_err = Some(e); continue; }
            } else if verbose {
                eprintln!("[liftover] using cached chain: {:?}", path);
            }
            // We consider the first successfully cached file as usable
            return Ok(Liftover { chain_path: Some(path), chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() });
        }
        Err(last_err.unwrap_or_else(|| "Failed to load any liftover chain candidate".into()))
    }

    /// Special identity liftover used for mitochondrial rCRS ↔ build-specific chrM when
    /// coordinate systems are effectively identical. Keeping a dedicated constructor allows
    /// us to switch to a proper mapping later if needed without touching call sites.
    pub fn identity() -> Self {
        Liftover { chain_path: None, chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() }
    }

    /// Map a 1-based query position on a chromosome name to a 1-based target position.
    /// Returns None if the position cannot be mapped.
    pub fn map_pos(&self, query_chrom: &str, query_pos_1based: u32) -> Option<u32> {
        // Identity shortcut when no chains loaded
        if self.chain_path.is_none() {
            return Some(query_pos_1based);
        }
        let chain_path = self.chain_path.as_ref().unwrap();

        let verbose = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().filter(|v| !v.is_empty() && v != "0").is_some();

        // Open chain file (gz or plain) using niffler for transparency ONCE
        if verbose { eprintln!("[liftover] using chain file: {}", chain_path.display()); }
        let file = match fs::File::open(chain_path) {
            Ok(f) => f,
            Err(e) => {
                if verbose { eprintln!("[liftover] failed to open chain file: {}", e); }
                return None;
            }
        };
        let (reader, _comp) = match niffler::get_reader(Box::new(file)) {
            Ok(t) => t,
            Err(e) => {
                if verbose { eprintln!("[liftover] failed to create decompressor: {}", e); }
                return None;
            }
        };
        let buf = BufReader::new(reader);
        let rdr = chain::Reader::new(buf);
        // Build liftover machine
        let machine = match liftover::machine::builder::Builder.try_build_from(rdr) {
            Ok(m) => m,
            Err(e) => {
                if verbose { eprintln!("[liftover] failed to build liftover machine: {}", e); }
                return None;
            }
        };

        // Helper to construct an interbase interval string in the format expected by chainfile/omics: "chr:+:start-end"
        let mut make_interval = |chrom: &str, start0: u64, end0: u64| -> Option<Interval> {
            let q = format!("{}:+:{}-{}", chrom, start0, end0);
            if verbose { eprintln!("[liftover] parse interval: {}", q); }
            q.parse::<Interval>().ok()
        };

        let mut try_map = |chrom: &str| -> Option<u32> {
            // Build an interbase interval [start, start+1) from 1-based coordinate
            let start0 = (query_pos_1based as u64).saturating_sub(1);
            let end0 = start0 + 1;
            if verbose { eprintln!("[liftover] query base: {}:{} (0-based {}..{})", chrom, query_pos_1based, start0, end0); }
            let interval: Interval = match make_interval(chrom, start0, end0) {
                Some(iv) => iv,
                None => {
                    if verbose { eprintln!("[liftover] failed to build interval for {}:{}-{}", chrom, start0, end0); }
                    return None;
                }
            };
            let results_opt = machine.liftover(interval);
            if results_opt.is_none() {
                if verbose { eprintln!("[liftover] no results for {}:{}-{}", chrom, start0, end0); }
            } else {
                let results = results_opt.unwrap();
                if verbose {
                    let mut first = true;
                    for r in &results {
                        if first { eprintln!("[liftover] got results:"); first = false; }
                        eprintln!("  -> {}", r);
                    }
                }

                // Prefer first result; parse start from Display form "chr:+:start-end"
                for r in results {
                    let s = r.query().to_string(); // format: chrom:+:start-end
                    let mut parts = s.split(':');
                    let _chrom = parts.next();
                    let _strand = parts.next();
                    if let Some(coords) = parts.next() {
                        if let Some((start_s, _end_s)) = coords.split_once('-') {
                            if let Ok(start0_tgt) = start_s.parse::<u64>() {
                                return Some((start0_tgt as u32).saturating_add(1));
                            }
                        }
                    }
                }
            }

            // Window search around the query base (helps at block edges, esp. chrY)
            let max_win: i64 = 10;
            let q0 = query_pos_1based as i64 - 1;
            for delta in 1..=max_win {
                for sign in [-1i64, 1i64] {
                    let pos0 = q0 + sign * delta;
                    if pos0 < 0 { continue; }
                    let s0 = pos0 as u64;
                    let e0 = s0 + 1;
                    let q = format!("{}:{}-{}", chrom, s0, e0);
                    if verbose { eprintln!("[liftover] window try delta {}: {}", sign*delta, q); }
                    let iv: Interval = match q.parse() {
                        Ok(iv) => iv,
                        Err(_) => continue,
                    };
                    if let Some(rs) = machine.liftover(iv) {
                        for r in rs {
                            let s = r.to_string();
                            if let Some((_, coords)) = s.split_once(':') {
                                if let Some((start_s, _end_s)) = coords.split_once('-') {
                                    if let Ok(start0_tgt) = start_s.parse::<i64>() {
                                        let adj = start0_tgt + (q0 - pos0);
                                        if adj >= 0 && adj <= i64::from(u32::MAX) {
                                            if verbose { eprintln!("[liftover] window mapped: base {} via {} -> start {} (adj {})", query_pos_1based, q, start0_tgt, adj); }
                                            return Some((adj as u32).saturating_add(1));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            None
        };

        // Try the provided chromosome first
        if let Some(p) = try_map(query_chrom) { return Some(p); }
        // Fallback aliases (chrY <-> Y, chrM <-> MT/M)
        if query_chrom.starts_with("chr") {
            let bare = &query_chrom[3..];
            if bare == "Y" {
                if let Some(p) = try_map("Y") { return Some(p); }
            } else if bare == "M" {
                if let Some(p) = try_map("MT") { return Some(p); }
            }
            // also try bare name directly
            if let Some(p) = try_map(bare) { return Some(p); }
        } else {
            if query_chrom == "Y" {
                if let Some(p) = try_map("chrY") { return Some(p); }
            } else if query_chrom == "MT" || query_chrom == "M" {
                if let Some(p) = try_map("chrM") { return Some(p); }
            } else {
                let prefixed = format!("chr{}", query_chrom);
                if let Some(p) = try_map(&prefixed) { return Some(p); }
            }
        }

        if verbose { eprintln!("[liftover] no mapping for {}:{}", query_chrom, query_pos_1based); }
        None
    }
}

fn chain_url_and_name(_source: &ReferenceGenome, _target: &ReferenceGenome) -> Option<(String, String)> { // legacy, unused
    None
}

fn ucsc_dir_code(genome: &ReferenceGenome) -> &'static str {
    match genome {
        ReferenceGenome::GRCh37 => "hg19",
        ReferenceGenome::GRCh38 => "hg38",
        ReferenceGenome::CHM13v2 => "hs1",
    }
}

fn ucsc_file_label(genome: &ReferenceGenome) -> &'static str {
    match genome {
        ReferenceGenome::GRCh37 => "Hg19",
        ReferenceGenome::GRCh38 => "Hg38",
        ReferenceGenome::CHM13v2 => "Hs1",
    }
}

fn chain_url_candidates(source: &ReferenceGenome, target: &ReferenceGenome) -> Vec<(String, String, String)> {
    // Simplified: use the canonical UCSC naming. Note that UCSC filenames are case-sensitive.
    // Special-cases for known pairs to avoid case/casing issues.
    // GRCh38 -> CHM13v2.0
    if let (ReferenceGenome::GRCh38, ReferenceGenome::CHM13v2) = (source, target) {
        let filename = "hg38ToHs1.over.chain.gz".to_string();
        let url = format!("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/{}", filename);
        return vec![(url, "hg38".to_string(), filename)];
    }
    // CHM13v2.0 -> GRCh38
    if let (ReferenceGenome::CHM13v2, ReferenceGenome::GRCh38) = (source, target) {
        // UCSC uses lowercase "hs1" in the filename here
        let filename = "hs1ToHg38.over.chain.gz".to_string();
        let url = format!("https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/{}", filename);
        return vec![(url, "hs1".to_string(), filename)];
    }

    // Fallback generic construction for other directions
    let src_dir = ucsc_dir_code(source);
    let tgt_label = ucsc_file_label(target);
    let src_label = ucsc_file_label(source);

    let filename = format!("{}To{}.over.chain.gz", src_label, tgt_label);
    let url = format!("https://hgdownload.soe.ucsc.edu/goldenPath/{}/liftOver/{}", src_dir, filename);
    vec![(url, src_dir.to_string(), filename)]
}

fn capitalize_first(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

fn get_cache_path(filename: String) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let proj_dirs = ProjectDirs::from("com", "decodingus", "decodingus-tools")
        .ok_or("Failed to determine project directories")?;
    let cache_dir = proj_dirs.cache_dir().join("liftover");
    fs::create_dir_all(&cache_dir)?;
    Ok(cache_dir.join(filename))
}

fn is_cache_valid(path: &Path) -> bool {
    if !path.exists() { return false; }
    if let Ok(metadata) = fs::metadata(path) {
        if let (Some(modified),) = (metadata.modified().ok(),) {
            if let Ok(elapsed) = modified.elapsed() {
                return elapsed.as_secs() < 14 * 24 * 3600;
            }
        }
    }
    false
}

fn download_chain(url: &str, dest: &Path) -> Result<(), Box<dyn std::error::Error>> {
    let progress = ProgressBar::new_spinner();
    progress.set_style(ProgressStyle::default_spinner().template("{spinner:.green} {msg}").unwrap());
    progress.set_message(format!("Downloading liftover chain: {}", url));

    let resp = reqwest::blocking::get(url)?;
    if !resp.status().is_success() {
        progress.finish_and_clear();
        return Err(format!("Failed to download liftover chain (status {}): {}", resp.status(), url).into());
    }
    let bytes = resp.bytes()?;
    if bytes.len() < 100 {
        // Likely an error page or empty file; avoid caching
        progress.finish_and_clear();
        return Err(format!("Downloaded liftover chain is unexpectedly small ({} bytes): {}", bytes.len(), url).into());
    }
    fs::write(dest, &bytes)?;
    progress.finish_with_message("Liftover chain cached");
    Ok(())
}

