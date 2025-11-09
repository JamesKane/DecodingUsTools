use directories::ProjectDirs;
use indicatif::{ProgressBar, ProgressStyle};
use niffler;
use std::collections::HashMap;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::types::ReferenceGenome;

// External liftover machinery
use chainfile as chain;
use chain::liftover;
use omics::coordinate::interval::interbase::Interval;

#[derive(Debug)]
pub struct Liftover {
    // Path to cached chain file; None means identity mapping
    chain_path: Option<PathBuf>,
    // Cached liftover machine built once per instance; None for identity mapping
    machine: Option<Arc<liftover::machine::Machine>>,
    // Retain legacy fields for compatibility with existing code paths/logging; kept empty
    chains_by_query: HashMap<String, Vec<()>>,
    intervals_by_query: HashMap<String, Vec<()>>,
}

impl Liftover {
    pub fn load_or_fetch(source: ReferenceGenome, target: ReferenceGenome) -> Result<Self, Box<dyn std::error::Error>> {
        if source == target {
            return Ok(Liftover { chain_path: None, machine: None, chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() });
        }

        let verbose = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().filter(|v| !v.is_empty() && v != "0").is_some();
        let refresh = std::env::var("DECODINGUS_LIFTOVER_REFRESH").ok().filter(|v| !v.is_empty() && v != "0").is_some();

        // Resolve chain file path (overrides first)
        let chosen_path: PathBuf = if let Ok(path) = std::env::var("DECODINGUS_LIFTOVER_FILE") {
            let p = PathBuf::from(&path);
            if verbose { eprintln!("[liftover] override: using local chain file {}", p.display()); }
            p
        } else if let Ok(url) = std::env::var("DECODINGUS_LIFTOVER_URL") {
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
            path
        } else {
            let candidates = chain_url_candidates(&source, &target);
            if candidates.is_empty() {
                return Err(format!("No chain mapping known from {} to {}", source.name(), target.name()).into());
            }
            let mut chosen: Option<PathBuf> = None;
            let mut last_err: Option<Box<dyn std::error::Error>> = None;
            for (url, cache_tag, filename) in candidates {
                const CACHE_VERSION_SALT: &str = "v2";
                let cache_name = format!("{}__{}__{}", CACHE_VERSION_SALT, cache_tag, filename);
                let path = get_cache_path(cache_name)?;
                if refresh { let _ = fs::remove_file(&path); }
                if !is_cache_valid(&path) {
                    if verbose { eprintln!("[liftover] downloading chain: {} -> {:?}", url, path); }
                    if let Err(e) = download_chain(&url, &path) { last_err = Some(e); continue; }
                } else if verbose {
                    eprintln!("[liftover] using cached chain: {:?}", path);
                }
                chosen = Some(path);
                break;
            }
            if let Some(p) = chosen { p } else { return Err(last_err.unwrap_or_else(|| "Failed to load any liftover chain candidate".into())); }
        };

        // Build liftover machine once
        if verbose { eprintln!("[liftover] building liftover machine from {}", chosen_path.display()); }
        let file = fs::File::open(&chosen_path)?;
        let (reader, _comp) = niffler::get_reader(Box::new(file))?;
        let buf = BufReader::new(reader);
        let rdr = chain::Reader::new(buf);
        let machine = liftover::machine::builder::Builder.try_build_from(rdr)?;
        Ok(Liftover {
            chain_path: Some(chosen_path),
            machine: Some(Arc::new(machine)),
            chains_by_query: HashMap::new(),
            intervals_by_query: HashMap::new(),
        })
    }

    /// Special identity liftover used for mitochondrial rCRS ↔ build-specific chrM when
    /// coordinate systems are effectively identical. Keeping a dedicated constructor allows
    /// us to switch to a proper mapping later if needed without touching call sites.
    pub fn identity() -> Self {
        Liftover { chain_path: None, machine: None, chains_by_query: HashMap::new(), intervals_by_query: HashMap::new() }
    }

    /// Map a 1-based query position on a chromosome name to a 1-based target position.
    /// Returns None if the position cannot be mapped.
    pub fn map_pos(&self, query_chrom: &str, query_pos_1based: u32) -> Option<u32> {
        // Identity shortcut when no chains loaded
        if self.machine.is_none() {
            return Some(query_pos_1based);
        }
        let machine = self.machine.as_ref().unwrap();
        let verbose = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().filter(|v| !v.is_empty() && v != "0").is_some();

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

    /// Map a position and return both the target position and whether strand is reversed.
    /// Returns (target_pos_1based, is_reverse_strand) or None if mapping fails.
    /// Use this when you need to know if alleles should be reverse complemented.
    pub fn map_pos_with_strand(&self, query_chrom: &str, query_pos_1based: u32) -> Option<(u32, bool)> {
        // Identity shortcut when no chains loaded
        if self.machine.is_none() {
            return Some((query_pos_1based, false));
        }
        let machine = self.machine.as_ref().unwrap();
        let verbose = std::env::var("DECODINGUS_VERBOSE_LIFTOVER").ok().filter(|v| !v.is_empty() && v != "0").is_some();

        let mut try_map = |chrom: &str| -> Option<(u32, bool)> {
            let start0 = (query_pos_1based as u64).saturating_sub(1);
            let end0 = start0 + 1;
            if verbose { 
                eprintln!("[liftover] query base: {}:{} (0-based {}..{})", chrom, query_pos_1based, start0, end0); 
            }
            
            let q = format!("{}:+:{}-{}", chrom, start0, end0);
            let interval: Interval = match q.parse() {
                Ok(iv) => iv,
                Err(_) => {
                    if verbose { 
                        eprintln!("[liftover] failed to parse interval {}", q); 
                    }
                    return None;
                }
            };
            
            let results_opt = machine.liftover(interval);
            if let Some(results) = results_opt {
                if verbose {
                    eprintln!("[liftover] got {} results for {}:{}", results.len(), chrom, query_pos_1based);
                }
                
                for r in results {
                    let s = r.to_string();
                    if verbose { 
                        eprintln!("[liftover]   result: {}", s); 
                    }
                    
                    // Parse format: "chrom:strand:start-end"
                    let parts: Vec<&str> = s.split(':').collect();
                    if parts.len() >= 3 {
                        let strand = parts[1];
                        let is_reverse = strand == "-";
                        
                        if let Some((start_s, _end_s)) = parts[2].split_once('-') {
                            if let Ok(start0_tgt) = start_s.parse::<u64>() {
                                let pos = (start0_tgt as u32).saturating_add(1);
                                if verbose {
                                    eprintln!("[liftover]   mapped to pos {} (strand {})", pos, strand);
                                }
                                return Some((pos, is_reverse));
                            }
                        }
                    }
                }
            } else if verbose {
                eprintln!("[liftover] no results for {}:{}", chrom, query_pos_1based);
            }
            
            None
        };

        // Try the provided chromosome first
        if let Some(result) = try_map(query_chrom) { return Some(result); }
        
        // Fallback aliases
        if query_chrom.starts_with("chr") {
            let bare = &query_chrom[3..];
            if let Some(result) = try_map(bare) { return Some(result); }
        } else {
            let prefixed = format!("chr{}", query_chrom);
            if let Some(result) = try_map(&prefixed) { return Some(result); }
        }

        if verbose { 
            eprintln!("[liftover] no mapping for {}:{}", query_chrom, query_pos_1based); 
        }
        None
    }

    /// Reverse complement a DNA base
    fn reverse_complement_base(base: char) -> char {
        match base.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        }
    }

    /// Reverse complement a DNA sequence (allele)
    pub fn reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| Self::reverse_complement_base(c))
            .collect()
    }

    /// Batch-map many 1-based positions on a chromosome; preserves input order.
    pub fn map_many(&self, chrom: &str, positions_1based: &[u32]) -> Vec<Option<u32>> {
        positions_1based.iter().map(|p| self.map_pos(chrom, *p)).collect()
    }

    /// Batch-map by chromosome: queries is a map chrom -> positions; returns mapped positions in the same order per chrom.
    pub fn map_many_by_chrom(&self, queries: &HashMap<String, Vec<u32>>) -> HashMap<String, Vec<Option<u32>>> {
        let mut out: HashMap<String, Vec<Option<u32>>> = HashMap::new();
        for (chrom, v) in queries.iter() {
            out.insert(chrom.clone(), self.map_many(chrom, v));
        }
        out
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

