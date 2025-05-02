use clap::Parser;
use rust_htslib::{bam, bam::Read};

/// Program to analyze BAM file callability
#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the BAM file
    #[arg(required = true)]
    bam_file: String,
}


// Define thresholds for callable regions
const MIN_DEPTH: u32 = 10;  // Minimum depth for callable regions
const MAX_DEPTH: u32 = 250; // Maximum depth before considering excessive coverage
const MIN_MAPPING_QUALITY: u8 = 20; // Minimum mapping quality for callable regions

#[derive(Default)]
struct CallableStats {
    no_coverage: usize,
    low_coverage: usize,
    excessive_coverage: usize,
    poor_mapping_quality: usize,
    callable: usize,
}

#[derive(Default)]
struct ReferenceInfo {
    length: usize,
    coverage: Vec<u32>,
    mapping_quality: Vec<u8>,
    all_mapping_qualities: Vec<u8>,
}

impl ReferenceInfo {
    fn new(length: usize) -> Self {
        ReferenceInfo {
            length,
            coverage: vec![0u32; length],
            mapping_quality: vec![0u8; length],
            all_mapping_qualities: vec![],
        }
    }
}


pub fn main() {
    let args = Args::parse();

    // Open BAM file with error handling
    let mut bam = match bam::Reader::from_path(&args.bam_file) {
        Ok(reader) => reader,
        Err(e) => {
            eprintln!("Error opening BAM file '{}': {}", args.bam_file, e);
            std::process::exit(1);
        }
    };

    let header = bam::Header::from_template(bam.header());
    let bam_header = bam.header().clone();
    
    let mut ref_info: std::collections::HashMap<String, ReferenceInfo> = std::collections::HashMap::new();

    // Collect reference sequence names and lengths
    for (key, records) in header.to_hashmap() {
        if key == "SQ" {
            for record in records {
                if let (Some(sn), Some(ln)) = (record.get("SN"), record.get("LN")) {
                    if let Ok(length) = ln.parse::<usize>() {
                        ref_info.insert(sn.to_string(), ReferenceInfo::new(length));
                    }
                }
            }
        }
    }


    // Set up pileup
    let mut pileup = bam.pileup();
    pileup.set_max_depth(1000000);

    // Process pileups to calculate depth and mapping quality
    for p in pileup {
        let pileup = p.unwrap();
        let tid: i32 = pileup.tid().try_into().unwrap();
        let tid_u32: u32 = tid.try_into().unwrap();
        let ref_name = String::from_utf8(bam_header.tid2name(tid_u32).to_owned()).unwrap();
        let pos = pileup.pos() as usize;

        if let Some(ref_info) = ref_info.get_mut(&ref_name) {
            if pos < ref_info.length {
                let depth = pileup.depth();
                ref_info.coverage[pos] = depth;

                // Collect all non-zero mapping qualities
                ref_info.all_mapping_qualities.extend(
                    pileup
                        .alignments()
                        .map(|aln| aln.record().mapq())
                        .filter(|&q| q > 0)
                );

                // Calculate average mapping quality for this position
                let mapq_sum: u32 = pileup
                    .alignments()
                    .map(|aln| aln.record().mapq() as u32)
                    .sum();

                let avg_mapq = if depth > 0 {
                    (mapq_sum / depth) as u8
                } else {
                    0
                };

                ref_info.mapping_quality[pos] = avg_mapq;
            }
        }
    }


    // Calculate and print coverage and callable statistics
    for (ref_name, ref_info) in ref_info {
        let mut stats = CallableStats::default();

        // Analyze each position
        for i in 0..ref_info.length {
            let depth = ref_info.coverage[i];
            let mapping_quality = ref_info.mapping_quality[i];

            match depth {
                0 => stats.no_coverage += 1,
                1..=MIN_DEPTH => stats.low_coverage += 1,
                d if d > MAX_DEPTH => stats.excessive_coverage += 1,
                _ => {
                    if mapping_quality < MIN_MAPPING_QUALITY {
                        stats.poor_mapping_quality += 1;
                    } else {
                        stats.callable += 1;
                    }
                }
            }
        }

        // Calculate median mapping quality
        if !ref_info.all_mapping_qualities.is_empty() {
            let mut sorted_mapq = ref_info.all_mapping_qualities.clone();
            sorted_mapq.sort_unstable();
            let median_mapq = if sorted_mapq.len() % 2 == 0 {
                let mid = sorted_mapq.len() / 2;
                (sorted_mapq[mid - 1] + sorted_mapq[mid]) as f32 / 2.0
            } else {
                sorted_mapq[sorted_mapq.len() / 2] as f32
            };

            println!("Median mapping quality: {:.1}", median_mapq);
        } else {
            println!("Median mapping quality: N/A (no mapped reads)");
        }

        
        // Calculate percentages
        let total_bases = ref_info.length as f64;
        println!("\nCallability statistics for {}", ref_name);
        println!("Total length: {}", ref_info.length);
        println!("NO_COVERAGE: {:.2}% ({} bases)",
                 (stats.no_coverage as f64 / total_bases) * 100.0, stats.no_coverage);
        println!("LOW_COVERAGE: {:.2}% ({} bases)",
                 (stats.low_coverage as f64 / total_bases) * 100.0, stats.low_coverage);
        println!("EXCESSIVE_COVERAGE: {:.2}% ({} bases)",
                 (stats.excessive_coverage as f64 / total_bases) * 100.0, stats.excessive_coverage);
        println!("POOR_MAPPING_QUALITY: {:.2}% ({} bases)",
                 (stats.poor_mapping_quality as f64 / total_bases) * 100.0, stats.poor_mapping_quality);
        println!("CALLABLE: {:.2}% ({} bases)",
                 (stats.callable as f64 / total_bases) * 100.0, stats.callable);

        // Calculate average depth (from original code)
        let total_depth: u64 = ref_info.coverage.iter().map(|&x| x as u64).sum();
        let avg_depth = total_depth as f64 / ref_info.length as f64;
        println!("Average depth: {:.2}X", avg_depth);
    }
}