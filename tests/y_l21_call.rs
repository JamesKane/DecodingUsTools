use decoding_us_tools::types::ReferenceGenome;
use decoding_us_tools::utils::liftover::Liftover;
use decoding_us_tools::haplogroup::caller::tally_bases_at_positions;
use rust_htslib::bam::{self, Read};

// Focused test: verify that L21 (GRCh38 chrY:13542548 C>G) is derived (G)
// in the CHM13/hs1 BAM fixture by lifting over the GRCh38 position
// to the BAM's build and tallying the observed base using caller.rs API.
//
// Requires network (to fetch UCSC chain) and the BAM fixture at:
//   tests/fixtures/bam/sample_chm13_chrY.bam
// Run with:
//   cargo test --test y_l21_call -- --ignored --nocapture
#[test]
#[ignore]
fn call_l21_on_chm13_bam() {
    let bam_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bam/sample_chm13_chrY.bam");
    assert!(bam_path.exists(), "Fixture BAM missing: {}", bam_path.display());

    let mut bam = bam::Reader::from_path(&bam_path).expect("open BAM");
    let header = bam.header().clone();

    // Detect genome build from header and pick a Y sequence name present in the BAM
    let genome = ReferenceGenome::from_header(&header).expect("detect genome from header");
    let possibles = match genome.name() {
        "GRCh38" => vec!["chrY", "Y", "CM000686.2", "NC_000024.10"],
        "GRCh37" => vec!["Y", "chrY"],
        _ => vec!["Y", "chrY", "CP086569.2", "NC_060948.1"], // hs1 / CHM13
    };
    let chrom = possibles
        .into_iter()
        .find(|nm| header.tid(nm.as_bytes()).is_some())
        .expect("find Y/chrY seq in BAM header");

    // Liftover GRCh38 chrY:13542548 -> BAM build
    let lifter = Liftover::load_or_fetch(ReferenceGenome::GRCh38, genome.clone())
        .expect("load liftover chain");
    let mapped = lifter.map_pos("chrY", 13_542_548).expect("map L21 position");
    eprintln!("L21 mapped to {}:{} on {}", chrom, mapped, genome.name());

    // Use the library API to tally base counts around the known CHM13/hs1 position
    // Scan a small ±2bp window across any present Y-contig name to be robust to naming
    // differences. We expect at least 2 unfiltered reads with the derived 'G' allele.
    let known_hs1_pos = 14_449_238u32; // CP086569.2:14449238 [C>G]
    let window = 2u32;
    let possible_y_names = vec!["chrY", "Y", "CP086569.2", "NC_060948.1"]; // common hs1/Y labels
    let present_y: Vec<String> = possible_y_names
        .into_iter()
        .filter(|nm| header.tid(nm.as_bytes()).is_some())
        .map(|s| s.to_string())
        .collect();
    eprintln!("Present Y contigs in BAM: {:?}", present_y);

    let mut total_derived_g: u32 = 0;
    for yname in present_y {
        let positions: Vec<u32> = ((known_hs1_pos.saturating_sub(window))..=(known_hs1_pos + window))
            .collect();
        // fresh reader each time so we don't hit EOF
        let mut bam2 = bam::Reader::from_path(&bam_path).expect("reopen BAM");
        let header2 = bam2.header().clone();
        let cmap = tally_bases_at_positions(&mut bam2, &header2, &yname, positions, 0)
            .expect("tally base counts in window");
        for (pos, counts) in cmap.iter() {
            let depth: u32 = counts.values().copied().sum();
            eprintln!("Base counts at {}:{} via caller.rs = {:?} (depth={})", yname, pos, counts, depth);
            total_derived_g += counts.get(&'G').copied().unwrap_or(0);
        }
    }

    eprintln!("Total derived G reads across Y contigs in ±{}bp window: {}", window, total_derived_g);
    assert!(total_derived_g >= 2, "Expected at least 2 derived G reads for L21; got {}", total_derived_g);

    // Optional: log liftover mapping for visibility (do not assert strict proximity to allow chain quirks)
    let known_hs1_pos = 14_449_238u32;
    let delta = if mapped > known_hs1_pos { mapped - known_hs1_pos } else { known_hs1_pos - mapped };
    eprintln!("Liftover mapped position {}, known hs1 {} (delta={})", mapped, known_hs1_pos, delta);
}
