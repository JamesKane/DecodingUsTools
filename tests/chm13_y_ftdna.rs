use std::fs;
use tempfile::NamedTempFile;

#[test]
#[ignore]
fn chm13_chrY_ftdna_produces_r_fgc29071() {
    // This test requires network access to fetch the FTDNA tree and uses a BAM fixture.
    // Run with: cargo test --test chm13_y_ftdna -- --ignored --nocapture
    let bam_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/bam/sample_chm13_chrY.bam");
    assert!(bam_path.exists(), "Fixture BAM missing: {}", bam_path.display());

    let mut tmp_out = NamedTempFile::new().expect("create temp output file");
    let out_path = tmp_out.path().to_path_buf();

    // Use provider FTDNA and YDNA tree type
    let provider = decoding_us_tools::cli::TreeProvider::FTDNA;
    let tree_type = decoding_us_tools::utils::cache::TreeType::YDNA;

    // Run analysis. Reference file path is not needed for BAM; pass empty string.
    decoding_us_tools::haplogroup::analyze_haplogroup(
        bam_path.display().to_string(),
        String::new(),
        out_path.display().to_string(),
        1,   // min_depth
        20,  // min_quality
        tree_type,
        provider,
        false, // show_snps
    ).expect("haplogroup analysis finished");

    let contents = fs::read_to_string(&out_path).expect("read output file");
    let mut lines = contents.lines();
    let _header = lines.next().expect("has header line");
    let first = lines.next().expect("has at least one result line");
    let first_field = first.split('\t').next().unwrap_or("");

    assert_eq!(first_field, "R-FGC29071", "Top haplogroup should be R-FGC29071, got: {}", first_field);
}
