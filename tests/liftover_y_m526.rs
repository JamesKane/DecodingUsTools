use decoding_us_tools::types::ReferenceGenome;
use decoding_us_tools::utils::liftover::Liftover;

// Focused liftover test for the Y-SNP M526
// Issue-provided coordinates:
// - GRCh38 chrY:21389038 [A>C]
// - CHM13v2 chrY:22264299 [A>C]
//
// This test downloads UCSC chain(s) on first run and caches them; mark ignored by default.
// Run with:
//   cargo test --test liftover_y_m526 -- --ignored --nocapture
#[test]
#[ignore]
fn liftover_m526_grch38_to_chm13v2() {
    let lo = Liftover::load_or_fetch(ReferenceGenome::GRCh38, ReferenceGenome::CHM13v2)
        .expect("load liftover GRCh38->CHM13v2");
    let mapped = lo.map_pos("chrY", 21_389_038);
    assert!(mapped.is_some(), "M526 position did not liftover");
    let mapped_pos = mapped.unwrap();
    eprintln!("M526 GRCh38 chrY:21389038 -> CHM13v2 chrY:{}", mapped_pos);
    assert_eq!(mapped_pos, 22_264_299, "Expected CHM13v2 chrY:22264299 for M526");
}
