use decoding_us_tools::types::ReferenceGenome;
use decoding_us_tools::utils::liftover::Liftover;

// Focused liftover test for the well-known Y-SNP L21 (aka R-L21 / R-S145),
// to validate GRCh38 <-> CHM13v2 mappings on chrY.
// Coordinates provided by the issue description:
// - GRCh38 chrY:13542548 [C>G]
// - CHM13v2 chrY:14449238 [C>G]
//
// This test fetches UCSC chain files on first run and caches them under the
// standard DecodingUs cache directory. Because it requires network access,
// it is marked as ignored by default. Run with:
//   cargo test --test liftover_y_l21 -- --ignored --nocapture

#[test]
#[ignore]
fn liftover_l21_grch38_to_chm13v2_and_back() {
    // Forward: GRCh38 -> CHM13v2
    let lo_fwd = Liftover::load_or_fetch(ReferenceGenome::GRCh38, ReferenceGenome::CHM13v2)
        .expect("load liftover GRCh38->CHM13v2");
    let mapped = lo_fwd.map_pos("chrY", 13_542_548);
    assert!(mapped.is_some(), "Liftover did not map the L21 GRCh38 position");
    let mapped_pos = mapped.unwrap();
    eprintln!("GRCh38 chrY:13542548 -> CHM13v2 chrY:{}", mapped_pos);

    // Reverse: CHM13v2 -> GRCh38 using whatever forward mapping produced
    let lo_rev = Liftover::load_or_fetch(ReferenceGenome::CHM13v2, ReferenceGenome::GRCh38)
        .expect("load liftover CHM13v2->GRCh38");
    let mapped_back = lo_rev.map_pos("chrY", mapped_pos);
    assert!(mapped_back.is_some(), "Reverse liftover did not map the CHM13v2 L21 position");
    let mapped_back_pos = mapped_back.unwrap();
    eprintln!("CHM13v2 chrY:{} -> GRCh38 chrY:{}", mapped_pos, mapped_back_pos);
    assert_eq!(mapped_back_pos, 13_542_548, "Expected GRCh38 chrY:13542548 after reverse mapping");
}
