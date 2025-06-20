use std::path::PathBuf;

fn main() {
    let generated_dir = PathBuf::from("src/generated");
    std::fs::create_dir_all(&generated_dir).unwrap();

    protobuf_codegen::Codegen::new()
        .input("proto/vg/gam.proto")
        .include("proto")
        .out_dir(&generated_dir)
        .run()
        .expect("Protobuf codegen failed");

    // Tell cargo to rerun if proto files change
    println!("cargo:rerun-if-changed=proto/vg/gam.proto");
}