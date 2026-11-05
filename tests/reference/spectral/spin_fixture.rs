//! Reuse the exact P13B test generator; print only to a temporary-file consumer.
#[allow(dead_code)]
#[path = "../../../crates/dft-cli/tests/common/mod.rs"]
mod common;
fn main() {
    print!("{}", common::gaussian_with(1.0, 1.0, false, 401)
        .replace("SLA PZ", "SLA PW")
        .replace("element=\"Xx\"", "element=\"H\""));
}
