//! Bootstrap command-line entry point; no calculation is implemented.
fn banner() -> &'static str {
    "rust-dft: P0 bootstrap (no DFT calculation implemented)"
}

fn main() {
    println!("{}", banner());
}

#[cfg(test)]
mod tests {
    #[test]
    fn bootstrap_banner() {
        assert!(super::banner().contains("no DFT calculation implemented"));
    }
}
