use std::process::Command;
fn capture(program: &str, args: &[&str]) -> String {
    Command::new(program)
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_owned())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".into())
}
fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // Worktrees keep HEAD and branch refs in different git directories.
    for name in [
        "HEAD".to_owned(),
        capture("git", &["symbolic-ref", "-q", "HEAD"]),
        "packed-refs".to_owned(),
    ] {
        let path = capture("git", &["rev-parse", "--git-path", &name]);
        if path != "unknown" {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    println!(
        "cargo:rustc-env=DFT_GIT_COMMIT={}",
        capture("git", &["rev-parse", "HEAD"])
    );
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".into());
    println!(
        "cargo:rustc-env=DFT_RUSTC_VERSION={}",
        capture(&rustc, &["--version"])
    );
}
