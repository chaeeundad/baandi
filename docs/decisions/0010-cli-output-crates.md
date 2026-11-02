# ADR 0010: Crates for CLI output (checksums, JSON)

Status: Accepted (2026-11-02, Claude Code). Applies ADR 0003 and 0005.

- `sha2` 0.10 (RustCrypto, MIT OR Apache-2.0): SHA-256 of the input TOML and of every UPF file, recorded in results.json for reproducibility. No hand-written hash.
- `serde_json` 1.x (MIT OR Apache-2.0): results.json and density metadata, using the serde derive already adopted in ADR 0005.
- Peak memory: not measured inside the engine (would need FFI or a platform crate). The CLI records array sizes (basis, grid, bands) and an estimated working-set size; measured peak RSS comes from `/usr/bin/time -l` in scripts/benchmark.sh.
- Build metadata (git commit, rustc version) is captured at compile time by a std-only build script with an "unknown" fallback.
