# Direct dependencies

Exact resolved versions below were read from Cargo.lock; license expressions
were checked in the corresponding vendored Cargo.toml files. The lockfile,
workspace manifests, vendor, and Cargo source configuration are unchanged by
P1A. Use `CARGO_NET_OFFLINE=true` for verification. Rust 1.99.1 is pinned in
rust-toolchain.toml (edition 2024, resolver 3).

| External direct dependency | Locked version | License | Workspace consumers |
| --- | --- | --- | --- |
| nalgebra | 0.35.0 | Apache-2.0 | dft-geometry, dft-basis |
| thiserror | 2.0.21 | MIT OR Apache-2.0 | dft-geometry, dft-numerics, dft-basis, dft-cli |
| ndarray | 0.17.2 | MIT OR Apache-2.0 | dft-numerics |
| num-complex | 0.4.6 | MIT OR Apache-2.0 | dft-numerics, dft-basis |
| rustfft | 6.4.1 | MIT OR Apache-2.0 | dft-numerics (also used through ndrustfft) |
| ndrustfft | 0.6.2 | MIT | dft-numerics |
| clap | 4.6.7 | MIT OR Apache-2.0 | dft-cli |
| serde | 1.0.229 | MIT OR Apache-2.0 | dft-cli |
| toml | 1.1.8+spec-1.1.0 | MIT OR Apache-2.0 | dft-cli |

There are currently no inter-workspace path dependencies or build dependencies.
dft-units uses std only. faer was selected in ADR 0004 for future dense linear
algebra but is not declared in the current workspace; it is not claimed as a
current dependency. LAPACK and the Accelerate source selectors were spike-only.
The hand-written Accelerate adapter and native tests were removed in P1A.

Adoption/maintenance review: ADR 0003 records general library policy and release
checks; ADR 0004 records the measured ndrustfft/rustfft backend choice and Apple
Silicon results; ADR 0005 records nalgebra and utility-library checks. Current
foundation crates build and test on the installed aarch64-apple-darwin toolchain;
no new releases, dependencies, native installs, or vendor changes are introduced.
The historical backend spike documents internal allocation and performance
limits. This direct-dependency inventory is not a fresh transitive license audit.
The project's own LICENSE remains undecided.
