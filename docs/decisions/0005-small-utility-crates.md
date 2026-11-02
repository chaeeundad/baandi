# ADR 0005: Fixed-size geometry and utility crates

Status: Accepted (2026-11-02, Claude Code). Extends ADR 0004 (applies ADR 0003).

| Need | Crate | Version checked (crates.io) | License |
| --- | --- | --- | --- |
| 3x3 matrices and 3-vectors for lattice, reciprocal lattice, G-vectors | `nalgebra` | 0.35.0 (2026-05-24) | Apache-2.0 |
| Error types | `thiserror` | 2.0.21 | MIT OR Apache-2.0 |
| TOML input and typed config | `serde` (derive), `toml` | 1.0.229, 1.1.x | MIT OR Apache-2.0 |
| CLI argument parsing | `clap` (derive) | 4.6.7 | MIT OR Apache-2.0 |

Rules: pinned by Cargo.lock; declared once in `[workspace.dependencies]`; the transitive tree is license-audited when vendored (no GPL-only crates). nalgebra stays limited to small fixed-size geometry; large dense linear algebra uses faer (ADR 0004). Dense `nalgebra` types do not cross crate boundaries where a plain array or ndarray view suffices for performance-critical data.
