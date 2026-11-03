# FFT backend comparison: P10B stage 1

Status: implementation and Mac declaration checks; dev08 numerical and timing
acceptance pending. No measured speedup for this implementation is claimed.
The default remains pure Rust. The optional user-selected backend is the existing
ADR decision, not a new recommendation from measurements that have not run.

## Prior evidence motivating implementation

ADR 0019 records one-core dev08 (2 x Xeon Gold 6230, MKL 2026.1) microbenchmarks:
MKL/rustfft speed ratios 3.9 at 24^3, 3.0 at 25^3, 2.1 at 27^3, 2.8 at 36^3,
4.1 at 48^3; GEMM 2.5–4.2, EVD approximately equal. These are historical kernel
measurements, not this adapter's results and not a complete-SCF speedup.
The profile motivating the work is FFT ~54% and GEMM ~9% (ADR 0019).
The existing optimization log remains authoritative; no previously rejected
experiment is repeated and no new optimization is marked accepted here.
The QE call-structure note is context only: stage 1 follows the ADR full-transform
path; note-derived batched/fused optimizations await P10C and new measurements.

## Implementation and evidence matrix

| Gate | Evidence | State |
| --- | --- | --- |
| Pure Rust default build/test/fmt/clippy | `scripts/check.sh`, Mac, four PASS lines, exit 0 | Passed locally |
| MKL declaration check | `cargo check -p dft-numerics --features mkl`, exit 0 | Passed locally |
| All workspace targets with MKL | `cargo check --workspace --all-targets --features dft-numerics/mkl`, exit 0 | Passed locally |
| Unsafe documentation and feature linting | `cargo clippy --workspace --all-targets --features dft-numerics/mkl -- -D warnings`, exit 0 | Passed locally |
| README evidence consistency | generator plus `test_readme_figures.py`, 2 tests, exit 0; no generated diff | Passed locally; no MKL reference added |
| Dev08 shell syntax | `bash -n crates/dft-mkl/tests/dev08.sh`, exit 0 | Passed locally; execution pending |
| Native ABI/link | `abi_probe.c`, linked release tests | Pending dev08 |
| Forward/inverse accuracy, round trip, Parseval, pruned contracts | linked FFI tests, MKL API tests, `fft_compare` | Pending dev08 |
| GEMM correctness | all four adjoint combinations vs scalar oracle; existing numerics/solver tests | Pending dev08 |
| Thread bit identity | FFT 1/2/4 workers; synthetic engine/physics tests; real SCF 1/2/4/8 budgets | Pending dev08 |
| Full Hamiltonian/SCF equivalence and timing | dev08 protocol, fresh Rust vs MKL binaries | Pending dev08 |
| Real QE agreement | fresh MKL results vs existing provenance-backed QE bundles | Pending dev08 |
| General recommendation | P10G multi-system suite and above gates | Pending; defaults unchanged |

## Harness and interpretation

Run the exact protocol in [howto-mkl](../howto-mkl.md). The harness lives under
`crates/dft-mkl` to stay within ticket ownership. Its dev-dependencies on the
existing numerics/ndarray crates provide the independent Rust comparison; they
add no production cycle or third-party package. The lockfile updates only local
package edges. The script selects numerics explicitly when enabling its feature
alongside downstream packages; Cargo rejects a transitive feature on a lone CLI
or engine package.

Grids: cube edges 24, 25, 27, 32, 36, 48, 64, 96, 128, 192. Report skips and byte
estimates. Plan and pattern construction are separate from forward/inverse
execution. Timings include copies through the safe adapter. The Rust pruned and
MKL full-plus-mask measurements use the same radius-n/4 reciprocal support;
MKL mask membership is precomputed outside timing. The low-level MKL benchmark
matches the production arithmetic/copy path; ndarray shape checks add small
untimed overhead relative to the physics-facing wrapper. Kernel timing pins to
one core. Multi-core thread correctness also runs in the ordinary linked suites.

The harness always transforms a fixed finite input, rather than repeatedly
feeding an unscaled inverse/forward result back into itself. Three independent
runs establish spread. Every element is checked against rustfft with the strict
1e-13 relative gate, without a near-zero floor. Unit tests cover both directions;
large-grid benchmark comparison is forward. Round-trip absolute error and
relative Parseval defect each must be <=1e-13. Independent compensated direct DFT
covers all coefficients of small noncubic grids and four coefficients of every
large grid. It uses f64 phase evaluation, so it is not a long-double or arbitrary-
precision oracle. Keep raw maximum errors if the strict relative gate fails.

No allocator instrumentation or complete new-system benchmark evidence is
available locally. The conservative in-place buffer choice adds two copies and
may lose some raw vendor speedup; out-of-place and batched variants require
measurement. The architecture keeps stage 2 open without changing P10C's files.
Sequential-library linkage is guaranteed by construction, but native runtime
behavior, numerical accuracy and determinism cannot be established by type-check.

## Recommendation rule

Do not choose automatic selection or advertise an end-to-end speedup until the
returned logs establish accuracy, fixed-point/charge/eigenpair checks and full
SCF timing on multiple systems. If correctness gates fail, retain the Rust default
and investigate the units/normalization/support/rounding cause. If only Si wins,
record that scope and keep defaults unchanged. No optimization-log acceptance or
README MKL result should be written before those measurements exist.

## Local default-backend regression numbers

`cargo test -p dft-numerics -- --nocapture` on the Mac: exit 0; 16 unit tests
passed, one unit benchmark and one integration benchmark ignored. These are Rust
backend regression results, not MKL measurements. Dimensionless FFT errors:

| Grid | Reciprocal round trip | Real round trip | Parseval relative defect |
| --- | ---: | ---: | ---: |
| 4 x 5 x 6 | 4.2533482651865464e-16 | 4.5899321840152809e-16 | 2.3544627737869695e-16 |
| 15 x 20 x 24 | 1.0327443942301632e-15 | 1.0007328703583267e-15 | 2.2735942985868243e-15 |
| 16 x 16 x 16 | 7.1101711388830003e-16 | 7.9007922800082654e-16 | 3.5528748979513846e-15 |
| 27 x 27 x 27 | 1.7644611924521969e-15 | 1.7433712570974990e-15 | 3.6967131021476877e-16 |

Plan-reuse differences were exactly zero on all four grids. Existing faer GEMM
oracle errors for identity/identity, identity/adjoint, adjoint/identity and
adjoint/adjoint were respectively 4.47545209131181e-16, 4.440892098500626e-16,
4.503738912626819e-16 and 2.482534153247273e-16 (dimensionless). No new MKL evidence
was added to the README reference dataset.
