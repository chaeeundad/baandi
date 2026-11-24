TASK_ID: P3-SETUP2
GOAL: |
  Reference-driver setup is 3.38 s on dev08 si128 MPI8 (one run, rank 0) against QE init_run 0.83 s (wfcinit 0.68). Breakdown before REFSTART_MKL (si128, engine 4.81 s): atomic_orthonormalize_complete 1.50, atomic_h_apply_gram 0.95, nonlocal_diagonal 0.58, atomic_orbital_fill 0.40, local tables + structure factors 0.41, atomic_diagonalize_rotate 0.21, density/field FFT plans 0.24/0.17, basis 0.13. REFSTART_MKL moved the start Gram/rotations to MKL (setup 4.83 -> 3.36 s). Cut the rest.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (REFSTART_MKL, P3-SETUP, CholQR rows), docs/design/reference-scf.md sections "Start" and "Setup and output cost", docs/research/qe-study/qe-scf-cycle.md and qe-stage-by-stage.md. Never read QE source.
TASK: |
  1. Re-measure the setup split on the Mac fixtures and account each item (FLOPs, bytes, collectives) as in P3-DAVREM; find what scales badly: e.g. atomic_orthonormalize_complete (two CholQR passes on the atomic block, its Gram reductions, the MGS fallback path), atomic_h_apply_gram (H applied to all atomic orbitals: number of columns vs needed bands), nonlocal_diagonal (built for all G on every rank?), orbital fill (radial interpolation per G per orbital: table lookup vs recomputation), FFT plan creation (cached across runs?).
  2. Remove excess behind switches (default off; bit-identical ones may be proposed default on). Keep the start state identical within rounding and the reference driver's trajectory counts unchanged on fixtures.
  3. Mac checks only: cargo check --workspace --tests; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings.
  4. docs/performance/lanes/RSCF/setup2-dev08-request.md: ONE profiled run per arm, si64/si128/elongated/rskew MPI8, setup_s table next to QE init_run/wfcinit.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/atomic_start.rs   # only if refstart calls into it
  - crates/dft-engine/tests/refscf.rs   # new tests only
  - crates/dft-numerics/src/switches.rs   # append at the end (P3-DAVREM also appends)
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/design/reference-scf.md   # "Setup and output cost" section only
  - docs/performance/lanes/RSCF/setup2-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: setup items accounted; removable work behind switches; start state within rounding; dev08 one-run request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
