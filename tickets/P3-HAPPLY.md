TASK_ID: P3-HAPPLY
GOAL: |
  The reference driver's H application is 16.5-17.0 s on dev08 si128 MPI8 (rank max, SCF + verification solves, 14203 H columns) against QE h_psi 11.14 s for 82 calls (QE split: vloc_psi 7.98, fft 6.70, calbec 1.41, add_vuspsi 1.52). BAANDI has no split today (timings.fft_hpsi "not separately measured", parallel.stages zero). First measure, then cut.
READ FIRST: docs/performance/AGENT-BRIEF.md, docs/performance/attempts-ledger.md (many FFT/H-apply ideas were already tried on the old path; do not repeat a rejected idea without new evidence), docs/performance/backlog.md, docs/research/qe-study/qe-hpsi-and-fft.md, fft-call-structure.md, qe-gamma-real-path.md. Never read or request QE source.
TASK: |
  1. Instrument the H apply used by the reference driver (crates/dft-physics/src/distributed_real_gamma.rs, DistributedRealHamiltonian's DistributedRealOperator impl and what it calls) with lightweight per-solve accumulators recorded only at the full diagnostic level (BAANDI_DAVIDSON_STAGES=on): pack/unpack to the FFT grid, forward and inverse 3D FFT (split into local 1D/2D work and the transposes/alltoall), local potential multiply, kinetic, nonlocal projections (calbec-like GEMM) and back-application, plus columns processed per call and the pairing of two real bands per complex FFT if used. Serialize the totals into the reference driver's davidson_stages JSON (apply_h sub-object), rank max and min.
  2. Check the GEMM routes inside H apply (nonlocal products) the same way the Davidson rotations were found on the pure-Rust whole-output fallback (crates/dft-numerics/src/gemm_dispatch.rs routes; see REF_DAV_WIDE_ROTATE in switches.rs): any call on a pure_* route with a large share is a candidate; route it through the wide/dispatching entry behind a switch (default off).
  3. Mac checks only (user rule 2026-11-24): cargo check --workspace --tests; cargo test -p dft-physics; cargo test -p dft-engine --test refscf; cargo test -p dft-numerics switches; fmt; clippy -D warnings on touched crates. Prove default results bit-identical (instrumentation off and on) on a small fixture.
  4. Write docs/performance/lanes/RSCF/happly-dev08-request.md: one profiled si64/si128 MPI8 run per arm with the new split next to QE's h_psi lines, then A/A + crossed A/B for each new switch.
OWNED_FILES:
  - crates/dft-physics/src/distributed_real_gamma.rs
  - crates/dft-physics/src/   # only files reached from the reference driver's H apply; old-path behaviour bit-identical
  - crates/dft-solvers/src/ref_davidson.rs   # passing the apply_h split into the stage record only
  - crates/dft-solvers/src/ref_davidson/accounting.rs
  - crates/dft-cli/src/run.rs   # davidson_stages apply_h JSON only
  - crates/dft-numerics/src/switches.rs   # append new switches at the end of the list
  - crates/dft-numerics/tests/all/switches.rs   # count only
  - crates/dft-cli/tests/all/run.rs   # switch count only
  - docs/performance/lanes/RSCF/happly-dev08-request.md
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo test -p dft-numerics switches
  - cargo fmt --check && cargo clippy -p dft-physics -p dft-solvers -p dft-cli --all-targets -- -D warnings
  - (orchestrator, on dev08) scripts/check.sh full, profiled runs, A/B
ACCEPTANCE_CRITERIA: apply_h split recorded at full diagnostics, zero cost and bit-identical results when off; pure-route candidates identified and switchable; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
