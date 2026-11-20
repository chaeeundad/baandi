# Phase 0 correctness and plumbing evidence

Status: **PROFILE_REQUEST**. No optimization remains in the patch. The Mac
was used for compilation and identity/plumbing checks after the revised user
instruction. No Mac time supports any performance conclusion. The local raw
files necessarily contain clocks for testing the clock-output contract, but
tracked `correctness.json` deliberately contains no timing values.

## Identity

Baseline commit `db25e70176bd5d111dab6fe1dc2a110766fcb005`; baseline and current
release CLI and UPF hashes are stored in `correctness.json`. These are internal
same-code comparisons, not QE references or absolute-accuracy measurements.

Five real-PP inputs, seed 601, fractional displacements +/-0.001, one core:
Si8, Si32, skew Si, slab Si, shifted k222 Si. Each compares the baseline binary
with current default-off, PHASE_PROFILE-only, lane count-only and lane-clock
modes. **20 comparisons across 25 completed arms**: full density bytes, energy,
eigenvalues, occupations, scientific histories, certificates, verification
metrics/schedules and iteration/work records match exactly (clocks/metadata
excluded). Maximum energy difference **0 Ha**. Four inputs converge. Slab
exhausts the same iteration budget in all five arms (exit 3, no verification
attempts); it is an unchanged failure control, not a converged slab reference.
The synthetic slab below does converge.

Release identity matrix: 48 SCF executions / 30 paired mode comparisons,
cubic/skew/slab synthetic nonlocal cells, bands 3/6, k blocks 1/2, two forced
failed verification attempts, plus simulated MPI 2/3 ranks. Entire density,
eigenvalues/residuals/orbitals and scientific history/certificate/report fields
are identical off/count-only/clocks. Profile-disabled + lane-enabled produces
no records. Count-only records have all new seconds zero. Timed records test
finite clock plumbing only. No comparison or assertion uses a speed ratio.
Simulated ranks test arithmetic identity at each rank count and agreement among
ranks; they are not MPI performance measurements.

The transport specification test covers unequal per-rank record counts and
empty records using three simulated ranks. The parser's five specification
tests cover per-call critical rank, whole-run critical rank (never sum maxima),
count-only/schedule separation, invalid/nonfinite inputs and empty records.
It also parses actual Mac Si8 and multi-k/slab `results.json` files; Si8 has two
attempt events and one scheduler event. All record counts match the existing
verification report. `correctness.json` gives the other case counts.

## Commands and results

- `CARGO_INCREMENTAL=0 cargo build --release -p dft-cli`: exit 0
  (`target/scratch/profile-build-v2.log`).
- `CARGO_INCREMENTAL=0 cargo test --release -p dft-engine --lib
  verification_lane_profile_identity -- --ignored --test-threads=1`: exit 0,
  1 matrix test passed (`target/scratch/profile-identity-final.log`).
- `python3 -I docs/performance/lanes/VERIFY/validate.py --binary
  target/release/dft-cli --baseline target/scratch/dft-cli-baseline --upf
  ../pseudo/Si.pz-vbc.UPF --out target/scratch/phase0-identity`: initial driver
  exit 1 after the first three cases passed, because it initially demanded
  convergence from the already-failing baseline slab.
- The driver now checks equality of successful **or exhausted** outcomes.
  `python3 -I docs/performance/lanes/VERIFY/validate.py --binary
  target/release/dft-cli --baseline target/scratch/dft-cli-baseline --upf
  ../pseudo/Si.pz-vbc.UPF --out target/scratch/phase0-extra --cases
  si-slab,si-k222s`: exit 0. All five arms for both cases match; the original
  successful cases are preserved under `phase0-identity` without unnecessary
  reruns. The revised driver reproduces the full matrix in a fresh directory.
- `python3 -I docs/performance/lanes/VERIFY/test_parser.py`: exit 0, 5 tests.
- `python3 -I scripts/lane_verify_parse.py
  target/scratch/phase0-identity/si8-gamma/clock/results.json --json-out
  target/scratch/parser-plumbing.json --markdown-out
  target/scratch/parser-plumbing.md`: exit 0. Empty slab and multi-k files also
  parse (same command with `phase0-extra/{si-slab,si-k222s}/clock/results.json`).
- `git diff --check`: exit 0.
- Required fast gate: result appended below when complete.

Development failures retained: first release build exit 101 used nonexistent
`PlaneWaveBasis::len`; fixed to `vectors().len`. First matrix test exit 101
assumed every timed scheduler call must exceed one clock tick; corrected to
accept finite nonnegative scheduler times (attempts still require nonzero
elapsed time). No scientific tolerance was changed. An explicit Python bytecode
compilation command was blocked by the sandbox's external cache path; AST and
actual parser checks pass without writing bytecode outside the workspace.

## Scope and limits

The revised Phase 0 request requires the new profiler module, its `profile.rs`
export, minimal post-run CLI JSON/gather glue and `scripts/lane_verify_parse.py`,
in addition to the original lane's session/tests/switch files. No numerical
algorithm file, README result, artifact registry or dashboard was changed.
This is instrumentation plumbing, not a new published physics result.
The orchestrator owns registry/dashboard updates and dev08 execution.

New clocks and counters are off by default. No new MPI collective occurs in a
numerical loop; enabled profile records are gathered after numerical execution.
Real MPI+MKL collection and performance attribution remain untested locally and
must use `profile-request.md`. Clock reads and post-run collection are profiling
overhead, not accepted production performance. Stop after this phase.

Final required gate:
`CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" CARGO_INCREMENTAL=0
scripts/check.sh fast dft-engine dft-numerics dft-cli` exited **1**.
Build, fmt and clippy **PASS**. **471 passed / 3 failed / 136 ignored**, 51
binaries. Failures are the known in-repository TMPDIR pseudopotential-cache
restrictions: `pseudo::cache_resolution_and_no_network_side_effects`,
`pseudo::cached_fetch_checks_hash_and_provenance`, and
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
The gate's ancillary `ps` load inventory is sandbox-blocked; all four actual
gates still ran. Full log: `target/scratch/phase0-check-final.log`.

The first gate also found `chunks_exact_to_as_chunks` in the new profile
transport. It was fixed to the compiler-recommended `as_chunks` API; the final
gate's transport test and clippy pass. Final release rebuild exits 0
(`phase0-release-final.log`); a fresh Si8 five-arm identity/plumbing check exits
0 (`phase0-final-smoke.log`), and the final parser command/tests exit 0 (5 tests).
Both the full-matrix binary hash and final binary hash are preserved in
`correctness.json`. There is no accepted performance result and no default flip.
