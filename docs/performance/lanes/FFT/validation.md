# LANE-FFT Phase-0 validation (Mac; correctness only)

Status: PROFILE_REQUEST. The approved instrumentation wiring is complete.
No performance conclusion is drawn from this host; no optimization is retained.

## Actual CLI equivalence and parser

Built `cargo build --offline --release -p dft-cli` (exit 0). Generated Si8 input:

```sh
python3 -I scripts/fast_screen_dav.py --case si8 --upf "$DFT_SI_UPF" \
  --threads 1 --seed 601 --amp 0.001 --emit-input "$PWD/target/scratch/fft-si8.toml"
```

Same binary/input in three fresh output directories (every CLI exit 0):

| Arm | FFT_PROFILE | PHASE_PROFILE | STAGE_TIMERS | DAVIDSON_STAGES |
|---|---|---|---|---|
| off | off | off | off | off |
| timed | on | on | on | off |
| counts | on | on | off | off |

Each switch above has the `BAANDI_` prefix. Command in every arm:
`target/release/dft-cli run target/scratch/fft-si8.toml --output-dir target/scratch/fft-final-<arm>`.

`python3 -I docs/performance/lanes/FFT/verify_cli.py target/scratch/fft-final-off target/scratch/fft-final-timed target/scratch/fft-final-counts`: exit 0.
Energy = -31.321776563202498 Ha in every arm; difference 0 Ha, 11 fine SCF
iterations, verification passed. Energies, eigenvalues, occupations, scientific
history and verification fields are compared by f64 bits. Only `timings` and
`wall_s` reporting fields are excluded. Density SHA-256 agrees exactly.
Integer FFT counters in timed/count-only output agree; the timed mode has
populated clocks and count-only FFT seconds are zero.

`python3 -I scripts/lane_fft_parse.py --json docs/performance/lanes/FFT/mac-counts-results.json`: exit 0, 144 rows. The retained file is actual,
unmodified CLI output. No Mac time in that file is performance evidence.
`cli-equivalence.json` records input, binary, UPF and source hashes plus exact
identity results. `mac-si8.toml` retains the input (external UPF path only).

## Kernel, transport and SCF correctness

- `cargo test --offline --release -p dft-numerics -p dft-solvers -p dft-engine --lib lane_fft -- --nocapture`: exit 0 (3 tests). Serial apply/kept-forward bits agree, density transforms do not contribute lane counters. M2 OFF/count-only/timed outputs agree at ranks 1/2/3/8 across cubic, skew-support and slab fixtures; the final fast gate also checks empty support. The 3-rank profile-gather test preserves exact u64 values above 2^53 and each rank's own map counts.
- `BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=off cargo test --offline --release -p dft-physics --test distributed_real_gamma`: exit 0, 11 tests. Includes multi-projector cubic/skewed/slab Hamiltonians, ranks 1–8, band counts 1/2/3/16, kept planes and density consumers.
- With the same flags, `cargo test --offline --release -p dft-engine --test distributed_scf m2c_stick_row_scf_matches_threads_only -- --nocapture`: exit 0. Existing synthetic cubic/skew SCF accuracy/rank-agreement checks pass. This is an independent existing SCF oracle, not a new claim of bit identity across different rank counts.
- `TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/FFT/test_parser.py`: exit 0, 3 tests. Rejects missing/disabled data, invalid rank sets, invalid metrics, unsupported multiple threads and unknown timer mode; checks units, per-call/per-map averages and count-only null times.

## Required gate and known failures

`CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-numerics dft-solvers dft-engine dft-cli`: exit 1. **Build, fmt and clippy pass**; 613 tests pass, 160 ignored, 3 known in-repository TMPDIR cache failures:

- `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`
- `pseudo::cache_resolution_and_no_network_side_effects`
- `pseudo::cached_fetch_checks_hash_and_provenance`

The optional load-report `ps` call is denied by the sandbox; the script still
runs every gate. No failed test is reported as passing. The unchanged README
validator previously exited 1 (one of 9 tests: an obsolete `standard profile
is unsupported` scope string). README, generator and reference fixtures are
untouched; this task adds instrumentation evidence, not a scientific result.

## Limits and handoff

One thread per rank only. Serial paired/pruned Gamma and M2 paths are covered;
complex/full-grid fallback FFT backends are not instrumented and zero calls
must be treated as unavailable. Map averages are same-rank bucket totals,
not individual-map traces. Timers perturb the measured code and byte/flop
counts are explicit logical lower-bound models, not hardware counters.
Real MPI/MKL execution and all performance interpretation belong to dev08.
No remote action or commit was made. See profile-request.md for the exact
commands, JSON paths, repetition protocol and return bundle.
