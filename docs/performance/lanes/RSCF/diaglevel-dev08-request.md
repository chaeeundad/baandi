# I26-2-DIAGLEVEL: diagnostic collection and delivered-wall A/B

Status: implementation and Mac validation; dev08 measurement pending through
the master session. No remote action, QE source access, numerical default
promotion or speed claim. The reference driver remains opt-in. This ticket
implements the requested lightweight diagnostic default, not a new solver or
convergence policy. Strict 1e-10 Ha and standard 5e-9 Ha classes are unchanged.

## Change and trace contract

The reference driver resolves `RunConfig::current().on(DAVIDSON_STAGES)` once
and passes it to ordinary, redo, retry and verification solves. Off avoids
the per-band Ritz matrix and inner stage clocks; on records all bands for all
corrections. Ritz matrices move to driver records rather than being cloned.
Work counters and coarse phase clocks are retained in either mode. Full-mode
inner clocks now reach JSON, with a run aggregate covering verification too.
The historical duplicate iteration paths are preserved deliberately, so full
traces remain readable and their entire write cost stays charged.
See [Diagnostics level](../../../design/reference-scf.md#diagnostics-level).

`refscf-check` and `refscf-map` consume QE binary density/stdout, not the driver
JSON. No migration of those readers is necessary: existing field locations,
types and complete full-mode rows remain intact. The additive descriptor
`reference-scf-diagnostics-1` distinguishes full from intentionally empty Ritz
arrays. Null inner clocks mean unmeasured. No completion/verification key is
conditional on the diagnostic setting.

## Mac reproduction and byte definitions

Run from the worktree root, with an external UPF (never committed):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
mkdir -p target/tmp/diaglevel
TMPDIR="$PWD/target/tmp" cargo build --release -p dft-cli
python3 -I crates/dft-cli/tests/diaglevel_measure.py \
  --binary target/release/dft-cli --upf ../pseudo/Si.pz-vbc.UPF \
  --output target/tmp/diaglevel/paired \
  --compare-baseline target/tmp/diaglevel/baseline
```

Before edits, build the parent and copy its executable inside the worktree to
`target/tmp/diaglevel/baseline-cli`; invoke the same script with `--baseline`
and `--output target/tmp/diaglevel/baseline`. Omit `--compare-baseline` when
only checking on/off. The script rejects stale output directories.

Inputs reproduce the existing reference-driver conventional Si8/Si16/Si32
fixtures (1x1x1, 2x1x1, 2x2x1, deterministic 1e-4 fractional displacement),
30 Ry cutoff, occupied+4 bands, Gamma, LDA-PZ, atomic start, beta 0.7/history 8,
standard 5e-9 Ha, independent verification on. Each case uses identical input
bytes, one binary and one serial process; only the diagnostic environment
override changes. Sequence: lightweight/full, then full/lightweight. Builds
must finish before timed pairs. These two Mac pairs are observations only,
not a significance test or evidence of general speed improvement.

The runner compares density bytes, complete scientific history, all energy
terms/eigenvalues/occupations and every verification attempt exactly, excluding
only clocks and optional diagnostics. It also checks every full trace has
`davidson_iterations` rows of `bands` values, and all light traces are empty.

- Ritz payload bytes: 8 times the number of recorded f64 changes, counting
  each SCF/verification call once, excluding Vec overhead and JSON aliases.
- Ritz JSON value bytes: compact UTF-8 JSON of every `ritz_changes_ha` value
  in results.json, **including all aliases**. Empty arrays cost two bytes each;
  this definition excludes key names and pretty-print whitespace.
- Output bytes: actual file sizes for every delivered file in the output
  directory plus redirected stdout/stderr; results.json is reported separately.
  Inputs, executable and the measurement harness's evidence.json are excluded.
- External wall: launch through child exit and stdout/stderr close; includes
  startup, every failed verification/recovery, hashing, density/log/JSON writes
  and all required console output. Internal phase times are not the verdict.

## Measured Mac results

Stored [machine-readable evidence](../../../../crates/dft-cli/tests/fixtures/diaglevel-mac.json); raw outputs and logs remain in `target/tmp/diaglevel/` (scratch, not committed).

Source base: `c2a143d524abcfafa860e78c3817d054b2589454`. Host: `macOS-27.0.1-arm64-arm-64bit`.
Candidate binary SHA-256: `a9aca4bd8ad8e61e96699871784539265ef78d1bdedc50600ddbafe41501f87e`.
Pre-edit binary SHA-256: `70317e333ab4655ddfd5d7eed90e827171e20ad082ce57b52f7993487875ab13`.
External Si UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

All 12 runs converged in 12 SCF steps, passed independent verification, and
passed `scripts/parse_results.py` with the candidate binary hash. Every
scientific value and density byte is identical on/off **and to the pre-edit
reference executable**: energy delta exactly 0 Ha; density max delta exactly
0 e/Bohr^3. Full Ritz matrices also match the pre-edit traces exactly.

| Case | Energy (Ha), both modes | Verification attempts (failed + passed) | Final Coulomb error (Ha) |
|---|---:|---:|---:|
| si8 | -31.321801997468832 | 2 + 1 | 1.776502657457e-09 |
| si16 | -63.129238894587246 | 1 + 1 | 3.087484818755e-09 |
| si32 | -126.703985958871115 | 3 + 1 | 1.029302384745e-09 |

Failed verifications and subsequent recovery remain present, with identical
counts/values. They are inside every wall below; no attempt was subtracted.

| Case | Ritz f64 payload bytes off / on | Ritz JSON value bytes, all copies, off / on | results.json bytes off / on (range across two samples) | Total output bytes off / on (range) |
|---|---:|---:|---:|---:|
| si8 | 0 / 9,920 | 56 / 48,411 | 212,813–212,858 / 313,970–314,102 | 593,652–593,699 / 694,818–694,943 |
| si16 | 0 / 17,856 | 58 / 84,414 | 219,116–219,149 / 384,616–384,738 | 974,275–974,308 / 1,139,775–1,139,897 |
| si32 | 0 / 52,224 | 62 / 236,815 | 247,593–247,650 / 684,788–684,837 | 1,751,253–1,751,308 / 2,188,449–2,188,497 |

| Case | Pair 1: off then on (s) | Pair 2: on then off (s) |
|---|---:|---:|
| si8 | 1.928281, 1.761756 | 1.330987, 1.738064 |
| si16 | 5.136840, 4.406416 | 4.159584, 3.004397 |
| si32 | 16.732941, 13.544444 | 17.211211, 25.832670 |

**No Mac speed conclusion:** measured 1-minute load was 65.95–108.11; only
two crossed pairs per case. Candidate pairs ran after this task's builds
finished, but other host work was not controlled. The pre-edit capture
overlapped the candidate build, so its wall is not used in any comparison.
Output-size ranges reflect clocks/provenance text, not changed science.

The synthetic skew-cell CLI test additionally covers input/default selection,
environment-over-input precedence, full trace aliases and verification
identity. A synthetic tridiagonal solver test checks zero-capacity Ritz storage
and zero stage clocks when off, exact eigenvectors/HX/flags/work identity,
and both exhausted and converged calls. These do not establish broad material
performance or real-MPI behavior.

## Verification commands and remaining gate failures

Commands run with `PATH=/opt/homebrew/opt/rustup/bin:$PATH` and scratch inside
`target/tmp`; logs are in `target/tmp/diaglevel/`.

| Command | Exit | Result |
|---|---:|---|
| `cargo build --release -p dft-cli` (pre-edit) | 0 | Saved parent executable for identity comparisons |
| `cargo test --release -p dft-cli --test reference_driver broyden_reference_inputs -- --nocapture` | 0 | Synthetic full/light/input/environment contract and science identity |
| `cargo test --release -p dft-solvers --lib diagnostic_tests -- --nocapture` | 0 | Optional allocation/clocks, exhaustion/convergence, exact science/work |
| `DFT_UPF=../pseudo/Si.pz-vbc.UPF cargo test --release -p dft-cli --test refscf_map -- --include-ignored --nocapture` | 0 | 4 tests; existing QE-reader/map contracts preserved |
| `cargo test --workspace --bin dft-cli refscf_ -- --nocapture` | 0 | 3 reader/parser/replay tests |
| `diaglevel_measure.py` baseline and crossed-pair commands above | 0 | 3 parent captures + 12 candidate deliveries; exact on/off/parent science |
| `scripts/parse_results.py` for all 12 candidate outputs, expected binary hash | 0 each | Completion, provenance and clock consistency accepted |
| `bench_scf.py --cases gamma,k222s,si8 --threads 1 --repeats 1 --skip-legacy-si8-anchor`, pre-edit capture then candidate `--baseline ... --require-identical` | 0 | Legacy/default standard driver: exact energy and density identity on three cases |
| `bench_scf.py --cases gamma --threads 1 --repeats 1 --baseline docs/performance/baseline.json --require-identical` | **3** | Historical strict baseline is incompatible with current standard default; not an identity oracle |
| `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/tmp" scripts/check.sh full` | **1** | build/fmt/clippy PASS; 829 tests pass, 3 fail, 191 ignored; 348 s |
| `python3 -I scripts/make_readme_figures.py` | 0 | No tracked README/figure diff |
| `python3 -I scripts/tests/test_readme_figures.py` | **1** | 8 pass, existing obsolete scope-text expectation fails |

The three full-gate failures are unchanged pseudo-cache tests:
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`,
`pseudo::cache_resolution_and_no_network_side_effects`, and
`pseudo::cached_fetch_checks_hash_and_provenance`. They require caches outside
repositories, while this task explicitly requires scratch inside the worktree.
The first was rerun alone (exit 101), reproducing
`pseudopotential cache must be outside repositories`; both pseudo source/test
files have no diff from the pre-edit HEAD. No cache policy or test was weakened.
The gate's `ps` load probe is also denied by this sandbox; it continues to the
four gates. The README test still expects `standard profile is unsupported`,
contradicting its already-updated generated text; those files are unchanged.
Full verification is therefore **not green**, despite the diagnostic tests
passing. Dashboard/session registry, ledger/backlog and new README result types
are outside this ticket's owned paths and were not changed. No remote results
or default speed adoption are claimed.

The legacy benchmark captures are stored in the evidence JSON and scratch
`legacy-baseline.json` / `legacy-candidate.json`. Both use the absolute UPF
path and `TMPDIR="$PWD/target/tmp"`; the standard run explicitly skips only
the historical strict/complex Si8 absolute-energy anchor. Parent energy and
density comparisons remain enabled. No mixed-profile timing ratio is used.

## dev08 pre-registration (master session only)

Cases: established Si64 and Si128 inputs, **MPI np8 on both arms**, one thread
per rank on the same eight physical cores/socket. Use identical source build,
compiler/backend, cell, positions, cutoff, explicit FFT grid, UPF bytes, bands,
start, seed, threshold, binding, output options and verification. Record input,
UPF and binary SHA-256, reviewed source commit and patch hash, rustc/features,
MPI version, CPU binding, node/load and filesystem. Keep the existing solver,
preconditioner and reference-driver selection unchanged in both arms.

Build a reviewed candidate in its own target directory with the established
pure MPI command (`cargo build --release -p dft-cli --features mpi`); record
the real exit code. If the established benchmark uses MKL, keep that backend
fixed for both arms and record the exact build instead of mixing backends.
No installs or server configuration changes are requested.

1. A = `BAANDI_DAVIDSON_STAGES=on` (full), B = `off` (lightweight). Leave
   unrelated instrumentation off in both arms. Warm up each arm once, retaining
   all files and statuses separately; warmups are explicitly not timed samples.
2. **A/A floor first:** 20 pairs per case of the same full arm under separate
   labels/directories, using the exact delivery command/output below. Calculate
   `d_AA = log(wall_A2 / wall_A1)`, its sample mean and `SE=sample_std/sqrt(20)`.
   Define floor `F=max(abs(mean(d_AA)), SE(d_AA))`. Keep every raw wall and
   output byte count; record thermal/load drift. A/A must have identical science.
3. **Crossed A/B:** preregister 30 pairs per case: AB, BA, AB, BA, ... . Each
   pair uses fresh output directories and the identical input, with no build,
   profiling session, or other benchmark overlapping. Record `d=log(wall_B /
   wall_A)` consistently regardless of execution order. Do not stop early on
   favorable significance or silently delete slow/failed samples.
4. Report raw paired walls, each arm's median, `mean(d)`, SE, geometric ratio
   `exp(mean(d))`, and `abs(mean(d))/SE`. A measurable reduction requires
   `mean(d)<0`, `abs(mean(d))/SE >= 3`, and `abs(mean(d)) >= 3*F`, with no
   correctness regression. If SE is zero, report that explicitly and do not
   fabricate a finite ratio; establish a resolvable timing floor first. Below
   threshold is inconclusive, not a speedup. Do not pool Si64 and Si128 samples.
5. Check identical density hashes, energies, eigenvalues, occupations, all
   SCF/work counters, every verification attempt and completion markers for
   every pair. Keep the full arm's complete traces, including failed checks
   and recovery. No relaxed tolerances or changed stopping criteria. A crashed
   or unverified run is a failed delivery, not a sample to drop; retain its
   wall/output and report the failure before any acceptance decision.

Timing command shape (substitute the established dev08 launcher/binding;
record the exact argv and environment for every sample):

```sh
mkdir -p "$sample_dir"
BAANDI_DAVIDSON_STAGES="$level" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -f '%e %x' -o "$sample_dir/external-time.txt" \
  mpirun -np 8 --bind-to core --map-by core \
  "$binary" run "$input" --output-dir "$sample_dir/output" \
  >"$sample_dir/stdout" 2>"$sample_dir/stderr"
```

Measure the entire launcher, through MPI finalization. Do not redirect to
`/dev/null`, omit density/log/provenance files, time only the SCF call, subtract
verification/recovery/output, or replace production output with a dry run.
Use a monotonic higher-resolution external wrapper if `/usr/bin/time` rounds
away the A/A floor; apply it identically to A/A and A/B. Sum all output files
and stdout/stderr for bytes; keep timing-wrapper metadata separately.
Run `scripts/parse_results.py` after timing to verify completion/provenance;
postprocessing comparisons do not belong inside delivered SCF wall.

Return raw data, hashes, exact commands/exits, A/A statistics, crossed-pair
statistics, trace/results/total bytes and per-case identity evidence. This is
BAANDI full versus BAANDI lightweight, not a QE comparison. Any future solver
default promotion still requires the independent strict oracle and diverse
cell/PP/material generality gates; silicon-only timing cannot authorize it.
