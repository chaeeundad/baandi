# P3-LADDER: reference cutoff ladder, dev08 request

Base: `09f6628be2bfb1d788166b4dd3540088fdb5d79f`, plus the uncommitted
P3-LADDER diff. No QE source was accessed. This request is for the master
session and `baandi-dev08`; no remote command was executed by this task.
`REFSCF_LADDER` stays **off by default**. Mac work is correctness and counters,
not performance acceptance.

## Mechanism and prior evidence

The task's relayed si128/si64 map matching and per-map gap motivate removing
full-cutoff maps. C1-DESIGN, C1-LADDER and C1-LADDER-DEV08 establish the prior
legacy experiment, including its fixed fraction 0.5 and sensitivity to
realistic displacements. ALGEBRA-STREAM-B1-MAC concerns a different mechanism
and is not evidence for reference-driver ladder speed. This implementation
uses the same reference driver twice with a fixed coarse/fine estimator ratio
100, coarse budget 30, exact owned prolongation and fresh fine Broyden history.
No Si tuning, S2 seeding, fine stopping change or verification change.
See [the contract](../../../design/reference-scf.md#cutoff-ladder).

## One run per arm, MPI8 only

Use a single candidate binary and compare **D** (`BAANDI_REFSCF_LADDER=off`)
with **L** (`BAANDI_REFSCF_LADDER=on`). Run each arm exactly once per row and
verification setting; no seed sweep, warmup SCF, repeat or parameter scan.
Alternate the D/L execution order by row to limit drift. Keep all other
switches, build features, inputs and backend settings identical; RMM remains
off. Freeze the current verification judgment mode explicitly in both arms.
This is a screen, not a statistical speed claim or authority to promote the
default. A follow-up measurement needs its own master-session work order.

| Geometry | Displacement | Parallel model | Verification settings | Runs |
| --- | --- | --- | --- | ---: |
| si64 | registered unperturbed input | MPI8, one thread per rank | on, off | 4 |
| si128 | registered unperturbed input | MPI8, one thread per rank | on, off | 4 |
| elongated | registered generality input | MPI8, one thread per rank | on, off | 4 |
| rskew | registered generality input | MPI8, one thread per rank | on, off | 4 |
| si32 | fractional uniform +/-0.001, seed 601 | MPI8, one thread per rank | on, off | 4 |
| si32 | fractional uniform +/-0.01, seed 601 | MPI8, one thread per rank | on, off | 4 |
| si64 | fractional uniform +/-0.001, seed 601 | MPI8, one thread per rank | on, off | 4 |
| si64 | fractional uniform +/-0.01, seed 601 | MPI8, one thread per rank | on, off | 4 |

Total: 32 candidate executions. Generate the displaced coordinates once with
the established C1 displacement generator, archive them, and copy those exact
coordinates to D, L and the strict oracle. These are fractional amplitudes,
not Bohr amplitudes. Use the existing matching strict reference where available;
otherwise obtain one strict-oracle run per geometry through the master session
and record that additional work separately (not another D/L repetition).
Do not substitute an oracle with different coordinates, cutoff, bands or PP.

Build with the current approved MPI/backend configuration on dev08. Archive
commit plus diff, binary SHA-256, build command, toolchain/backend/MPI versions,
input files and their hashes, PP SHA-256, complete switch manifests, host load
and binding. Use the same eight physical cores, same rank binding and one
thread per rank in both arms. Do not compare this row with threads-only QE or
BAANDI. No new QE run is requested; if an existing QE comparator is quoted,
it must be MPI8 at exactly matching conditions with full run provenance.

Illustrative invocation (retain the site's approved launcher/binding options):

```sh
BAANDI_REFSCF_LADDER=off OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 --bind-to core /path/to/candidate run D/input.toml --output-dir D/results
BAANDI_REFSCF_LADDER=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -p mpirun -np 8 --bind-to core /path/to/candidate run L/input.toml --output-dir L/results
/path/to/candidate refscf-check --stop-to-verified L/results/results.json
```

Inputs select `[scf] driver="reference"`, the registered fine cutoff/grid,
bands, atomic start, beta/history and threshold; set `verify="fixed_point"` or `"off"`
for the two settings. Use separate output directories. The `time` wall must
include launch, all coarse/fine setup, maps, verification/recovery, output and
MPI process exit. Respect the 600 s fast-job limit; report any timeout rather
than silently reducing accuracy or extending the job.

## Required table and gates

One row per geometry, amplitude and verification mode, with columns:

- D and L complete external wall (s), and raw L/D screen ratio;
- D fine maps; L coarse maps and fine maps; redo/retry/verification call counts;
- L coarse stage seconds (setup + maps + export), fine setup/handoff seconds,
  fine-map seconds and verification seconds; state rank-local/max convention;
- D and L delivered energy (Ha), E_D-E_strict and E_L-E_strict (Ha), and L-D;
- fine estimator stop index, first fine stop-to-verified seconds, verification
  verdict/attempts, coarse threshold and `coarse_stopped` flag;
- completion/exit code, input/PP/binary hashes and exact commands.

`reference_scf.cutoff_ladder` contains the L coarse records and counters;
`reference_scf.iterations` and `convergence.history` remain fine only. Each
record has a stage field. `wall_s.start` contains coarse time, so
`wall_s.coarse_in_start` is a nested subtotal, not an additional contribution.
The stage's `map_start_s` clocks share one driver-entry origin. Verification
cost includes failures and recovery; never drop a failed/slow arm. The run's
energy audit is at the fine cutoff.

Hard checks: all verify-on arms converge and pass the unchanged fine
verification, |E-E_strict| <= 5e-9 Ha, and tight D/L comparisons <= 1e-10 Ha
where both arms use the registered tight threshold. Verify-off arms must
report verification off, never a passing certificate. Check finite energies,
charge conservation and complete records. Report negative speed results and
additional fine maps without tuning fractions or thresholds after seeing them.
The requested one-run screen cannot establish median gain or run-to-run noise.

## Local evidence

Mac synthetic results (not QE references), from
`cargo test -p dft-engine --test refscf -- --nocapture` (exit 0):
42 passed, 5 pre-existing ignored. All four new fixture tests passed.
The matrix runs 48 reference SCFs: four systems, simulated ranks 1/3/8,
tight D/L at 1e-12 Ry with verification, and standard L at 1e-8 Ry with
verification on/off. Four independent strict legacy solves supply the oracle.

| Fixture | Ranks | Tight D fine maps | Tight L coarse / fine maps | Standard L coarse / fine maps (verify on/off) | Max tight abs(L-D), Ha | Max standard abs(L-strict), Ha |
| --- | --- | ---: | --- | --- | ---: | ---: |
| cubic | 1/3/8 | 9 | 9 / 9 | 5 / 6 | 1.394441e-13 | 1.628367e-10 |
| skew | 1/3/8 | 11 | 9 / 9 | 5 / 6 | 1.243450e-14 | 5.042811e-11 |
| elongated | 1/3/8 | 8 | 7 / 8 | 5 / 5 | 2.309264e-14 | 1.012728e-10 |
| synthetic nonlocal | 1, 3, 8 | 26 | 18 / 15, 18 / 16, 18 / 15 | 11 / 8 | 2.553513e-14 | 1.330873e-9 |

All verify-on runs pass the unchanged fine verification. Standard delivered
energies satisfy 5e-9 Ha against the independently certified strict result;
tight D/L agree within 1e-10 Ha. Nonlocal tight counts differ at three ranks,
so rank-count work identity is not claimed. Cross-rank energies remain within
1e-10 Ha; no cross-rank bit-identity claim is made. Two small local-only cases
save no tight fine maps; these counts do not demonstrate a universal gain.

The first fixture attempt failed in the **legacy strict oracle**, before any
nonlocal ladder run: its atomic initial subspace stagnated and triggered the
existing density-growth guard. The oracle now explicitly uses the existing
deterministic random start (`start_wavefunctions=FreeElectron`,
`free_electron_initial_guess=false`) for all four systems. Its tolerances and
certificate are unchanged. The reference D/L atomic starts and ladder rule
were not changed. The initial failing test exit was 101 and is not passing
evidence.

`cargo test -p dft-engine cutoff_ladder` (exit 0) includes the new owned
prolongation test plus unchanged legacy tests: 21 passed, 4 ignored across
the matching targets. The new test exercises ranks 1/3/8, exact copied orbital
bits, bitwise Gamma conjugacy, orthonormality within 1e-12, charge within
1e-13 electrons, and analytic/serial Fourier padding within 1e-15 Bohr^-3.
These are numerical tolerances, not claims that FFT arithmetic is exact.

Additional local gates:

- `cargo check --workspace --tests`: exit 0 on the final source.
- `cargo test -p dft-cli reference`: exit 0; 14 passed, 9 pre-existing ignored,
  including ladder JSON/clock accounting with verify on/off and stage-reader
  compatibility/rejection. Existing fine stop-to-verified tests also run in
  the combined gate below.
- `cargo test -p dft-numerics switches`: exit 0; 14 passed.
- `cargo fmt --check`: exit 0.
- `cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings`:
  exit 0. An earlier large-enum warning was corrected by boxing the coarse
  handoff payload; no numerical change.
- `python3 -I scripts/tests/test_readme_figures.py`: exit 0; 9 passed. Existing
  generated README evidence is fresh. No new README result type, source
  evidence JSON, dashboard/session JSON or ledger row was written outside
  this ticket's OWNED_FILES; this section is the task's local attempt record.
- `BAANDI_TEST_TMP="$PWD/target/tmp" TMPDIR_990="$PWD/target/tmp" CHECK_SCCACHE=0
  scripts/check.sh fast dft-engine dft-cli dft-numerics`: exit 1; build, fmt and clippy pass;
  tests report **570 passed, 3 failed, 152 ignored** across 58 binaries.
  All three failures are existing pseudopotential-cache fixtures rejecting
  the sandbox's in-checkout scratch root: `fetch_boundary_verifies_before_publication_and_never_overwrites`,
  `cache_resolution_and_no_network_side_effects`, and
  `cached_fetch_checks_hash_and_provenance`. They are not passing evidence.
  The sandbox also denies the script's advisory `ps`/`nice` operations.
  The task-specific commands above pass without weakening cache protections.

Reproduction commands:

```sh
cargo check --workspace --tests
cargo test -p dft-engine --test refscf
cargo test -p dft-engine cutoff_ladder
cargo test -p dft-cli reference
cargo test -p dft-numerics switches
cargo fmt --check
cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings
```
