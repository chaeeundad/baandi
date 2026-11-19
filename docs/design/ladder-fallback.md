# G-LADDER-FALLBACK: coarse basis feasibility

2026-11-19. Base commit `a8347ddfbb9df5c32f6947de720266b1282ddf76`.
Mac Apple Silicon, pure Rust release, one thread, no MPI or QE execution.
Default behavior is unchanged. The opt-in change is a planning/diagnostic fix,
not a new SCF algorithm or a measured speed improvement.

## Finding and reproduction discrepancy

The fallback is intentional protection of the full-cutoff calculation, but
**discovering a predictable dimension mismatch by constructing an invalid
coarse problem is a defect in ladder planning and diagnostics**. No design
contract deliberately excludes these named materials. There are too few
coarse plane waves for the requested bands, not too few for an unconditionally
required 48-column subspace.

Eight of the nine reported errors reproduce with the actual checked-in inputs
in `benchmarks/cases/`, reached through `docs/performance/benchmark-suite.md`
and `scripts/bench_suite.py`. The ninth, **syn-odd-shift**, does not reproduce
that error: its checked-in input has a shifted 3x1x1 k mesh and six bands.
It takes the existing planned `multiple k points` fallback before any coarse
basis is constructed. This is true on the pristine binary and the patch.
Removing only the mesh to make a diagnostic Gamma variant actually permits
its ladder (seven coarse PW, six bands). Thus the reported dev08 ninth error
cannot be attributed to this fixture without the exact dev08 input/overrides.
A Gamma variant with eight bands would fail, but that is not the checked-in
input and is not asserted to be the remote input.

### Actual dimensions at the requested half cutoff

Counts below are full-complex PW dimensions. Gamma ORC packs the conjugate
pairs into the same number of real coordinates; it does not halve the
available dimension. B is the requested band count and the Davidson block
size. C is the requested maximum subspace, `max(48, 6 B)`; the usable maximum
is `min(C, Npw)`. Cutoffs are in Hartree (fixture inputs are in Rydberg).

| Fixture | Fine / coarse cutoff (Ha) | Fine / coarse PW | B = block | C / usable coarse capacity | Default decision |
|---|---:|---:|---:|---:|---|
| syn-cubic-degenerate | 4 / 2 | 27 / 7 | 8 | 48 / 7 | construction error, ordinary SCF |
| syn-d3 | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-fft-primes | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-p3 | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-s1 | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-s2 | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-spd | 1.5 / 0.75 | 7 / 1 | 6 | 48 / 1 | construction error, ordinary SCF |
| syn-z1 | 4 / 2 | 27 / 7 | 8 | 48 / 7 | construction error, ordinary SCF |
| syn-odd-shift (original mesh) | 4 / not constructed | first k: 23 / not constructed | 6 | not used | multiple k points |
| syn-gaussian-fcc (positive control) | 4 / 2 | 113 / 27 | 8 | 48 / 27 | ladder runs |
| syn-odd-shift (diagnostic Gamma variant) | 4 / 2 | 27 / 7 | 6 | 48 / 7 | ladder runs |

For the 4-Bohr cube the first nonzero kinetic shell is
`(2 pi / 4)^2 / 2 = 1.23370055013617 Ha`: 0.75 Ha holds only G=0,
whereas 2 Ha holds G=0 and six axial vectors. This independently explains
both rejection dimensions. The coarse density sphere and its FFT grid do
not supply additional wavefunction degrees of freedom.

`ScfProblem::new` invokes `bands_fit_subspace` inside `ordered_solves` while
constructing each k state. Its condition is `B <= Npw && (capacity > B ||
B == Npw)` with capacity-per-band constrained to 2..12. Equality is a valid
full-space solve; neither 48 PW nor two complete band blocks are required.
`ordered_solves` wraps the input error with `Numerical("k index 0: ...")`.
`run_ladder` catches the coarse construction/solve/transfer error and formats
`coarse stage: {error}`, then runs the untouched fine problem. The former
eligibility check covered mode restrictions, not dimensions. The former
`coarse_plane_waves=0` is an unfilled report field, **not a measured zero-PW
basis**. No coarse SCF map or H application occurred in these failures.

## Opt-in implementation

`BAANDI_LADDER_BASIS_CHECK=on` (or `[performance] ladder_basis_check = "on"`)
checks the effective coarse wavefunction basis before `coarse_problem`.
The registry default is **off**. It runs after the existing eligibility and
fraction-rule decisions, uses the same basis constructor and the same
`bands_fit_subspace` predicate, and uses the fine problem's actual band count
(the unchanged electrons/occupation/band options give the same coarse count).
It does not reduce bands, change a Davidson capacity, increase the cutoff, or
alter any stopping tolerance, verification, density or orbital initialization.

An incompatible plan records:

- Engine typed `CutoffLadderReport.coarse_basis_rejection`, containing PW,
  bands, block size, requested/effective subspace and capacity-per-band.
- Existing CLI `scf.switches.BAANDI_LADDER_BASIS_CHECK.path =
  "skipped_coarse_basis"`: a stable machine-readable decision, distinct from
  `compatible_coarse_basis` and `basis_check_unavailable`.
- Existing fallback fields carry, for example,
  `coarse_basis_incompatible: plane_waves=1 bands=6 block_size=6 requested_subspace=48 effective_subspace=1 capacity_per_band=6`.
- Existing ladder `coarse_plane_waves` is filled with the measured count;
  coarse work, time and grid remain zero because no coarse problem ran.

The CLI serializer is outside this ticket's edit scope: the typed rejection
is available to engine callers, while CLI consumers use the existing switch
path/fallback fields and coarse PW count. No JSON or logging enters a numerical
loop. Other basis-construction errors retain the old coarse-error fallback.
Ineligible multi-k/complex paths retain their original reason and do not
consult the new switch. Accepted plans proceed through the original code.

**Which nine would now run a ladder? None.** This is the explicitly permitted
preflight-and-skip fix. Eight become planned basis skips; syn-odd-shift remains
a planned multi-k skip. The gaussian-fcc control continues to run a ladder.
Automatically raising the cutoff is a different trajectory experiment: it
would need shell-boundary selection, nesting/grid checks and fresh accuracy
and work measurements. It is not necessary to fix this diagnostic defect and
is not implemented. Earlier rejected fraction/stop policies in
[ladder-fraction.md](ladder-fraction.md) are not repeated.

## Release measurements and accuracy

Pristine CLI saved before editing; candidate built with the same offline
release command. Three seeds 601, 602, 603; one thread; all other environment
`BAANDI_*` values cleared. Standard uses the production 5e-9 Ha threshold and
fixed-point verification. Strict uses `BAANDI_REFERENCE_STRICT=on`, input
energy tolerance 1e-11 Ha, density tolerance 1e-9, eigen tolerance 1e-10 Ha;
no tolerances were loosened. Inputs otherwise come from `case_input` and
`synthetic` in `scripts/bench_suite.py`, including the original mesh/FFT grid.

Total raw H is fine history H + fine certificate/verification H + fine atomic
start H + coarse history H + coarse atomic start H. Verification is already
inside history certificate records and is counted only once. Grid-weighted
H replaces the two coarse terms by `coarse_h_full_grid_equivalents`.
For rejected ladders coarse H is zero, so raw and weighted totals agree.
Do not read the old ladder-only `total_h_full_grid_equivalents` field as the
whole-run H for a fallback: its fine-stage counters were never populated.

### Original geometries, requested nine

Each number below is the same for seeds 601 / 602 / 603. The deterministic
free-electron start makes these repetitions, not three independent geometries.
All 54 conditions converge. The on, off, unset and no-ladder arms have exactly
the same scientific results: **max |dE| = 0 Ha, density SHA-256 identical,
identical eigenvalues, scientific histories, certificates and H counts**.
This satisfies the paired strict 1e-10 Ha and standard 5e-9 Ha classes.

| Fixture | Standard total H (each seed) | Strict total H (each seed) | On minus no ladder H |
|---|---:|---:|---:|
| syn-cubic-degenerate | 70 | 185 | 0 |
| syn-d3 | 13 | 91 | 0 |
| syn-fft-primes | 11 | 142 | 0 |
| syn-odd-shift | 186 | 311 | 0 |
| syn-p3 | 13 | 99 | 0 |
| syn-s1 | 14 | 138 | 0 |
| syn-s2 | 43 | 161 | 0 |
| syn-spd | 11 | 92 | 0 |
| syn-z1 | 51 | 105 | 0 |

Positive controls, also identical between on/off/unset:

| Control / profile | No ladder raw H | Ladder raw H | Ladder grid-weighted H | Absolute dE vs no ladder (Ha) |
|---|---:|---:|---:|---:|
| gaussian-fcc / standard | 366 | 345 | 250.920 | 1.9587886868066562e-11 |
| gaussian-fcc / strict | 731 | 769 | 618.472 | 9.769962616701378e-15 |
| odd-shift Gamma / standard | 58 | 99 | 89.9296875 | 2.1899815294545988e-10 |
| odd-shift Gamma / strict | 204 | 209 | 184.8125 | 5.329070518200751e-15 |

These control ladder effects predate the fix and are not credited to it.

### Additional independent geometry check

To avoid mistaking seed repetitions for generality, additionally use
`suite.variant(original_case, seed, 0, 0.001)` for the requested nine: independent
fractional displacements in +/-0.001, with the generator's per-case SHA-derived
random streams, seeds 601..603. Same four candidate arms and two pristine arms.
All on/off/unset/no-ladder scientific results still agree exactly, including
the one pre-existing nonconverged strict condition. No default or threshold
was changed to obtain convergence.

| Fixture | Standard H, seeds 601 / 602 / 603 | Strict H, seeds 601 / 602 / 603 | Paired H change |
|---|---|---|---:|
| syn-cubic-degenerate | 70 / 70 / 70 | 185 / 185 / 185 | 0 |
| syn-d3 | 13 / 13 / 13 | 92 / 91 / 91 | 0 |
| syn-fft-primes | 11 / 11 / 11 | 129 / 147 / 142 | 0 |
| syn-odd-shift | 185 / 186 / 186 | 310 / 310 / 311 | 0 |
| syn-p3 | 13 / 13 / 13 | 93 / 93 / 98 | 0 |
| syn-s1 | 15 / 14 / 14 | 113 / 113 / 138 | 0 |
| syn-s2 | 19 / 20 / 22 | 109 / 607* / 134 | 0 |
| syn-spd | 11 / 11 / 11 | 92 / 92 / 98 | 0 |
| syn-z1 | 49 / 49 / 49 | 131 / 115 / 106 | 0 |

`*` syn-s2 seed 602 strict exhausts the unchanged 100-map budget, CLI exit 3,
on the pristine default/no-ladder and all four candidate arms. Its retained
energy is 0.6363255217103463 Ha and all density hashes agree. **53/54 displaced
conditions converge; no accuracy-class claim for the unconverged condition.**

An exploratory extension to displaced gaussian-fcc stopped at pristine
seed 602 standard/no-ladder with `Diverged: density residual exceeded growth
limit` (CLI exit 4; confirmed by a separate repeat). Its scratch parser initially raised KeyError because the
failure JSON has no total_energy; the parser was hardened and the final
additional matrix was scoped to the requested nine. This exploratory failure
is retained under `target/scratch/displaced-baseline/`, not hidden or claimed
as a passing measurement.

Independent absolute accuracy is not implied by equality to no ladder.
For the original nine excluding syn-s2, standard versus the tight strict
no-ladder reference differs by at most 1.4557688388094903e-10 Ha (<5e-9).
Syn-s2 differs by 0.005042673044315826 Ha even in the pristine no-ladder arms:
this is the already documented synthetic fixed-point multiplicity, not a
change introduced here (see [multiplicity](../validation/scf-fixed-point-multiplicity.md)).
It prevents an absolute cross-profile accuracy claim for that fixture;
paired identity remains exact. Displaced converged cases other than syn-s2
have a maximum corresponding difference of 1.4823964278321e-10 Ha.

Across original geometries plus controls: 132 pristine/candidate scientific
identities and 66 on/off/unset identities. Across the displaced nine: 108 and
54 respectively (including identical failed outcomes). In total **240 pristine
identities and 120 on/off/unset identities**. Density files are compared by
SHA-256; all scientific history/verification/energy fields are compared after
removing durations. Registry/provenance/diagnostic fields intentionally differ.

## Reproduction, checks, and scope

All scratch, generated synthetic UPFs and binaries are uncommitted under
`target/scratch`; no UPF is added to the repository. The complete input, stdout,
stderr, density and results files are retained per arm. Scratch runners:
`target/scratch/survey.py` and `target/scratch/analyze.py`.

```sh
export TMPDIR_990="$PWD/target/scratch"
export TMPDIR="$TMPDIR_990"
mkdir -p "$TMPDIR_990"
# Before applying the patch:
cargo build --offline --release -p dft-cli
cp target/release/dft-cli target/scratch/baseline-cli
python3 -I target/scratch/survey.py baseline
# After applying the patch:
cargo build --offline --release -p dft-cli
python3 -I target/scratch/survey.py candidate
python3 -I target/scratch/analyze.py
python3 -I target/scratch/survey.py displaced-nine-baseline
python3 -I target/scratch/survey.py displaced-nine
python3 -I target/scratch/analyze.py displaced-nine
CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli
python3 -I scripts/tests/test_readme_figures.py
```

Each scratch survey runner returns 0 when the survey completes, including
findings; each row records the actual CLI exit code. Underlying execution is
`target/{scratch/baseline-cli,release/dft-cli} run <input.toml> --output-dir <out>`.
Canonical inputs are reproducible directly using the checked-in suite runner,
`--cases syn-cubic-degenerate,syn-d3,syn-fft-primes,syn-odd-shift,syn-p3,syn-s1,syn-s2,syn-spd,syn-z1,syn-gaussian-fcc --seed 601 --threads 1 --repeats 1`,
with new `--work-dir` and `--output` under scratch. Repeat for seeds 602/603;
`--variants 1 --displacement 0.001` generates the additional geometries.
The strict inputs use the three tighter values specified above; arms set
`BAANDI_LADDER_BASIS_CHECK` unset/off/on or `BAANDI_CUTOFF_LADDER=off`.

Binary SHA-256:

- Pristine: `753ab6aedd17247c98d2fe489f6aed074bc786694f364294f950be1d494a7c18`.
- Candidate: `1c45285ddfe0617f6a1f2138f9fc891bfba57885f59fac50ed57bc52ae1ac03f`.

Scratch summary SHA-256 (raw rows include input/output paths and CLI exits):

| File under target/scratch | SHA-256 |
|---|---|
| baseline/summary.json | 269b4fbdad4b6c577ed27d2b751ef4c4e651505f97f43f5e03417015e780c6ab |
| candidate/summary.json | 4aa609953875380316f45e8dda8fe3ba659ac0d16c3ba726531bef02db55c6a7 |
| displaced-nine-baseline/summary.json | 5adb28e65419bcec9ee451c21f00c0afd084d4ee99b4d18e3acafc28a48df96a |
| displaced-nine/summary.json | 7fcdbf13394e26dd60df0d7307f8253da5e2027be20c8643f1a8385bfb663a3b |

Both release builds exit 0. Both final analysis commands exit 0. Required fast
gate exits **1**: build/fmt/clippy pass; **470 tests pass, three fail, 134
ignored**, 51 binaries, 164 s total. The failures are the known CLI pseudo
cache tests rejecting in-repository TMPDIR; the load-only `ps` probe also
reports sandbox EPERM but the gate continues. New tests cover one/seven PW,
bands below/equal/above PW, clipped subspaces smaller than 48, correspondence
with real coarse construction, and strict/standard fallback result identity.
No gate or test was skipped or relaxed. README check exits **1**, eight pass,
one pre-existing structural phrase expectation (`standard profile is
unsupported`) fails; generated-content freshness passes.

No wall-time gain is claimed: these tiny shared-Mac runs overlapped gate
compilation, and H work changes by zero on the requested nine. No dev08, MKL,
MPI, external-potential or QE accuracy/performance evidence was produced.
The broader UPF-dependent `bench_scf.py` suite was not rerun; this ticket's
synthetic release matrix supplies the before/after oracle. README generator,
evidence bundles, dashboard and session registry are outside OWNED_FILES;
no edits made there. The copied registry has no G-LADDER-FALLBACK entry and
this worktree has no `.delegate/G-LADDER-FALLBACK/report.json`; orchestration
state and dashboard publication remain with the master session. Changes are
uncommitted as requested.

Dev08 follow-up through the master: retain the exact survey input and
results.json for syn-odd-shift, compare on/off on the original nine at seeds
601..603 with the same single-core model, verify the structured decisions and
scientific identity, and rerun the fast gate with its permitted external cache
scratch. There is no proposal to promote the switch from these Mac counts.
