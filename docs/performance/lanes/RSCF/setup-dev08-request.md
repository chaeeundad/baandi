# P3-SETUP: setup and delivery accounting, dev08 request

Status: local implementation; remote measurement pending through the master
session and baandi-dev08. No QE source was accessed. Only the own-words study
notes `qe-scf-cycle.md` and `qe-stage-by-stage.md` informed the phase boundaries.
No Mac observation establishes a speedup or default promotion.

## Freeze the experiment

Reuse the exact Si64 and Si128 MPI8 inputs/UPFs/bands/grid/geometry from
`/data/sig-dft/vc/<case>/setup.md` and its manifests (verify-cost A/B). Reuse the
existing matched QE runs, with their command/version/input/UPF SHA-256, MPI8,
start, cutoff, grid, bands, convergence threshold and CPU binding recorded.
Do not substitute a fresh geometry or compare with QE OpenMP. The ticket's
relayed main 2aee083 medians (three samples) motivate the work only:

| case, MPI8 | BAANDI to first map (s) | matched QE init_run (s) | QE wfcinit (s) |
| --- | ---: | ---: | ---: |
| Si64 | 1.091 | 0.190 | not supplied |
| Si128 | 4.865 | 0.750 | 0.680 |

These are user-relayed values, not a new QE reference. Retrieve the original
raw timing lines before using them in the final table. QE init_run excludes
some process/input costs that BAANDI's main-to-first-map includes; show the
engine and CLI splits side by side, with these scope differences, not a ratio
of unlike intervals. Si128 output before serialization 0.772 s, unaccounted
0.624 s and total map-energy evaluation 1.295 s are likewise relayed observations.

Record reviewed main/candidate source tree and binary SHA-256, toolchain,
MPI/MKL versions, build commands/features, input/UPF hashes and complete switch
manifests. Use the existing release MPI/MKL build and physical core allocation,
8 ranks, one thread per rank, explicit binding, OMP_NUM_THREADS=1 and
MKL_NUM_THREADS=1. Master runs `scripts/check.sh full` and the release correctness
matrix on dev08; Mac checks deliberately exclude release/full-suite work.
Keep REFSCF_VERIFY_GATE=on, REFSCF_VERIFY_REUSE=off and
REFSCF_VERIFY_GATE_FIRST=off in both arms, with fixed-point verification on.

## Measurement sequence

1. One profiled full-convergence run per arm/case, with
   BAANDI_DAVIDSON_STAGES=on and BAANDI_PHASE_PROFILE=on. Preserve root and
   every `rank*/results.json` file.
   Read setup_s and output_s, not the historical init_lane zeros. Capture
   external launcher wall and per-rank RSS. Verify schema, positive used spans,
   scope/closure, and successful independent verification before timing.
2. A/A noise floor: 20 labelled pairs total per case, bounded batches in one
   measurement window. Then 30 crossed A/B pairs per switch/case, AB/BA
   alternating, with recurring A/A controls. Detailed profiles OFF;
   the new coarse setup/output clocks remain ON in every arm. One case/model
   per job, maximum 600 s including warmups/output/analysis; Si128 starts with
   one pair per job. Failed/censored runs remain in the delivered-cost record.
3. Report same-window medians for both arms, paired log(B/A), SE, A/A floor F,
   and all raw samples. Resolved means |mean|/SE >=3 and |mean| >=3F. Any
   reproducible gain is eligible; no percentage cutoff. Inconclusive is not a
   rejection or promotion. Cold and warm cache runs are separate strata.
4. Widen survivors to the existing rskew and elongated verification-cost
   fixtures and synthetic non-Si/multi-projector shapes before promotion.
   Preserve strict 1e-10 Ha / standard 5e-9 Ha classes, passing independent
   verification, thread-count consistency and exact science for overhead-only
   changes. Compare each rank count only against itself for bit identity.

## Required table and closure

For each arm/case show per-rank raw seconds and min/max (not sum-of-maxima):
bootstrap; CLI parse/hash/preflight/remaining; engine validation/basis/Ewald;
field FFT; ownership; complex projectors; operator FFT/real tables; nonlocal
diagonal; local/atomic tables and phases; atomic-density transforms; initial
potential/distribution; atomic fill/QR/completion/H-Gram/EVD-rotation;
density plan/mixer preparation; first-map boundary. Put the matched QE
init_run and wfcinit lines in adjacent columns and state unavailable subtimers.
Collectives are nested; atomic_block_total is the parent of its four children.
A rank maximum for each phase need not come from the same rank and cannot be
summed into a critical-path total.

Follow with output density/orbital gather, result assembly, CLI encode/hash,
writes, report collectives, binary identity, JSON body serialization/write and
remaining external wall. The JSON trailer/close/MPI-exit gap is explicitly
outside measured body-write cost. Keep per-map energy_s as a separate row;
QE has no matching separate audit timer, so leave its counterpart unavailable.
Preserve input/output hashes, eigenvalues, every scientific history/verification
record, map/H counts and density.bin bytes. Default versus baseline and each
cache A/B must be identical excluding clocks/provenance/configuration fields.

Dashboard/session registry, README generated evidence and central
ledger/optimization-log updates belong to the master session outside this
worktree's OWNED_FILES. This request records the local attempt for integration.

## Implemented candidate and explicit arms

Only one new switch, appended to the registry to avoid P3-DAV-DENSE overlap:

| arm | BAANDI_REFSCF_SETUP_PHASE_REUSE | purpose |
| --- | --- | --- |
| A | off | identical existing direct phase calculations |
| B | on | cache direct atom/G phase across charge/ionic and orbital m channels |

The default remains **off**. No previously rejected phase-factorization,
radial interpolation, QR change, eigensolver change or binary hash cache is
reintroduced. The latter cache is already adopted; keep its configuration
identical in both arms. No change is made to the separately measured map-energy
evaluation. The raw phase is reused without changing arithmetic or summation
order. Optional scratch <=16*max(owned rows, atoms) bytes/rank is included in
the conservative peak. This is not a full atom-by-G replicated cache.

Mechanism: atomic local tables currently evaluate each atom/G phase twice
(charge then ionic); B evaluates it once. Atomic orbital fill currently repeats
that phase for each of 2*l+1 magnetic channels; B evaluates once per atom/G and
radial channel. Radial transforms/shell lookup, raw orbital entries, QR, H,
small solves, gathers and SCF work counts must stay identical. Inspect
`atomic_local_tables_structure_factors` and `atomic_orbital_fill`, plus total
setup and delivered wall, for the effect. Wider timing of radial quadrature
may still dominate these combined spans; do not infer its cost from this A/B.

Reproduction skeleton (the master supplies exact recorded paths/build):

```sh
BAANDI_REFSCF_SETUP_PHASE_REUSE=off BAANDI_DAVIDSON_STAGES=off \
BAANDI_PHASE_PROFILE=off BAANDI_REFSCF_VERIFY_GATE=on \
BAANDI_REFSCF_VERIFY_REUSE=off BAANDI_REFSCF_VERIFY_GATE_FIRST=off \
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --bind-to core --report-bindings "$CANDIDATE" run "$MATCHED_INPUT" --output-dir "$ARM_OUTPUT"
```

Use B with only SETUP_PHASE_REUSE=on changed. Preserve `rank1/results.json`
through `rank7/results.json` under the arm directory, plus root results.json.
Do not change profile flags between acceptance arms. Startup records are
always present; detailed profiling is only the separate diagnostic pair.

## Local observation and correctness record

Mac ARM64, Cargo dev/test profile (repository opt-level settings, no release
build), one small synthetic atom, s/p atomic channels, 161 radial samples,
2 Ha cutoff, five bands, one map. First instrumented pre-optimization run:

| shape | engine setup (s) | local/atomic tables (s) | atomic fill (s) | atomic block total (s) |
| --- | ---: | ---: | ---: | ---: |
| cubic | 0.001819500 | 0.000274625 | 0.000046333 | 0.000235375 |
| skew | 0.002751917 | 0.001682709 | 0.000176083 | 0.000238958 |
| elongated | 0.001139334 | 0.000294000 | 0.000065583 | 0.000125334 |

Single samples, observation only: local tables dominate the skew fixture;
cubic field FFT planning is 0.000279250 s (similar to local tables). This
justifies examining repeated phase work, not claiming a measured speed gain.
Reproduce: `cargo test -p dft-engine --test refscf setup_phase_clocks -- --nocapture`.
Raw first output: `/private/tmp/p3-observe.log` (local scratch, not portable
committed evidence). No QE measurement was run.

`setup_phase_reuse_preserves_owned_start_bits` compares 32 per-rank off/on
pairs: four cell shapes (cubic/skew/elongated/slab), local/nonlocal synthetic
potentials, two atoms, nine bands including kinetic completion, simulated
rank counts 1/3 and worker counts 1/2. Every orbital, Ritz value, density,
reciprocal density, ionic/potential value and static diagonal is bit-identical.
All detailed profile switches are off there, while setup clocks stay present.
A separate fully verified synthetic fixture checks identical complete science
and histories for switch off/on at both diagnostic levels. The initial new
full-convergence test used a four-electron nonlocal fixture that did not
converge even in baseline; replaced it with the existing two-electron local
verification fixture. No tolerance was loosened and no production default
was changed. The nonlocal fixture remains in the start-state equality matrix.

Four actual CLI Si2 runs (30 Ry, eight bands, one thread) also pass: phase
reuse off/on crossed with lightweight/full diagnostics (full additionally
PHASE_PROFILE=on; STAGE_TIMERS=off in all four). All exit 0, six maps, verified,
E = -7.295530467525706 Ha, maximum off/on energy difference **0 Ha**, identical
energy components, eigenvalues, scientific iteration/verification records and
density bytes. Density SHA-256:
`6cdcf83e464a410d201ab05930bcb96c5a1f1d6b39e0e5ed8cdfe47fca14e81e`.
This also exercises the existing launch-profile append after the new JSON
trailer and parses the resulting file independently with Python's json module.
No optional volumetric output was requested; its writer span is still present.

Local reproduction/evidence: `target/tmp/p3-setup/{si2.toml,observe.py,observations.json}`;
`python3 target/tmp/p3-setup/observe.py` exits 0. The external UPF stays outside
version control. Provenance:

- input SHA-256: `ab8f15c571dde59d11f8ca8fa51c2afd2d9fd52233383ac98410fb7b9e12be72`
- UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`
- observed debug binary SHA-256: `37e4287d21b973e82c9fde14834014303f764a0b338ded7c9857a4dc76a78205`

First cold off/light sample: engine setup 0.017509250 s, local/atomic tables
0.004278125 s, atomic block 0.002647875 s; CLI before JSON 0.003962625 s,
JSON body serialization 0.000216250 s, body write 0.000142083 s. Other samples
are retained in observations.json; cold/warm and profiler differences forbid
using these four runs as a speed comparison. No timing-based default decision.

## Local gates and integration boundary

Base HEAD `07389881f37ea1d42811884ec33306c16a91d1d7`; changes uncommitted.
Commands were run on Mac, with their actual outcomes:

| command | exit | result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | workspace/tests type check |
| `cargo test -p dft-engine --test refscf` | 0 | 19 passed, 1 ignored |
| `cargo test -p dft-cli reference` | 0 | 9 passed, 9 ignored across harnesses |
| `cargo test -p dft-numerics switches` | 101 | 13 passed, 1 failed: registry count assertion |
| `cargo fmt --check` | 0 | clean |
| `cargo clippy -p dft-engine -p dft-cli -p dft-numerics --all-targets -- -D warnings` | 0 | clean |
| `python3 target/tmp/p3-setup/observe.py` | 0 | four verified CLI runs, exact science/density |
| `git diff --check` | 0 | clean |

The numerical switch failure is introduced by the authorized appended switch:
`crates/dft-numerics/tests/all/switches.rs:301` asserts `TABLE.len() == 136`,
where the new table contains 137. That test file is outside this ticket's
OWNED_FILES; it has not been modified or bypassed. Integration requires the
scope owner to update the expected count and assert the new default is off.
All other switch checks pass. A narrow ownership extension was requested;
until granted, this remains an explicit failed gate, not a passing check.

Earlier exploratory engine test exit 101 was the unsuitable nonconverging
synthetic fixture described above; the corrected existing convergent fixture
passes the final suite. No production tolerance or gate was weakened.
`scripts/check.sh full`, real MPI/MKL, dev08 A/A and A/B were not run here,
as required by the task's Mac-only validation limit. No dashboard/README or
central registry file outside ownership was edited; the local registry has no
P3-SETUP entry and no `.delegate/P3-SETUP/report.json` at inspection time.
