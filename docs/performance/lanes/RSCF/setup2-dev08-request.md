# P3-SETUP2: exact setup reuse, work accounting and dev08 request

Status: implemented, default off, Mac correctness checked; dev08 pending.
Base: `0e5cfa6f03cfeff35aa4cfc849a0345f6265b8d8`, uncommitted task diff.
Derived from the initialization and shell-cache principles in the own-words
`docs/research/qe-study/qe-scf-cycle.md` and `qe-stage-by-stage.md`.
No QE source accessed. No tolerances, band counts, QR guards or starts changed.
The central ledger/log, registry/dashboard and generated README are outside
OWNED_FILES; the master integrates this record after review. No dev08 login,
remote changes, installs or new QE reference run was performed here.

## Diagnosis and implemented arms

The ticket relays Si128 MPI8 rank-zero setup 3.38 s, matched QE init_run
0.83 s and wfcinit 0.68 s. Before REFSTART_MKL, the relayed setup 4.81 s
included atomic QR/completion 1.50, H/Gram 0.95, nonlocal diagonal 0.58,
orbital fill 0.40, local tables/phases 0.41, diagonalize/rotate 0.21,
density/field plans 0.24/0.17 and basis 0.13 s. These motivate the inspection;
they are not independently reproduced QE references or fresh measurements.
The REFSTART_MKL ledger records its setup 4.83 -> 3.36 s screen and adoption.
Keep that and the adopted phase-reuse switch ON in every remote arm.

| Arm | REFSTART_DIAGONAL_REUSE | REFSTART_ORBITAL_REUSE | Removed work |
| --- | --- | --- | --- |
| A | off | off | Current main arithmetic and work |
| D | on | off | Repeated owned-projector binary searches and duplicate ORC partner diagonal |
| O | off | on | Atom-independent exact radial-times-harmonic evaluations/lookups |
| DO | on | on | Both independent changes |

Both new switches default OFF. D keeps every coupling term, including zeros,
and its accumulation order. It caches at most one atom's complex projectors
and the last representative's scalar diagonal. O caches one radial channel's
owned amplitudes after the first atom, preserving multiplication order and
quadrature; one-atom species bypass the cache. No phase factorization or radial
interpolation is introduced. These are additional removals beyond P3-SETUP's
already-adopted direct-phase reuse, not another claim for that same saving.

The full FLOP/byte/collective ledger is in
[reference-scf.md, Setup and output cost](../../../design/reference-scf.md#setup-and-output-cost).
Key remaining costs:

- CholQR2: 8 L J^2 GEMM FLOPs, two 8 J^2-byte Gram reductions; the replicated
  complex scalar factor/inverse is O(J^3). New children separate Gram,
  factor/inverse and rotation. The second pass and fallback remain unchanged;
  this does not repeat A52's rejected change of orthogonalization algorithm.
- MGS fallback: J(J+1) scalar reductions; completion: 2q+1 per trial column.
  The new fallback clock is present only when used. None of the recorded
  valid fixtures entered it; no fallback speedup is claimed.
- H/Gram applies K=max(J,B) columns and reduces K-square projected H in
  B-column chunks. Reducing K to B would change the initial Ritz subspace.
  For Si128 J=512/B=260 *if confirmed by its manifest*, this means 512 H
  columns, not 260; this shape is illustrative until the remote run confirms it.
- The nonlocal diagonal already uses owned G representatives, not all G on
  every rank. Its repeated searches, rather than a missing MPI partition,
  are removed. All-G replicated layout metadata still exists elsewhere.
- Orbital quadratures already have an exact-bit shell cache. O removes the
  repeated map lookups/harmonics across atoms, not the shell integrations.
  Low-symmetry cells still have more distinct shell keys; approximate tables
  are outside the identical-start contract.
- FFT field/operator/density plans are separate layouts with mutable scratch.
  The rustfft planner is local to each construction; same-length axes share
  within that planner, not between runs. Cross-run immutable plan caching
  belongs in the numerics adapter and is outside this ticket's ownership.

## Mac evidence and reproduction

Machine: arm64 macOS 27.0.1 (26A434), cargo 1.99.1. Pure backend, dev profile
(workspace opt-level=1, dependencies opt-level=2), no release or real MPI/MKL
run. Single observations, sequential A/D/O/DO order, no noise floor, no timing
acceptance or extrapolation to dev08. The first large cubic A sample has a
visible cold/load outlier in the diagonal; it is retained, not selected away.
Some build activity overlapped the local observation window.

```sh
cargo test -p dft-engine --test refscf setup_phase_clocks_cover_small_start -- --nocapture
cargo test -p dft-engine --test refscf setup2_ -- --nocapture --test-threads=1
cargo test -p dft-engine --test refscf setup2_larger_start_profile -- --nocapture
```

The first command was executed before the implementation. Original one-atom,
5-band, 2 Ha cubic/skew/elongated fixtures confirm the existing split; the
larger post-instrumentation A arm below isolates the new nested QR/H phases.
The first `setup2_` run contained the two small tests; the large test was
added afterward and run separately. Final `refscf` runs include all three.

Small start identity matrix: cubic, skew, elongated, slab; one and two species;
nonlocal s/s/p projector structure with off-diagonal coupling and local-only
potentials; J>B, J=B and kinetic completion J<B; 1/3/8 simulated ranks,
1/2 workers, existing phase reuse both on and off. Each of three candidate
arms matches its one-worker A baseline on every owned orbital, Ritz value,
static diagonal, density, reciprocal density, ionic and initial potential bit.
This is 72 arm/rank/shape/worker comparisons (plus 12 baseline starts),
not a rank-count bit-identity claim. Full trajectories use the default
threshold and passing independent fixed-point verification, with no relaxed
tolerance. All four arms match every recorded scientific value and count:

| Synthetic two-atom fixture | Energy (Ha) in every arm | Maps | Davidson iterations incl. verification | H columns incl. verification |
| --- | ---: | ---: | ---: | ---: |
| Cubic | -1.16726646059305894 | 21 | 96 | 356 |
| Skew | -1.16960710991236394 | 13 | 47 | 237 |
| Elongated | -1.01375345105461845 | 9 | 43 | 169 |

All starting-state and trajectory energy differences are **0 Ha**, with
bit-identical density/orbitals. These are synthetic regression oracles, not
QE or material validation. Large fixtures below use 16 synthetic atoms,
64 atomic columns, 34 bands, 4 Ha cutoff, 161 radial samples, 80 projector
columns (C2=400), serial pure backend. Both candidates individually and
combined retain every start-state bit across the multi-panel GEMM boundary.

| Large fixture | G=L | Grid points F | QR GEMM FLOPs | H-Gram FLOPs | Final rotation FLOPs | D: baseline -> cached searches | O: baseline -> cached amplitudes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| large-cubic | 1021 | 15625 | 33,456,128 | 8,364,032 | 4,443,392 | 816,800 -> 40,880 | 65,344 -> 4,084 |
| large-skew | 1043 | 15625 | 34,177,024 | 8,544,256 | 4,539,136 | 834,400 -> 41,760 | 66,752 -> 4,172 |
| large-elongated | 937 | 18000 | 30,703,616 | 7,675,904 | 4,077,824 | 749,600 -> 37,520 | 59,968 -> 3,748 |

All three: QR reductions 2 x 32,768 bytes; H-Gram 2 reductions totaling
32,768 bytes; projector sums (34 + 30) x 80 doubles, totaling 40,960 bytes.
The operator reports 4 transpose all-to-all calls and 2 projector tree sums
for the two H chunks at one rank (local copies, not network traffic). No
candidate changes these counts. O scratch peaks at 24 L bytes; D at 80 bytes.

| Fixture | Arm | start total (s) | nonlocal diagonal (s) | orbital fill (s) | QR/completion (s) | H/Gram (s) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| large-cubic | A | 0.041476791 | 0.011602500 | 0.002190916 | 0.001333542 | 0.007987542 |
| large-cubic | D | 0.019405500 | 0.000468375 | 0.001947916 | 0.001199708 | 0.007996750 |
| large-cubic | O | 0.020279917 | 0.004020583 | 0.000751750 | 0.001105292 | 0.007194209 |
| large-cubic | DO | 0.016369333 | 0.000430250 | 0.000740500 | 0.001090667 | 0.007190625 |
| large-skew | A | 0.037922625 | 0.004440792 | 0.003877083 | 0.001134583 | 0.007238291 |
| large-skew | D | 0.033898208 | 0.000456584 | 0.003894333 | 0.001175667 | 0.007247291 |
| large-skew | O | 0.036178208 | 0.004608041 | 0.002115917 | 0.001118292 | 0.007218042 |
| large-skew | DO | 0.032683875 | 0.000463584 | 0.002233083 | 0.001175666 | 0.007653250 |
| large-elongated | A | 0.024216875 | 0.004302333 | 0.002282209 | 0.001212416 | 0.006700916 |
| large-elongated | D | 0.020834958 | 0.000481416 | 0.002223584 | 0.001301042 | 0.007053167 |
| large-elongated | O | 0.024590958 | 0.004700542 | 0.001069375 | 0.001362917 | 0.007240209 |
| large-elongated | DO | 0.020594583 | 0.000470583 | 0.001096917 | 0.001295333 | 0.007270583 |

Complete A-arm split, seconds (nested children are labelled; never add them
to their parents). These `start()` totals exclude the driver density-plan
and mixer setup; the verified-driver fixture above includes those.

| Phase | Large cubic | Large skew | Large elongated |
| --- | ---: | ---: | ---: |
| atomic_block_total | 0.012126291 | 0.012617833 | 0.010617833 |
| atomic_completion | 0.000000000 | 0.000000042 | 0.000000000 |
| atomic_density | 0.001749375 | 0.000677125 | 0.000704792 |
| atomic_diagonalize_rotate | 0.000612375 | 0.000366833 | 0.000420542 |
| atomic_h_apply | 0.007755125 | 0.007018376 | 0.006465708 |
| atomic_h_apply_gram | 0.007987542 | 0.007238291 | 0.006700916 |
| atomic_h_gram | 0.000212958 | 0.000200625 | 0.000209000 |
| atomic_local_tables_structure_factors | 0.007502375 | 0.015048417 | 0.004359500 |
| atomic_orbital_fill | 0.002190916 | 0.003877083 | 0.002282209 |
| atomic_orthonormalize_complete | 0.001333542 | 0.001134583 | 0.001212416 |
| atomic_qr_factor_inverse | 0.000464251 | 0.000388958 | 0.000445624 |
| atomic_qr_gram | 0.000396750 | 0.000340792 | 0.000337958 |
| atomic_qr_rotate | 0.000425625 | 0.000396791 | 0.000418958 |
| basis_gvectors | 0.001305125 | 0.000530750 | 0.000549625 |
| distribution_owned_rows | 0.000241792 | 0.000121583 | 0.000131375 |
| ewald | 0.000584041 | 0.000196000 | 0.000286500 |
| field_fft_plan | 0.002242542 | 0.000796458 | 0.001185000 |
| initial_potential | 0.001166084 | 0.000415125 | 0.000485834 |
| nonlocal_diagonal | 0.011602500 | 0.004440792 | 0.004302333 |
| nonlocal_tables | 0.001338709 | 0.002411750 | 0.000949416 |
| operator_fft_plans_and_tables | 0.001495167 | 0.000637208 | 0.000629750 |
| potential_distribution | 0.000030792 | 0.000026333 | 0.000010125 |
| validate_species | 0.000023958 | 0.000001458 | 0.000001584 |

No `atomic_qr_mgs_fallback` entry occurred. The replicated factor/inverse
is 0.00039 to 0.00046 s in these J=64 A-arm observations, while H dominates the atomic
block. Larger J has a cubic factor/EVD cost; do not infer its dev08 share
from the small Mac fixtures.

Local transient logs (not committed) and SHA-256:

- `target/tmp/setup2-baseline.log`: `415730351806bf359b0c78e34401fecca6793c68acc9f59de350283b51b8ac93`
- `target/tmp/setup2-tests.log`: `4276a23ce1b9c8475c51add7ea9abd2ad7cd376b4e347c63df11fea325ba6eaa`
- `target/tmp/setup2-large.log`: `31d3bba7463c539cbe47fd77953845df2350b193e50ca952249ddac750fbaabc`

## dev08: exactly one profiled full run per arm/case

Master session forwards this request to `baandi-dev08`; this implementer does
not access the server. Use si64, si128, elongated and rskew from the existing
P3 campaign manifests (`/data/sig-dft/p3b/<case>/` or their explicitly frozen
successor recorded by the master), with the identical geometry, UPF, grid,
cutoff, band count, start, standard threshold and independent verification.
All are **MPI8**, one thread per rank. Match QE MPI8 only. Archive QE's raw
init_run/wfcinit lines and the original command/version/input/UPF SHA-256.
Unavailable clocks stay unavailable; do not create a QE reference from the
relayed motivation numbers. Reuse matched QE evidence; this requests no new
QE campaign or additional timing repetitions.

Build the reviewed candidate once with the established dev08 release
`cargo build --release -p dft-cli --features mpi,mkl` and the existing MKLROOT,
MPI/runtime configuration; do not install or reconfigure libraries. Record
build exit code, reviewed source/diff hash, binary SHA-256, rustc/MPI/MKL
versions, input/UPF SHA-256 and complete switch manifests. Use the established
socket-0 eight-physical-core mapping, bind and record each rank; no oversubscription.

Execute A, D, O, DO once each on each of the four cases: **16 BAANDI runs**.
No warmup, repeated pair, A/A, median or statistical speedup campaign is
requested by this ticket. Each invocation is its own <=600 s fast job; retain
failures/timeouts with their elapsed cost rather than silently replacing them.
Do not truncate the SCF/verification tail for a timeout. This is a profiled
development screen, not a resolved-noise/default-adoption claim.

Launch template: the master fills the absolute manifest paths/core list and
sets `diagonal`/`orbital` according to the arm table; the scheduler submits
one invocation at a time. All other switches are frozen to the reviewed
baseline (in particular REFSTART_MKL=on, REFSCF_SETUP_PHASE_REUSE=on).

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on \
BAANDI_REFSTART_MKL=on BAANDI_REFSCF_SETUP_PHASE_REUSE=on \
BAANDI_REFSTART_DIAGONAL_REUSE="$diagonal" \
BAANDI_REFSTART_ORBITAL_REUSE="$orbital" \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c "$cores" \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$baandi_binary" run "$baandi_input" --output-dir "$unique_output"
```

Keep production verification ON and identical across arms. Existing gate/reuse
settings and all unrelated experiments must match the baseline manifest.
Capture all stdout/stderr, external launcher wall, exit code, start/end time,
CPU binding, per-rank RSS, each rank's results.json, full density/eigenvalue
artifacts and every scientific history/verification record. Setup clocks alone
are not total delivered wall. Detailed setup clocks have no new barriers;
rank-zero time contains waiting and is not the critical-rank maximum.

Required adjacent setup_s/QE table (pending means unmeasured):

| Case, MPI8 | A engine setup_s / rank max (s) | D (s) | O (s) | DO (s) | Matched QE init_run (s) | QE wfcinit (s) |
| --- | --- | --- | --- | --- | --- | --- |
| si64 | pending | pending | pending | pending | retrieve matched raw clock | retrieve matched raw clock |
| si128 | pending (ticket motivation: rank 0 3.38) | pending | pending | pending | ticket motivation 0.83; verify provenance | ticket motivation 0.68; verify provenance |
| elongated | pending | pending | pending | pending | retrieve matched raw clock | retrieve matched raw clock |
| rskew | pending | pending | pending | pending | retrieve matched raw clock | retrieve matched raw clock |

For each row also provide CLI-inclusive setup_s, engine-only total, density
plan/mixer setup and main_to_first_map with their boundaries. QE init_run and
wfcinit have different scopes; wfcinit is nested, never added to init_run.
Do not form a ratio between unlike boundaries or between parallel models.

For every arm, list every setup phase and the new QR/H children for rank 0
and min/max across ranks. Show which rank attains each maximum; do not sum
maxima of different ranks or nested children. Include L/J/B/K/Q/C2, unique
shell counts if already exposed, H chunk widths, Gram shapes, 2mnk FLOPs and
local reduction payload bytes; API collective counts and padded FFT exchange
payloads are separate from useful bytes. Confirm whether MGS fallback occurs.
Use the added factor clock to decide if real triangular/library factorization
is the next ticket; use H versus Gram children to decide if H kernels dominate.

Correctness gate: A versus D/O/DO must retain identical start and delivered
scientific values at the fixed rank count, density bytes, map/SCF counts,
Davidson iterations, H columns, restarts/repairs/drops, first redo and all
verification attempts/outcomes. Collect diagnostic initial spectra/residuals
already available; if raw starts are not emitted, run the existing start
identity tests on the MKL build and report that limitation explicitly.
A/D/O/DO rounding divergence is a failure of these exact-reuse claims, not a
reason to loosen the gate. Preserve strict 1e-10 Ha / standard 5e-9 Ha checks
against the campaign's independent reference, without fabricating a new one.
Mac bit equality is insufficient to promote a default on dev08.

## Local gate record

| Executed command | Exit code | Observed result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace test targets type-check; final rerun 0.44 s |
| `cargo test -p dft-engine --test refscf` | 0 | 24 passed, 1 pre-existing ignored; final run 0.36 s |
| `cargo test -p dft-numerics switches` | 0 | 8 library + 6 integration switch tests passed |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings` | 0 | Formatting and all-target lint gates passed |
| `BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh dft-numerics` | 0 | build/test/fmt/clippy PASS; 89 passed, 22 ignored, 38 s wall |
| `git diff --check` | 0 | No whitespace errors |

The check.sh load-only diagnostic could not invoke `ps` in the sandbox
(PermissionError); its four actual verification gates all ran and passed.
No full workspace test suite, release/MKL execution, real MPI or QE run was
performed. No README reference numbers or default promotion is implied.
The three new tests are below the ten-second debug budget (combined small
matrix/trajectory observations 0.20 s, larger start test 0.34 s in their
individual runs). The task is left uncommitted for master review.
