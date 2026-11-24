# P3-SETUP3: generalized atomic start, accounting and dev08 request

Base: `d5937bd33a209e8a1c978d17d2303bb56d85413e`, uncommitted task diff.
Both new switches default OFF. No dev08 execution, new QE reference, speed
claim, source access, dependency change or default promotion is made here.
Derived from the starting-subspace principles in the own-words notes
`docs/research/qe-study/qe-scf-cycle.md` and `qe-stage-by-stage.md`.

## Diagnosis and implementation

Ticket motivation, relayed at fa055919: Si128 MPI8 setup 2.75 s versus QE
init_run 0.76 s; Si64 0.74 versus 0.19 s, one run each. Si128 children include
H 0.64 s, H/Gram 0.73 s, QR/completion 0.53 s, QR 0.51 s, local tables/phases
0.39 s, density plan 0.29 s, field plan 0.17 s and basis 0.13 s. H is nested
inside H/Gram and QR inside QR/completion: these numbers are not additive.
They motivate inspection, not a fresh or independently reproduced QE comparison.

Full dimensions, FLOPs, reduction payloads and fallback contract are in
[Setup and output cost](../../../design/reference-scf.md#p3-setup3-generalized-atomic-ritz-start-and-compact-g-lookup).

- `REFSTART_GENERALIZED=on`: for J>=B form S=X^t X and A=X^t H X on raw
  L-by-J atomic vectors. Reuse the existing scaled guarded Cholesky inverse
  T, check T^t S T, solve the small T^t A T eigenproblem with the library
  EVD and return X(T U_B). No wide QR passes. H still acts on J columns,
  not B; truncating the raw input to B would change the Ritz space.
- J<B keeps legacy kinetic completion. Bad pivots or metric defect above
  1e-12 keep legacy CholQR2/MGS. No rank truncation or changed tolerance.
  Candidate is rounding class; compare subspaces rather than orbital signs.
- `REFSTART_COMPACT_G=on`: omit the duplicate density-grid map and owned
  BTreeMap. Use the public Miller mapping and binary search of sorted owned
  FFT positions. Numerical coefficients and their accumulation order are
  unchanged. Both arms already evaluate radial/phase data on used owned G
  only; this removes metadata, not previously nonexistent global value work.
- The generalized path saves one J-square Gram reduction and one factor;
  its net dense-GEMM saving is 6 L J^2 - 8 J^3 - 2 J^2 B FLOPs per rank.
  Small-rank rows can cost more. The explicit metric check is included in
  this accounting and all timing; it must not be silently excluded.

## Scope remaining outside OWNED_FILES

The task's FFT reuse and broader basis-construction removal cannot be safely
implemented in `refstart.rs` alone. `refscf/density.rs::DensityPlan::new`
constructs another density basis and full-grid DistributedPairFft after
start(). Its FFT and StartFields' FFT are private mutable owned objects;
neither has a constructor accepting transferred/shared plans or metadata.
The operator also owns a distinct sparse wavefunction FFT. Caching these by
dimensions alone is incorrect (rank, rank count, stick-policy switches,
backend, axis choice and active G indices matter). No global cache or
cross-run retained memory was introduced to bypass the API boundary.

Concrete next ownership set: `crates/dft-engine/src/refscf/density.rs` and
`driver.rs` to accept existing density metadata and a run-owned full-grid
transform; `crates/dft-solvers/src/distributed_fft.rs` and
`crates/dft-numerics/src/fused_pair_plan.rs` for immutable FFT lane-plan
sharing with per-consumer mutable buffers; their tests. A basis subset or
shared enumeration API additionally belongs to `crates/dft-basis/src/lib.rs`.
This ticket only removes the redundant density FftGridSpec (F reverse entries,
D forward entries) and S_r-entry lookup tree. Remaining plan/basis timers must
remain visible; no saving is claimed for them.

The central attempts ledger, optimization log, sessions/dashboard artifacts,
README generator/evidence and scripts are outside this task's OWNED_FILES.
This document is the integration record for the master. The local registry
contains no P3-SETUP3 entry and `.delegate/P3-SETUP3/report.json` is absent;
no new delegation/session process was launched from inside the delegated task.

## Mac correctness evidence

Pure Rust, default dev test profile; simulated collectives, not real MPI.
No Mac timing is used as performance evidence. Test data are synthetic,
never labelled QE reference. Four added tests cover:

- 36 shape/band/rank configurations (cubic/skew/elongated/slab, 1/2 species,
  local-only and coupled s/s/p projectors, J>B/J=B/J<B, simulated ranks
  1/3/8), each with baseline and three candidate arms at 1/2 workers:
  **216 candidate comparisons**, 252 starts. All non-orbital fields and
  metadata match exactly; C-only and completion outputs match all bits.
  Every generalized start passes 1e-12 Ha Ritz-eigenvalue, 1e-12
  orthonormality and 1e-12 linear subspace-residual checks. The latter is
  max_j ||B_j-A(A^t B_j)||_2, not a squared-angle test hiding O(sqrt(eps))
  errors. G/GC eigenvalues/orbitals are bit-identical across worker counts
  at fixed rank count; no cross-rank-count bit-identity claim.
- Nine fully converged shape/rank fixtures, each with four arms (36 runs):
  default threshold, passing independent verification, identical per-solve
  iterations/H columns/restarts/repairs/drops/reduced orders/active bands,
  per-map redo flags and verification-call ranges. C-only matches the full
  scientific trajectory bits. G and GC are identical to each other.
- Seven wave rows on eight ranks with empty owners; near-coincident centers
  exercise guarded fallback and reproduce legacy start bits; exactly
  coincident centers fail on both paths without silent rank truncation.
- Three 16-atom J=64/B=34 starts cross the 256-row GEMM panel boundary;
  compare baseline with G at 1/2 workers (nine starts), same 1e-12 bounds.

| Large synthetic shape | G=L | Ritz error (Ha) | Orthonormality max error | Linear subspace residual |
| --- | ---: | ---: | ---: | ---: |
| Cubic | 1021 | 1.687538997e-14 | 2.220446049e-15 | 5.142770367e-14 |
| Skew | 1043 | 2.042810365e-14 | 3.885780586e-15 | 4.400575486e-14 |
| Elongated | 937 | 3.952393968e-14 | 3.108624469e-15 | 3.576025570e-13 |

The maximum reported Ritz error is 3.952393968e-14 Ha, orthonormality error
3.885780586e-15, and linear subspace residual 3.576025570e-13. These are
pure-backend synthetic measurements, not silicon or dev08 guarantees.

G/GC versus A (C-only deltas are exactly zero). Counts include verification;
within each shape they match at all three tested rank counts:

| Shape | Simulated ranks | Maps / Davidson / H columns | Signed final delta E (Ha) | First-map dr2, G (Ry) | Signed delta dr2 (Ry) |
| --- | ---: | --- | ---: | ---: | ---: |
| Cubic | 1 | 21 / 96 / 356 | 2.442490654e-15 | 0.827133825166816439 | 5.491163080e-13 |
| Cubic | 3 | 21 / 96 / 356 | 1.239008895e-13 | 0.827133825166772696 | -2.186029135e-13 |
| Cubic | 8 | 21 / 96 / 356 | -7.016609516e-14 | 0.827133825166254333 | -4.847233726e-13 |
| Skew | 1 | 13 / 47 / 237 | 0 | 0.817974135888375198 | 3.393951786e-13 |
| Skew | 3 | 13 / 47 / 237 | -1.332267630e-15 | 0.817974135888305698 | -4.898303985e-13 |
| Skew | 8 | 13 / 47 / 237 | -2.220446049e-16 | 0.817974135889144138 | 1.387112647e-12 |
| Elongated | 1 | 9 / 43 / 169 | 1.332267630e-15 | 0.840466698536150081 | 4.996003611e-15 |
| Elongated | 3 | 9 / 43 / 169 | 1.554312234e-15 | 0.840466698536145085 | 9.103828802e-15 |
| Elongated | 8 | 9 / 43 / 169 | 0 | 0.840466698536140644 | -1.065814104e-14 |

Maximum final |delta E| = 1.239008895e-13 Ha. First-map dr2 changes are
reported rather than asserted equal: this is the loose first Davidson map,
not the final energy. An initial development assertion incorrectly applied
an invented 1e-12 Ry equality bound to this diagnostic; skew MPI8 exceeded
it (1.387112647e-12 Ry) while its final delta E was -2.220446049e-16 Ha
and every count matched. That unsupported diagnostic assertion was removed;
no start-state, final-energy, convergence, verification or count gate changed.
The failed development run is not hidden or relabelled as passing.

Reproduction:

```sh
cargo test -p dft-engine --test refscf setup3_ -- --nocapture --test-threads=1
cargo check --workspace --tests
cargo test -p dft-engine --test refscf
cargo test -p dft-numerics switches
cargo fmt --check
cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings
```

## dev08: ONE profiled full run per arm/case

The master sends this request to `baandi-dev08`. Use the existing P3i/P3j
frozen manifests for **si64, si128, elongated, rskew; MPI8**, one worker,
OpenMP thread and MKL thread per rank, physical socket-0 cores 0-7. Preserve
geometry, coordinates, UPF bytes, cutoff, FFT grid, B, atomic deterministic
start, LDA, standard threshold, verification and all unrelated switches.
Record the exact manifest paths/hashes rather than assuming a prior case
name identifies its bytes. No implementer login, installs or server changes.

| Arm | REFSTART_GENERALIZED | REFSTART_COMPACT_G |
| --- | --- | --- |
| A | off | off |
| G | on | off |
| C | off | on |
| GC | on | on |

Exactly **16 BAANDI invocations** (four arms x four cases), each its own
<=600 s fast job, e.g. `fast-p3-setup3-si128-G`. Full convergence and
independent verification are mandatory. Keep failures/timeouts and elapsed
cost; no silent retry, truncation, warmup, repeated pair, A/A or median
campaign is requested. If a case cannot finish in that limit, report it to
the master rather than change the run or launch an unrequested long job.
This is a diagnostic screen only, not statistical/default-adoption evidence.

Build the reviewed source once in its own target using the established
MPI/MKL installation and build flags:

```sh
cargo build --release -p dft-cli --features mpi,mkl
```

Record build command/exit, source/diff hash, binary SHA-256, compiler/MPI/MKL
versions and complete effective switch manifest. The existing REFSTART_MKL,
REFSCF_SETUP_PHASE_REUSE, REFSTART_DIAGONAL_REUSE and REFSTART_ORBITAL_REUSE
stay ON in all arms. Keep all other adopted switches identical. New switches
are off unless selected by the arm; input values must not mask the intended
settings (check effective manifest). Do not turn off production verification.

Launch template, one invocation per scheduler submission (master fills the
frozen absolute paths and sets generalized/compact from the arm table):

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
BAANDI_REFERENCE_SCF=on \
BAANDI_REFSTART_MKL=on BAANDI_REFSCF_SETUP_PHASE_REUSE=on \
BAANDI_REFSTART_DIAGONAL_REUSE=on BAANDI_REFSTART_ORBITAL_REUSE=on \
BAANDI_REFSTART_GENERALIZED="$generalized" \
BAANDI_REFSTART_COMPACT_G="$compact" \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_STAGE_TIMERS=on \
  /usr/bin/time -p taskset -c 0-7 \
  mpirun -np 8 --map-by core --bind-to core --report-bindings \
  "$baandi_binary" run "$baandi_input" --output-dir "$unique_output"
```

Retain full stdout/stderr, binding/load, start/end timestamps, exit status,
external delivered wall, per-rank results.json, scientific records, density
and eigenvalue artifacts, verification attempts and per-rank RSS. Record
first-map dr2 in Ry and the signed delta against A. Raw initial vectors are
not currently emitted by the CLI: run the new start-oracle tests on the
reviewed MKL build if needed, and explicitly report this observation limit.
Do not mistake absence of a start dump for bit identity.

Required adjacent clock table below is a **setup attribution diagnostic**,
not a speed-comparison table: profiled BAANDI and QE timers have different
boundaries/output costs, so no ratios or total-performance conclusions.
Reuse existing matched QE **MPI8** init_run/wfcinit raw output only; archive
its command, QE version/build, input and UPF SHA-256. Confirm it is the
fastest eligible matched configuration, including diagonalization layout.
If provenance or a clock is missing, write unavailable. No new QE run is
requested, and ticket motivation values are not substituted for raw output.

| Case (MPI8) | A setup_s rank 0 / max (s) | G (s) | C (s) | GC (s) | Matched QE init_run (s) | QE wfcinit, nested (s) |
| --- | --- | --- | --- | --- | --- | --- |
| si64 | pending | pending | pending | pending | retrieve raw clock | retrieve raw clock |
| si128 | pending | pending | pending | pending | retrieve raw clock | retrieve raw clock |
| elongated | pending | pending | pending | pending | retrieve raw clock | retrieve raw clock |
| rskew | pending | pending | pending | pending | retrieve raw clock | retrieve raw clock |

Beside every row include the full setup_s phase split for A/G/C/GC with
rank zero and min/max/argmax. Include engine-only setup, CLI-inclusive setup,
main_to_first_map, density_fft_plan, field_fft_plan, basis_gvectors,
atomic_local_tables_structure_factors, orbital fill, atomic_block_total,
QR Gram/factor/rotation/completion/fallback, overlap Gram/factor/metric,
H/Gram and its H/Gram children, and generalized diagonalize/rotate.
Missing QR clocks in a successful G run mean omitted work, not missing data.
Overlap clocks followed by QR clocks mean guarded fallback. Never sum nested
children or maxima on different ranks, or add wfcinit to init_run.

List actual G, each L_r, J, B, K, grid shape, projector count, chunk widths,
Gram dimensions, dense 2mnk FLOPs and reduction payload bytes. Check the
committed accounting against those shapes. For nominal Si64 J=256/B=132
and Si128 J=512/B=260, confirm rather than infer these numbers from names.
All atomic H columns remain; reduced H time is not automatically less H work.

Correctness: C must be bit-identical to A. G/GC must preserve start spectra,
orthogonality and selected subspaces within 1e-12, and unique-solution final
energy within 1e-12 Ha as an internal regression check. Retain strict
1e-10 Ha / standard 5e-9 Ha accuracy against the campaign's independent
reference and passing verification. Compare every SCF/map/redo, Davidson
iteration, H column, restart/repair/drop, reduced order and verification
attempt/outcome; report any changed count as failure of this ticket's
same-trajectory criterion, not a timing win. Report all first-map dr2 deltas.
No default promotion or tolerance relaxation follows from this one-run screen.

## Local gate record and handoff

| Executed command | Final exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | All workspace test targets type-check |
| `cargo test -p dft-engine --test refscf` | 0 | 28 passed, 1 pre-existing ignored |
| `cargo test -p dft-numerics switches` | 0 | 8 unit + 6 integration tests passed |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings` | 0 | Format and all-target warning-free lint gates |
| `cargo test -p dft-engine --test refscf setup3_ -- --nocapture --test-threads=1` | 0 | Four new tests passed; combined execution 1.10 s, below the debug budget |
| `BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 scripts/check.sh dft-numerics` | 0 | build/test/fmt/clippy PASS; 89 passed, 22 ignored |
| `git diff --check` | 0 | No whitespace errors |

The check.sh load diagnostic could not invoke sandbox-restricted `ps`; its
four actual gates ran and passed. Development failures before these final
gates: an extracted local dimension variable was initially missing
(`cargo check -p dft-engine --tests`, 101); a temporary test-edit delimiter
error (`cargo fmt` parse failure; test compile, 101); the unsupported first-dr2 assertion
noted above (test, 101); and clippy's collapsible-if lint (101, corrected).
No failed command is treated as a passing gate. No numerical acceptance
threshold or production setting was loosened to address these failures.

Local transient logs (not committed), SHA-256:

- `target/tmp-setup3-evidence.log`: `fe65610a6755cd6298d3b730a2dd07b8af01c30f5b51550935fd5f6213687563`
- `target/tmp-setup3-check.log`: `2c328ec6b2054f772d75a8056a7a62f919ebce2bbb45f00bafe3c9e0c1820924`
- `target/tmp-setup3-refscf.log`: `4282e2464b0404327e458534ac3f8d1b4234cd830235e06d25a3a90c3c67d1c2`
- `target/tmp-setup3-switches.log`: `166a2fceed1ce741d46ee1a12c8d0d635503c8efe1b2353860ad7d2ce1d9023d`
- `target/tmp-setup3-clippy.log`: `cf3de609f76223213ac30ed0ffaa5e22fbd410120178e2760be8b9e32bc3281d`
- `target/tmp-setup3-gate.log`: `468d2c759fabd1b36fa0398db0b8715e5037cba5469e3fd7d755b3606b3960a1`

The primary accounting/generalized-start/fixture-count/request criteria are
satisfied locally. Task item 3 is only partially implemented: compact G/grid
metadata removal is tested, but FFT-plan sharing and broader basis reuse
remain blocked by the ownership boundary documented above. No claimed speed
saving, real MPI/MKL validation, or default promotion is implied. The master
can review this uncommitted diff and assign that remaining API work separately.
