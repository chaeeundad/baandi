# P36B k-point hybrid budget

Starting commit: `4bce6acbe563fb7aec04dd3db876765b05b0110f`.
Implementation and local correctness gates complete after architect scope approval; no performance acceptance claim.

## Premise and frozen gates

Reopen A50 only because THR-4 in thread-width-analysis.md section 6 supplies
new evidence: a full-pool inner budget provides one Hamiltonian scratch slot per
worker, eliminating the modulo collisions of the rejected smaller-slot prototype.
Use the existing persistent pool, fixed arithmetic tiles, and ascending band/k
reductions. No spin policy, tolerance, default, real-path or phase overlap change.
The only duplicated-budget search match is ScfProblem construction in scf.rs;
session.rs and collinear.rs reuse it and need no budget edits.

The ticket's complete bit-identity and timing gates are frozen before changes.
Evidence lives in `P36B-evidence/`, outside the sole `target/` directory.
Baseline binary is saved before changing source. No artifacts files are edited.

## Implementation and ownership

- `ScfProblem::construct` gives each k Hamiltonian and band-density plan the full
  pool budget. `ordered_solves` still limits outer width to min(T,nk); all inner
  tasks use that same pool, not T new threads per k. The existing Hamiltonian
  implementation supplies per-worker scratch and fixed band jobs. No numerical
  kernel, tile shape, or density-reduction ordering was changed.
- `PreparedSymmetry` is immutable, with private geometry, mesh, operations and
  orbit-map fields. Config caches an Arc after successful geometry verification.
  Its lookup checks cell, atoms/species, full mesh, time reversal and tolerance;
  mutable configurations used by relaxation cannot reuse a stale preparation.
  Crystal density operations reuse the initial detection too. The engine checks
  geometry, irreducible points and FFT-grid compatibility before using the map.
- Validation, preparation, preflight, metadata and derivative adapter lookups
  share the config preparation. The ordinary SCF run passes it to the engine.
  Legacy engine `enable_symmetry` remains valid and computes its own preparation.
- CLI threading metadata now describes the actual shared pool: outer=min(T,nk),
  inner=T, threads_used=T. The product outer*inner is not a thread count.
  The old metadata expectation test was updated, retaining its complete
  scientific-output equality checks and adding nondivisible and larger budgets.
- No edits to session.rs or the concurrent real-path work. Budget selection remains centralized in scf.rs.

### Approved additive forwarding (architect decision)

The architect extended ownership to the spin and finite-temperature adapters.
The additional production files are exactly:

- `crates/dft-engine/src/collinear.rs`: new prepared-map forwarding method for
  both channels, preserving rejection of time-reversal reduction. The existing
  full-mesh method now prepares once and delegates, with its signature unchanged.
- `crates/dft-engine/src/finite_temperature.rs`: new forwarding method to the
  contained ScfProblem; the existing full-mesh method remains available.
- `crates/dft-cli/src/spin_run.rs`: forwards the Config-owned preparation.
- `crates/dft-cli/src/finite_run.rs`: forwards the same preparation in both
  spin and unpolarized FD dispatches.

No existing public signature, scientific option, checkpoint field or format
version changed. The new methods are additive. The CLI's four SCF dispatches
now share one preparation from validation through engine setup. No further
production file beyond these four was needed. `crates/dft-cli/tests/spin.rs`
adds a release-only baseline/thread matrix; the engine symmetry test adds a
regression for rejecting antiunitary preparation in collinear models.

The new CLI matrix covers fixed spin, unpolarized FD, common-mu spin FD and
fixed-moment spin FD (0.3 electrons), symmetry off/auto, and all nine budgets.
It compares all scientific JSON (round-trip float serialization) and binary
density with the original main executable. An initial invalid FD fixture with
fully polarized integer targets was correctly rejected; it was replaced by the
existing supported FD fixtures, with the failure log retained. No tolerance or
assertion was relaxed.

## Frozen identity protocol and raw evidence

`P36B-evidence/identity.py` preserves input, stdout/stderr, results.json and
binary density for each arm. The main matrix is 3 meshes (k222/k222s/k444s) x
2 profiles (strict/standard) x 3 symmetry choices (off/auto/auto+crystal) x
9 budgets (1,2,3,4,6,8,10,12,16), with gamma_real off. `audit.py` compares
IEEE-754 encodings of every energy/eigenvalue and the entire scientific SCF
history (including certificate and per-k work counts), excluding only timings.
Density SHA-256 is compared against T=1 and against unchanged HEAD at the same T.
The budget-only arm isolates the scheduling mechanism from symmetry reuse.

`hybrid_budget.rs` reuses the existing odd-electron LDA-PW spin fixture and
synthetic N=1 Fermi-Dirac metal (theta=0.01 Ha, 12 bands, LDA-PZ).
Snapshots include density, eigenvalues, eigenvectors, final solve counts,
certificates and scalar history. The source-only P14A accounting build also
compares aggregate counts of every solve, including certification. Baseline
and candidate are built sequentially in the single target directory; the
baseline script restores candidate source in a finally block. No installed
binary or numerical source is overwritten outside the task worktree.

Observed so far: all 486 main-matrix records agree bit for bit (162 in each of
three arms); all 18 spin/FD snapshots match unchanged HEAD. Across all budgets,
spin has 26 solves / 52 Davidson iterations / 189 H-vector applications;
FD has 12 solves / 26 Davidson iterations / 434 H-vector applications.
No arithmetic-dependent-budget case was found in this matrix.

The first input-harness attempt used Boolean gamma_real=false and exited 1:
the CLI correctly requires the string "off". The corrected matrix uses "off".
The first workspace check failed the old metadata expectations and clippy's
collapsible-if lint. Both were fixed without weakening a numerical assertion.
The original logs are retained alongside the final rerun.

## Timing protocol and limits

Uninstrumented release binaries, ordinary strict symmetric cases k222/k222s/
k444s, T=1,2,4,8,12,16. Each condition has one discarded warmup and five
retained samples per arm, alternating before/budget/shared ordering. Inputs,
outputs, load averages, individual walls and phase times are retained.
No task-owned build runs during timing. This is a shared heterogeneous Mac;
concurrent tasks and the 2 Super + 4 Performance + 6 Efficiency core mix prevent
an identical-core Amdahl or dev08 acceptance claim. Five samples cannot establish
small 1–3% effects. No default/tolerance/QE comparison is changed.

Phase values are enclosing wall times, not summed worker CPU times. Setup,
potential, density, total-energy and mixing are separately recorded. The
reported CLI-and-unattributed remainder is wall minus setup and all recorded
engine phase times; it includes process startup/output and unclassified glue.
These are outside-diagonalization costs, not proof that every operation in the
remainder is serial (some setup/audit work is already pooled).

## Completion gates after architect authorization

`python3 -I P36B-evidence/complete_gates.py` exits 0. Every command and actual
exit code is recorded in `complete-commands.json`. It runs:

1. The source-only accounting build of the release engine `hybrid_budget` test:
   all 18 new snapshots byte-match the original baseline, including total
   Davidson/H counts. No budget-dependent arithmetic was found.
2. Release CLI `spin::shared_symmetry_budget_matrix`: 144 runs, four model
   policies x off/auto x nine budgets x original-main/candidate. All scientific
   JSON and density bytes agree; emitted density SHA-256 agrees within each
   model/symmetry condition. The ordinary API/time-reversal rejection test also
   passes in the workspace gate.
3. Uninstrumented release CLI build, then `scripts/check.sh`: build, tests,
   formatting and clippy all pass with exit 0. Existing public signatures and
   checkpoint formats remain unchanged.
4. README regeneration and nine freshness/specification tests: exit 0, no
   generated README/figure changes.
5. Final 162-condition ticket matrix and five `bench_scf --require-identical`
   invocations covering 25 case/budget comparisons: all exit 0, energy delta
   0 Ha and identical density/eigenvalues/main and certificate work counts.
6. Complete matrix audit: 810 stored records across original baseline,
   budget-only, initial shared-map, previous final and completed candidate.
   All scientific fields match bit for bit; timings alone are excluded.

The initial invalid spin-FD test exited 101 because an integer fully polarized
sector has no interior finite-T target. This was a fixture input error, not a
solver failure or bit mismatch. The retained log is
`spin-fd-sharing-invalid-fixture.log`; the supported common-mu and 0.3-electron
fixed-moment fixtures both pass without changed numerical tolerances.

## Verification commands

- `cargo build --release -p dft-cli`: baseline, budget-only and final builds exit 0.
- `python3 -I P36B-evidence/identity.py before` and `budget`: exit 0;
  162 conditions in each arm. `shared` (initial cache implementation) also exits 0.
- `python3 -I P36B-evidence/audit.py`: initial three-arm audit exit 0,
  486 records with zero differing scientific fields. Final-binary rerun follows.
- `P14A_WORK_ACCOUNTING=1 python3 -I P36B-evidence/fixture_baseline.py`:
  exit 0; both cargo release tests exit 0 and all 18 snapshots byte-identical.
- `python3 -I P36B-evidence/bench_matrix.py`: exit 0, all ten commands exit 0
  (five budget-specific baseline records and five same-HEAD `--require-identical`
  comparisons, each covering gamma/k222s/k222/k444s/si8).
- `python3 -I P36B-evidence/bench_final.py`: exit 0, final-binary comparisons
  at T=1,2,4,8,12 all exit 0; energy delta 0 Ha, density and work unchanged.
- Initial `scripts/check.sh`: exit 1 (metadata expectation and collapsible-if).
  Final `scripts/check.sh`: exit 0; build/test/fmt/clippy all PASS. Original
  numerical assertions retained; metadata tests now check the shared-pool contract.
- `python3 -I scripts/make_readme_figures.py`: exit 0, generated files unchanged.
  `python3 -I scripts/tests/test_readme_figures.py`: exit 0, nine tests pass.
  No new scientific reference was added to tests/reference or README.
- `git diff --check`: exit 0. No new unsafe, dependencies, artifact writes or
  real-path edits. Changes remain uncommitted.

The final benchmark energies (Ha, unchanged HEAD and candidate) are Gamma
-7.29553046821076, k222s -7.9194234124690865, k222 -7.831230057121135,
k444s -7.9259761266480915, Si8 -31.32180360103127; SCF counts 10/10/10/10/15.
These are engine equivalence results, not new QE reference measurements.

## dev08 request

After review, compare unchanged HEAD and the reviewed candidate
on dev08 at equal core count and accuracy, with the ticket's nine-budget identity
matrix first. Then alternate baseline/candidate on symmetric k222/k222s/k444s at
T=1,2,4,8,12,16 with at least five samples (30 pairs for a 1–3% claim), CPU affinity,
node load, binary/UPF hashes, pure/MKL builds, and no build overlap. Preserve full
results, density SHA-256 and main/certificate work counts. Profile solver width
and the setup/audit remainder with at least 30 samples per claimed site; compare
memory usage from per-k T-slot scratch. The Mac observations cannot decide
acceptance or parity with QE's best configuration. No remote changes were made.



## Performance observations and spread qualification

Each cell is a median of five retained runs after one warmup. Wall cells show **median [min, max] seconds**. Spread is **100 × (max − min) / median**. A comparison row is **UNUSABLE** for a claim if any of its three arms exceeds 10%. Passing that filter alone does not establish a speedup, especially for 1–3% differences.

### Initial shared-host batch

Raw samples: `timing-raw.json`. All 324 outputs, including warmups, agree with unchanged HEAD on scientific histories, energies, eigenvalues and density SHA-256.

| Case | T | Before wall (s) | Budget only wall (s) | Shared map wall (s) | Spread before/budget/shared | Claim eligibility |
|---|---:|---|---|---|---|---|
| k222 | 1 | 0.342854 [0.316846, 0.454206] | 0.352146 [0.316812, 0.490827] | 0.369320 [0.315516, 0.488070] | 40.06% / 49.42% / 46.72% | UNUSABLE: >10% |
| k222 | 2 | 0.251303 [0.217885, 0.315268] | 0.217163 [0.204124, 0.453971] | 0.203124 [0.200279, 0.275261] | 38.75% / 115.05% / 36.91% | UNUSABLE: >10% |
| k222 | 4 | 0.157666 [0.155327, 0.311206] | 0.141955 [0.141142, 0.159783] | 0.149650 [0.138401, 0.162416] | 98.87% / 13.13% / 16.05% | UNUSABLE: >10% |
| k222 | 8 | 0.133779 [0.131763, 0.138656] | 0.127406 [0.126084, 0.128521] | 0.124757 [0.123895, 0.125541] | 5.15% / 1.91% / 1.32% | ≤10%; observation only |
| k222 | 12 | 0.125355 [0.124981, 0.127149] | 0.124864 [0.124384, 0.125492] | 0.122800 [0.122045, 0.123198] | 1.73% / 0.89% / 0.94% | ≤10%; observation only |
| k222 | 16 | 0.126474 [0.125795, 0.134371] | 0.126341 [0.125954, 0.128123] | 0.124408 [0.123921, 0.125399] | 6.78% / 1.72% / 1.19% | ≤10%; observation only |
| k222s | 1 | 0.188163 [0.187673, 0.189914] | 0.187458 [0.187200, 0.187555] | 0.185868 [0.184811, 0.186471] | 1.19% / 0.19% / 0.89% | ≤10%; observation only |
| k222s | 2 | 0.121805 [0.121342, 0.122841] | 0.121731 [0.121219, 0.126766] | 0.120367 [0.119611, 0.120662] | 1.23% / 4.56% / 0.87% | ≤10%; observation only |
| k222s | 4 | 0.105619 [0.105320, 0.107442] | 0.103127 [0.102972, 0.104703] | 0.101780 [0.101056, 0.102029] | 2.01% / 1.68% / 0.96% | ≤10%; observation only |
| k222s | 8 | 0.096435 [0.095839, 0.096879] | 0.094253 [0.093624, 0.094919] | 0.092798 [0.091841, 0.093481] | 1.08% / 1.37% / 1.77% | ≤10%; observation only |
| k222s | 12 | 0.096319 [0.095816, 0.100556] | 0.093206 [0.092305, 0.093467] | 0.091444 [0.090822, 0.092308] | 4.92% / 1.25% / 1.62% | ≤10%; observation only |
| k222s | 16 | 0.097329 [0.095539, 0.103991] | 0.095305 [0.094417, 0.108369] | 0.095565 [0.093462, 0.102122] | 8.68% / 14.64% / 9.06% | UNUSABLE: >10% |
| k444s | 1 | 0.757458 [0.756637, 1.220644] | 0.798176 [0.757462, 1.277644] | 0.774613 [0.744801, 1.286300] | 61.26% / 65.17% / 69.91% | UNUSABLE: >10% |
| k444s | 2 | 0.998898 [0.675427, 1.228947] | 1.117206 [0.731111, 1.354361] | 0.988572 [0.733032, 1.405982] | 55.41% / 55.79% / 68.07% | UNUSABLE: >10% |
| k444s | 4 | 0.384534 [0.377103, 0.449293] | 0.361699 [0.356788, 0.413620] | 0.354459 [0.346658, 0.380142] | 18.77% / 15.71% / 9.45% | UNUSABLE: >10% |
| k444s | 8 | 0.313165 [0.296755, 0.490736] | 0.290331 [0.270960, 0.315729] | 0.270459 [0.261337, 0.412483] | 61.94% / 15.42% / 55.88% | UNUSABLE: >10% |
| k444s | 12 | 0.240057 [0.234851, 0.242236] | 0.233817 [0.229084, 0.239902] | 0.222047 [0.215532, 0.231938] | 3.08% / 4.63% / 7.39% | ≤10%; observation only |
| k444s | 16 | 0.242049 [0.238360, 0.250117] | 0.247617 [0.238251, 0.323473] | 0.234434 [0.217558, 0.303328] | 4.86% / 34.42% / 36.59% | UNUSABLE: >10% |

One-minute load during runs: 12.27–18.26. Rows exceeding 10%: 9/18: `[['k222', 1], ['k222', 2], ['k222', 4], ['k222s', 16], ['k444s', 1], ['k444s', 2], ['k444s', 4], ['k444s', 8], ['k444s', 16]]`.

### Quieter repeat after scope completion

Raw samples: `timing-repeat-raw.json`. All 324 outputs, including warmups, agree with unchanged HEAD on scientific histories, energies, eigenvalues and density SHA-256.

| Case | T | Before wall (s) | Budget only wall (s) | Shared map wall (s) | Spread before/budget/shared | Claim eligibility |
|---|---:|---|---|---|---|---|
| k222 | 1 | 0.314850 [0.292347, 0.364174] | 0.307078 [0.291234, 0.338203] | 0.322219 [0.290534, 0.329729] | 22.81% / 15.30% / 12.16% | UNUSABLE: >10% |
| k222 | 2 | 0.244708 [0.219946, 0.245785] | 0.226287 [0.212392, 0.242465] | 0.221880 [0.218479, 0.250118] | 10.56% / 13.29% / 14.26% | UNUSABLE: >10% |
| k222 | 4 | 0.166874 [0.160494, 0.183485] | 0.161428 [0.151288, 0.164824] | 0.162773 [0.144066, 0.167786] | 13.78% / 8.38% / 14.57% | UNUSABLE: >10% |
| k222 | 8 | 0.159746 [0.150311, 0.283580] | 0.150915 [0.138541, 0.258696] | 0.147228 [0.141497, 0.156275] | 83.43% / 79.62% / 10.04% | UNUSABLE: >10% |
| k222 | 12 | 0.157025 [0.155023, 0.247337] | 0.157319 [0.155730, 0.190421] | 0.154953 [0.148088, 0.281045] | 58.79% / 22.05% / 85.81% | UNUSABLE: >10% |
| k222 | 16 | 0.133183 [0.127356, 0.271193] | 0.138486 [0.126451, 0.270475] | 0.158104 [0.126018, 0.259837] | 108.00% / 104.00% / 84.64% | UNUSABLE: >10% |
| k222s | 1 | 0.195731 [0.189261, 0.205046] | 0.190555 [0.188597, 0.193709] | 0.188200 [0.186795, 0.188563] | 8.06% / 2.68% / 0.94% | ≤10%; observation only |
| k222s | 2 | 0.136697 [0.135561, 0.138258] | 0.133985 [0.132692, 0.138631] | 0.132786 [0.129591, 0.137188] | 1.97% / 4.43% / 5.72% | ≤10%; observation only |
| k222s | 4 | 0.112763 [0.112387, 0.122025] | 0.109363 [0.103595, 0.109729] | 0.106903 [0.101888, 0.107872] | 8.55% / 5.61% / 5.60% | ≤10%; observation only |
| k222s | 8 | 0.102203 [0.095807, 0.103175] | 0.102278 [0.098321, 0.108517] | 0.101463 [0.100429, 0.110246] | 7.21% / 9.97% / 9.68% | ≤10%; observation only |
| k222s | 12 | 0.104451 [0.098520, 0.127877] | 0.101165 [0.100173, 0.115289] | 0.101693 [0.097414, 0.121976] | 28.11% / 14.94% / 24.15% | UNUSABLE: >10% |
| k222s | 16 | 0.101469 [0.101413, 0.102775] | 0.102097 [0.100501, 0.127610] | 0.109746 [0.100987, 0.118065] | 1.34% / 26.55% / 15.56% | UNUSABLE: >10% |
| k444s | 1 | 0.783033 [0.767557, 0.799934] | 0.768921 [0.766165, 0.786810] | 0.766117 [0.763066, 0.794649] | 4.13% / 2.68% / 4.12% | ≤10%; observation only |
| k444s | 2 | 0.496168 [0.477690, 0.611138] | 0.484767 [0.470392, 0.973340] | 0.470387 [0.460950, 0.519046] | 26.90% / 103.75% / 12.35% | UNUSABLE: >10% |
| k444s | 4 | 0.352344 [0.346514, 0.382892] | 0.345200 [0.317656, 0.357020] | 0.322762 [0.306956, 0.355809] | 10.32% / 11.40% / 15.14% | UNUSABLE: >10% |
| k444s | 8 | 0.313631 [0.294928, 0.320271] | 0.278012 [0.265403, 0.278791] | 0.265787 [0.262681, 0.272319] | 8.08% / 4.82% / 3.63% | ≤10%; observation only |
| k444s | 12 | 0.247967 [0.235959, 0.253492] | 0.232834 [0.209251, 0.256059] | 0.233831 [0.223573, 0.236840] | 7.07% / 20.10% / 5.67% | UNUSABLE: >10% |
| k444s | 16 | 0.246208 [0.235074, 0.248467] | 0.244255 [0.212675, 0.284287] | 0.229925 [0.225228, 0.247967] | 5.44% / 29.32% / 9.89% | UNUSABLE: >10% |

One-minute load during runs: 6.89–8.74. Rows exceeding 10%: 12/18: `[['k222', 1], ['k222', 2], ['k222', 4], ['k222', 8], ['k222', 12], ['k222', 16], ['k222s', 12], ['k222s', 16], ['k444s', 2], ['k444s', 4], ['k444s', 12], ['k444s', 16]]`.

All 18 repeat batches checked load before admission and waited for the one-minute average to fall below 8, below the initial batch's minimum 12.27. Admitted load range: 6.89–7.98; raw checks: `repeat-load-checks.json`. Every previously flagged condition was repeated; no original sample was discarded.

### Repeat outside-solve phases (seconds, median [min, max])

The remainder includes CLI startup/output and unclassified glue. These are enclosing phase spans, not a proof that the entire remainder executes serially.

| Case | T | Arm | Setup | Potential | Density | Energy | Mixing | CLI + unattributed |
|---|---:|---|---|---|---|---|---|---|
| k222 | 1 | before | 0.025508 [0.023604, 0.028626] | 0.023639 [0.020979, 0.027578] | 0.028108 [0.026935, 0.032577] | 0.011784 [0.011292, 0.013015] | 0.005869 [0.005369, 0.006746] | 0.006432 [0.005556, 0.007525] |
| k222 | 1 | complete | 0.022095 [0.021627, 0.024994] | 0.024125 [0.022169, 0.024373] | 0.030327 [0.026818, 0.031031] | 0.011940 [0.010672, 0.012659] | 0.006033 [0.005660, 0.006172] | 0.005707 [0.005669, 0.005956] |
| k222 | 2 | before | 0.025210 [0.022708, 0.026003] | 0.022931 [0.020846, 0.023546] | 0.023751 [0.021209, 0.024937] | 0.010840 [0.009575, 0.012123] | 0.005805 [0.005229, 0.005992] | 0.006364 [0.005483, 0.006542] |
| k222 | 2 | complete | 0.023814 [0.022658, 0.029982] | 0.022961 [0.022685, 0.027552] | 0.024092 [0.023103, 0.025895] | 0.010796 [0.010139, 0.012318] | 0.005873 [0.005501, 0.006443] | 0.005717 [0.005640, 0.007165] |
| k222 | 4 | before | 0.022662 [0.022059, 0.025498] | 0.022334 [0.020494, 0.026263] | 0.018853 [0.018372, 0.019885] | 0.009334 [0.008740, 0.011675] | 0.005547 [0.005162, 0.006303] | 0.005978 [0.005345, 0.006813] |
| k222 | 4 | complete | 0.023045 [0.020461, 0.028318] | 0.024645 [0.020512, 0.026033] | 0.020025 [0.017856, 0.020523] | 0.010766 [0.008658, 0.011613] | 0.006031 [0.005203, 0.006597] | 0.006014 [0.005106, 0.006384] |
| k222 | 8 | before | 0.025058 [0.023240, 0.039036] | 0.024899 [0.023573, 0.033645] | 0.019661 [0.018799, 0.035777] | 0.011373 [0.009704, 0.021088] | 0.006143 [0.005728, 0.009671] | 0.006859 [0.006054, 0.008014] |
| k222 | 8 | complete | 0.022503 [0.022327, 0.024386] | 0.024512 [0.022914, 0.026844] | 0.018350 [0.018059, 0.019495] | 0.011127 [0.010769, 0.012156] | 0.005998 [0.005664, 0.006543] | 0.006620 [0.006079, 0.007271] |
| k222 | 12 | before | 0.025508 [0.024573, 0.032516] | 0.026590 [0.025197, 0.030761] | 0.019677 [0.019457, 0.036150] | 0.012011 [0.011231, 0.022230] | 0.006813 [0.006020, 0.007811] | 0.006339 [0.005831, 0.008094] |
| k222 | 12 | complete | 0.024091 [0.022798, 0.036295] | 0.026195 [0.025125, 0.053553] | 0.019739 [0.018660, 0.025915] | 0.012001 [0.011266, 0.022930] | 0.006392 [0.006059, 0.007943] | 0.007012 [0.006845, 0.007729] |
| k222 | 16 | before | 0.024821 [0.022373, 0.042449] | 0.021294 [0.020780, 0.041805] | 0.017796 [0.016930, 0.032516] | 0.009021 [0.008793, 0.017271] | 0.005292 [0.005207, 0.008349] | 0.005800 [0.005632, 0.010082] |
| k222 | 16 | complete | 0.031378 [0.021149, 0.039617] | 0.024178 [0.020783, 0.036515] | 0.018641 [0.016732, 0.032299] | 0.009831 [0.008700, 0.026905] | 0.005574 [0.005177, 0.008394] | 0.007695 [0.005263, 0.011369] |
| k222s | 1 | before | 0.013560 [0.013145, 0.017127] | 0.015597 [0.015043, 0.016072] | 0.015170 [0.015056, 0.015722] | 0.009015 [0.008941, 0.009376] | 0.005221 [0.005172, 0.005508] | 0.005019 [0.004858, 0.005989] |
| k222s | 1 | complete | 0.011795 [0.011702, 0.011854] | 0.015104 [0.015028, 0.015257] | 0.015095 [0.015013, 0.015140] | 0.008930 [0.008893, 0.008958] | 0.005193 [0.005166, 0.005217] | 0.004890 [0.004651, 0.006094] |
| k222s | 2 | before | 0.012511 [0.012319, 0.013035] | 0.015957 [0.015426, 0.016333] | 0.013013 [0.012668, 0.013285] | 0.008771 [0.008745, 0.009252] | 0.005329 [0.005269, 0.005521] | 0.005117 [0.004866, 0.005190] |
| k222s | 2 | complete | 0.010719 [0.010492, 0.011161] | 0.015561 [0.015268, 0.015822] | 0.012807 [0.012474, 0.013204] | 0.008822 [0.008632, 0.008974] | 0.005213 [0.005177, 0.005223] | 0.004807 [0.004554, 0.005584] |
| k222s | 4 | before | 0.012352 [0.012020, 0.014932] | 0.015407 [0.015188, 0.017509] | 0.012196 [0.012103, 0.013012] | 0.008705 [0.008610, 0.009894] | 0.005208 [0.005138, 0.005366] | 0.005191 [0.004862, 0.005683] |
| k222s | 4 | complete | 0.010945 [0.010520, 0.011052] | 0.015219 [0.014917, 0.015547] | 0.010372 [0.010077, 0.010645] | 0.008622 [0.008319, 0.008783] | 0.005188 [0.005110, 0.005213] | 0.005062 [0.004653, 0.007166] |
| k222s | 8 | before | 0.012580 [0.012028, 0.012687] | 0.015306 [0.015054, 0.015502] | 0.010438 [0.010123, 0.010704] | 0.008788 [0.008443, 0.008855] | 0.005206 [0.005100, 0.005262] | 0.005138 [0.005001, 0.007495] |
| k222s | 8 | complete | 0.011616 [0.011503, 0.013263] | 0.015511 [0.015352, 0.018228] | 0.010431 [0.010319, 0.011040] | 0.009115 [0.008826, 0.010056] | 0.005283 [0.005243, 0.005773] | 0.004931 [0.004789, 0.005426] |
| k222s | 12 | before | 0.012894 [0.012723, 0.017909] | 0.015549 [0.015221, 0.020078] | 0.011187 [0.010819, 0.013689] | 0.008824 [0.008434, 0.011979] | 0.005326 [0.005226, 0.006246] | 0.005329 [0.005033, 0.005710] |
| k222s | 12 | complete | 0.012136 [0.011080, 0.031994] | 0.015861 [0.015448, 0.017336] | 0.010747 [0.010586, 0.011201] | 0.009007 [0.008799, 0.009515] | 0.005376 [0.005237, 0.005795] | 0.005200 [0.004907, 0.006021] |
| k222s | 16 | before | 0.013084 [0.012820, 0.013342] | 0.015718 [0.015366, 0.015792] | 0.011095 [0.010936, 0.011389] | 0.008898 [0.008794, 0.009231] | 0.005352 [0.005270, 0.005477] | 0.005231 [0.005159, 0.005275] |
| k222s | 16 | complete | 0.013284 [0.011625, 0.013636] | 0.016715 [0.015394, 0.018370] | 0.011523 [0.011074, 0.012503] | 0.009239 [0.008842, 0.010700] | 0.005488 [0.005242, 0.006138] | 0.005257 [0.005107, 0.005646] |
| k444s | 1 | before | 0.040618 [0.039994, 0.042688] | 0.016608 [0.016038, 0.017196] | 0.051651 [0.050080, 0.053510] | 0.015069 [0.014454, 0.015425] | 0.005253 [0.005165, 0.005642] | 0.006305 [0.006259, 0.006606] |
| k444s | 1 | complete | 0.032962 [0.031639, 0.033721] | 0.016242 [0.016132, 0.016764] | 0.050697 [0.050285, 0.053589] | 0.014629 [0.014586, 0.015318] | 0.005239 [0.005174, 0.005457] | 0.005310 [0.005098, 0.005596] |
| k444s | 2 | before | 0.033987 [0.032673, 0.042605] | 0.017358 [0.016030, 0.028324] | 0.037142 [0.035318, 0.047710] | 0.013618 [0.012734, 0.020661] | 0.005767 [0.005376, 0.006425] | 0.006797 [0.006545, 0.009759] |
| k444s | 2 | complete | 0.025625 [0.023307, 0.029647] | 0.016405 [0.015755, 0.019492] | 0.036957 [0.035495, 0.038700] | 0.013204 [0.012468, 0.014730] | 0.005451 [0.005235, 0.005937] | 0.006254 [0.005915, 0.014269] |
| k444s | 4 | before | 0.028523 [0.027911, 0.032812] | 0.017247 [0.016817, 0.022084] | 0.029158 [0.028228, 0.032637] | 0.011898 [0.011407, 0.014495] | 0.005990 [0.005511, 0.007435] | 0.007863 [0.006934, 0.008495] |
| k444s | 4 | complete | 0.020086 [0.018389, 0.022224] | 0.019146 [0.015865, 0.022948] | 0.030175 [0.024722, 0.034917] | 0.013688 [0.010812, 0.015609] | 0.006590 [0.005287, 0.006792] | 0.006389 [0.006104, 0.007655] |
| k444s | 8 | before | 0.029098 [0.025151, 0.030520] | 0.019397 [0.019213, 0.023092] | 0.025324 [0.025231, 0.026098] | 0.013018 [0.012410, 0.014045] | 0.006146 [0.005921, 0.006439] | 0.007489 [0.007349, 0.008730] |
| k444s | 8 | complete | 0.020198 [0.019662, 0.020993] | 0.023146 [0.023039, 0.025794] | 0.024526 [0.024053, 0.025211] | 0.014531 [0.014236, 0.015044] | 0.006866 [0.006737, 0.008122] | 0.007900 [0.007414, 0.009944] |
| k444s | 12 | before | 0.027540 [0.025789, 0.031483] | 0.023089 [0.020643, 0.025954] | 0.022776 [0.021411, 0.027169] | 0.014627 [0.012316, 0.015360] | 0.006617 [0.006242, 0.007505] | 0.008371 [0.007739, 0.008957] |
| k444s | 12 | complete | 0.018699 [0.017772, 0.021028] | 0.021641 [0.020803, 0.022542] | 0.023180 [0.022728, 0.023354] | 0.013730 [0.012591, 0.014200] | 0.006708 [0.006618, 0.007190] | 0.009311 [0.006905, 0.010181] |
| k444s | 16 | before | 0.028910 [0.025091, 0.031312] | 0.021801 [0.020681, 0.022288] | 0.022885 [0.022291, 0.023198] | 0.013826 [0.013480, 0.014271] | 0.006381 [0.006187, 0.007034] | 0.009093 [0.007110, 0.010486] |
| k444s | 16 | complete | 0.019638 [0.017690, 0.020766] | 0.023098 [0.022198, 0.025317] | 0.022918 [0.022555, 0.023141] | 0.014160 [0.012964, 0.014657] | 0.006751 [0.006649, 0.007206] | 0.009046 [0.007617, 0.009580] |

## Final interpretation, reproducibility and cleanup

All correctness gates and the requested measurement protocol are complete;
whole-task status is **DONE**, with no general performance acceptance claim.
The original nine noisy conditions were all repeated. The lower-load batch
still has 12/18 rows above 10% in at least one arm; those rows are explicitly
UNUSABLE for a claim. A lower admission load does not guarantee isolation during
a short run, so no outlier is silently dropped and no median from those rows is
promoted as a gain. The remaining six rows are observations, not dev08 acceptance.

For example, the usable k444s T=8 row records before/budget/complete medians
0.313631/0.278012/0.265787 s, ranges [0.294928,0.320271] /
[0.265403,0.278791] / [0.262681,0.272319] s, spreads 8.08/4.82/3.63%.
The final median is 15.25% shorter on this Mac observation. In contrast the T=16
three-arm row remains unusable because the budget-only spread is 29.32%, even
though the final median is lower. No claim of general scaling, equal-core QE
parity, or dev08 speedup follows. T-slot-per-k scratch increases memory and
should be profiled on the remote acceptance host.

The completed source changes after ownership expansion are listed above;
no other production files were needed. Existing API signatures, checkpoints,
SCF/Davidson arithmetic, reduction order, defaults and tolerances are preserved.
No `unsafe`, new dependencies, real-path edits, artifacts changes or git metadata
writes were introduced. Changes are uncommitted on the requested task branch.

Final evidence: 810 audited ordinary-SCF records (162 conditions per build arm),
18 instrumented before/after engine snapshots, 144 CLI spin/FD runs, 25 final
benchmark comparisons, and 648 timing-output scientific comparisons across the
two retained timing batches. The evidence is in `P36B-evidence/`, never target.
`complete-commands.json` records real exit codes; `provenance.json` records
binary, source and UPF hashes; `performance-spreads.json` carries all min/max
and median phase/wall numbers. README regeneration and freshness tests pass
without changes to existing scientific reference material.

Temporary executable copies `P36B-evidence/dft-cli-*` are removed after all
comparisons, as requested for disk hygiene. `binary-cleanup.json` retains their
hashes/sizes and names. Raw inputs/results/densities/logs remain; the final
executable remains in the single existing `target/release/` build directory.
No extra target directory was created, and no unrelated worktree files were
deleted. Reproduction after cleanup requires rebuilding the baseline and
candidate sequentially in that same target directory and supplying those
executables to the retained drivers. For the original baseline, restore the
eight changed production files from the recorded HEAD in a temporary local
source transaction, build, then restore the candidate in a finally block;
never use an extra target directory. For the budget-only arm, apply only the
recorded full-pool inner-budget change to that baseline. The engine-only
fixture baseline driver already restores the three affected engine files in
this fashion. The saved source hashes identify the completed candidate.
