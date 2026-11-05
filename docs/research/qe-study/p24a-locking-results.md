# P24A: hard locking and one-rectangle Davidson experiments

2026-11-06. Derived from `qe-solvers-and-subspace.md`, hypotheses 2 and 5;
BAANDI source and prose notes only. No QE source was read or copied, and no new
QE run was performed. Baseline/current HEAD is
`46dfeddd370f5aaf37c898aca41735269c36cc05`, branch `task/P24A-LOCKING`, uncommitted.

## Decision

**Keep `[scf] locking = "off"` as the default.** `"auto"` is an experimental
strict-Davidson path, combining prefix hard locking and one-rectangle reduced H.
It is not an accepted performance optimization. Neither component is enabled in
the off path. Unknown values are rejected; combining auto with the standard
Ritz-change profile or hybrid solver is rejected. Runtime hybrid transitions
are also rejected. Checkpoint v1 cannot retain this policy, so auto checkpoints
are rejected instead of silently resuming with locking off.

The original A56 zero-removable-work statement concerned *unlocked* roots.
This experiment changes that premise by explicitly freezing vectors and their
H products. It does not repeat the rejected capacity or tolerance experiments
A08/A13/A15/A21. Capacity remains `max(48, 6*bands)`, bounded by dimension.

The experiment fails multiple adoption gates:

- `syn-s2` changes converged energy by **0.004416513597413463 Ha**, despite
  agreeing with the dense spectrum of its own final H. Iterations change 16 to
  14. This is not an allowed toy iteration-only difference.
- Si8 ordinary final solve at SCF iteration **15**, sorted zero-based roots
  **17 and 19**, exceeds the frozen 1e-9 Ha all-root dense cap. Maximum buffer
  error is **3.7218033654973226e-9 Ha**. Occupied error is 1.40443e-14 Ha and the
  fresh certificate spectrum agrees to roundoff. Approximate unmatched buffer
  values are not proof of missing eigenstates; the all-root validation still
  fails. No cap is relaxed.
- P15C `spectral_guard="final"` detects this Si8 discrepancy, repairs once, and
  passes at SCF **16**, using **2265** additional guard H columns. Ordinary SCF
  counts are therefore not preserved with that guard enabled. The independent
  guard remains unchanged.
- Synthetic 20-atom skew diverges at SCF **8** (density-growth limit), compared
  with the known incomplete legacy six-iteration success. Dense checks on the
  candidate show missing roots **2,10,11,12** at iteration 1 and **12** at
  iterations 6, 7 and 8 (zero based). At iteration 8 the occupied discrepancy
  is **0.7084692247412576 Ha**, with fixed-H occupied band-sum excess
  **3.4444200673607286 Ha**. These are not corrected SCF total energies.
- Whole-run MAC reductions are only **6.00% Si8 / 12.51% k444s**, below the
  note's 25% target. H work rises **2265 to 2809 / 11181 to 13500** columns,
  exceeding the 2% work-increase gate.

## Implemented algorithm and allocation boundaries

Only a contiguous lowest-root prefix can lock, in whole clusters. Adjacent
Ritz values belong to one cluster when their gap is no greater than
`max(10 * strict_tolerance, 64 * epsilon * max(1, abs(next_value)))`, in the
operator's eigenvalue units (Ha for the engine). The entire reduced spectrum
is inspected so that a cluster crossing the requested-band boundary is never
partly locked. Every member must have residual no greater than its own strict
or buffer tolerance divided by `sqrt(number_of_requested_bands)`.

Locked vectors and H products are stored once. The active reduced problem
excludes them. Ritz rotations and residual formation operate on active roots
only. The combined projection basis retains a fixed locked prefix: subsequent
steps do not recopy that prefix. Corrections remain orthogonal to both locked
and active columns using the existing two-pass projection machinery.

A new lock contracts the active basis to its retained Ritz vectors; ordinary
expansions thereafter preserve the active basis until capacity requires a
restart. This conservative lock-transition contraction is a design limitation:
it discards accumulated active directions. It is distinct from a capacity
retune, and is one possible source of the increased step count.

At a candidate return, combine all requested vectors, construct an **independent
full reduced product without conjugate fill**, perform one all-root rotation,
and freshly apply H to every returned root. Every root must meet its original
threshold. On a failed audit, record `(iteration, sorted_root, residual)`,
unlock all locked roots, disable further locking for that solve, and continue
from the audited full block. No residual or SCF certificate is weakened.
A fault-injection test corrupts only an initial cached product and verifies
actual unlock, continuation, and final dense agreement.

Tall arrays, reduced matrices, EVD scratch, projection output, selection arrays,
history and bounded failure storage are reserved before the loop. The new path
uses reusable `HermitianWorkspace::solve_into`, `reduced_hamiltonian_into` and
`ProjectionWorkspace::orthonormalize_into`. Stalled random corrections reuse
reserved storage. Normal numerical iterations introduce no allocation, clock,
logging or I/O. Operator implementations retain their existing overwrite
contract; an external operator that itself replaces output storage is outside
this solver's allocation guarantee. Return-owned metadata and source-only
measurement collectors are outside the numerical work policy.

The existing allocating adapters remain available for the legacy path. A unit
test proves variable-dimension EVD output is bit-identical to its original
adapter while preserving output storage. Projection failure paths restore a
valid output matrix even when a backend check fails.

## One rectangle and independent validation

When the active basis grows from p to q columns at N plane waves, preserve the
old p by p block and compute new rows against all q columns. Conjugate-fill
old/new entries, copy the lower new/new triangle to the upper triangle, and
make the diagonal real. Incremental MAC work is `N*(q-p)*q`, compared with
`N*(q*q-p*p)` for the legacy two rectangles. Initial/restarted products still
cost `N*q*q`. Hermiticity of the filled matrix is not an independent check.

Tests compare every entry against an independently evaluated full product for
complex bases and a real diagonal Hermitian operator, dimensions 31/513/1025,
growing and unchanged prefixes, and T=1,2,4,8. The original non-Hermitian
incremental-product tests remain intact for the legacy function. The final
all-root raw product supplies a separate in-solve Hermiticity check.

## Work measurements, ablations and margin

[Generated tables](../../../tests/reference/perf/p24a/tables.md) and
[summary JSON](../../../tests/reference/perf/p24a/summary.json) contain site,
step, fixture, counter, energy and dense-oracle results. MACs mean **complex
multiply-adds**, one half of P14B's historical `2mnk` count; no wall speedup is
inferred from them. Source-only P24A rows add the solve and iteration index to
P14B call-site shapes. Step zero denotes setup. Tail means the last quarter of
numbered steps in each solve, including terminal audit/rotation.

| Variant | Si8 whole-run MACs | k444s whole-run MACs | Si8 SCF |
|---|---:|---:|---:|
| Off, two rectangles | 4,744,613,419 | 1,129,003,567 | 15 |
| Auto, one rectangle | 4,459,794,280 | 987,797,005 | 15 |
| Locking, two rectangles | 5,212,958,324 | 1,079,833,280 | 15 |
| Rectangle only, no locking | 4,855,582,763 | 1,029,000,526 | 16 |
| Auto, twice the margin divisor | 4,651,704,502 | 1,045,720,012 | 15 |

The rectangle-only variant changes Si8 SCF 15 to 16 and reaches the same
alternate syn-s2 energy as combined auto. Thus locking is not the sole source
of the syn-s2 trajectory difference. No cause beyond the measured arithmetic
sensitivity is established. Two-rectangle locking also fails that energy gate
(delta 0.0050426730310153545 Ha).

The stronger `1/(2*sqrt(m))` margin costs more MACs on both principal cases and
does not recover syn-s2's reference energy. The experimental option retains
`1/sqrt(m)`; neither variant justifies default adoption. The gap criterion is
an exploratory rule tested on exact/near-degenerate clusters, not a universally
calibrated physical gap. No threshold was relaxed to make a gate pass.

Prefix locking is rare on these trajectories: the mean tail Ritz width is
**19.717 versus 20** for Si8 and **7.729 versus 8** for k444s. Si8 tail mean MACs
actually rise **8,904,781 to 10,665,880**; k444s falls **573,523 to 455,064**.
The measured implementation does not realize the note's assumed large locked
fraction. All-root final audits and lock-transition restarts add work. Preserve
this negative result before considering a different deflation structure.

## Generality, determinism and limits

The five Si baseline cases retain SCF counts (10/10/10/10/15) and differ in
energy by at most **1.13687e-13 Ha** with auto and guard off. Off retains exact
baseline energies and density bytes at T=1,2,4,8. Auto converged trajectories
also agree bit for bit across these thread budgets within the pure backend.
Timers, input thread identity and diagnostic completion order are excluded.

The robustness harness completes with **16 converged / 3 existing exhausted**
cases for both policies. Syn-s2's energy fails; other converged energies meet
1e-12 Ha. Toy iteration changes are syn-gaussian-fcc 17 to 15, syn-odd-shift
9 to 10, syn-s1 14 to 13, syn-s2 16 to 14. The separate P15A elongated fixture
changes 9 to 8 iterations with roundoff energy difference; Gamma 16 remains 10.
The slab exhausts 100 iterations in both paths with different trajectories;
mesh 16 diverges at iteration 8 rather than exhausting 100. Failed trajectories
are not equal-state performance comparisons. Error logs and dense snapshots
are retained, including the known pre-existing skew completeness problem.

All six external PseudoDojo Si (20/30/40 Ha) and MgO (40/50/60 Ha) fixtures
pass the 1e-12 Ha energy gate; maximum difference is
**8.24230e-13 Ha** (MgO 60). Five retain SCF counts; MgO 60 changes
**19 to 18**, another failure of the realistic unchanged-count gate. Final ordinary and fresh-certificate dense samples
at k=0, floor(nk/2), nk-1 pass the frozen 1e-9 Ha all-root and occupied band-sum
caps. Maximum occupied eigenvalue error is **1.24345e-13 Ha**. Inputs, UPF hashes,
commands, binary hashes and comparisons are stored; no UPF contents are stored.
The large-cell oracle scope remains sampled k points, not a full-mesh dense
proof. Thread checks cover all six material cases as well.

The exact P13B/P15C collinear Gaussian fixture also retains 12 SCF iterations
and exactly the baseline energy. Its final ordinary/certificate spectra for
both spin channels pass the dense cap (maximum 1.49519e-10 Ha over all roots).
The final guard passes for both channels, and T=1,2,4,8 reports agree bit for bit.
See `perf/p24a/spin/results.json`. Together with the 84 unpolarized thread
combinations, this supplies 88 recorded combinations; failed skew/mesh runs
match in exit status and scientific error history, not in a fabricated final state.


Ordinary build/test/fmt/clippy, MKL all-target type checking and README freshness
are recorded under `perf/p24a/final-gates/`. MKL is not linked/executed on this
Mac; its numerical and thread gates require dev08. Early failures (import/wire
plumbing, reusable-output shape, Clippy and attempts to run one-thread replay
at T>1) are retained with their fixes. The replay thread assertion was preserved;
production binaries supply the actual multi-thread evidence.

No existing baseline or tolerance was rewritten. README regeneration produces
no new result block: no accepted speedup is claimed. The README generator,
canonical delegate/dashboard artifacts and upload workflow are outside this
ticket's owned files. No sub-agent, remote change, commit or external upload was
performed. A recurring orchestrator report timer was not installed.

## Reproduce and dev08 handoff

See `tests/p24a_validate.py`, `tests/p24a_suite.py`, `tests/p24a_summarize.py` and
stored command manifests. Build the diagnostic executable with
`P14A_WORK_ACCOUNTING=1 P14B_GEMM_ACCOUNTING=1 P15C_ORACLE=1 cargo build --release -p dft-cli --offline`.
Use T=1 for dense/replay captures; use an ordinary release executable for
T=1,2,4,8 state identity. Set `DFT_NC_DIR` to the external PseudoDojo directory.
Ablations use `P24A_RECTANGLE_ONLY=1`, `P24A_TWO_RECTANGLES=1`, or
`P24A_MARGIN_DIVISOR=2` at build time; they are measurement controls, not runtime
or default policies. `--locking off` selects the rectangle-only ablation.

Request dev08 through the orchestrator only after review: baseline HEAD versus
reviewed commit, ordinary and native builds, pure and MKL at T=1,4,8,16, pinned
to one socket, MKL_NUM_THREADS=1; record binaries, energy/density hashes, SCF,
steps, H columns and P24A MACs. Repeat known failing cases rather than hiding
them. A 30-sample caller-resolved profile should distinguish added all-root
audit/lock restart cost from reduced algebra. No default-on or QE-exceedance
claim is justified by this experiment, irrespective of future local timing.
