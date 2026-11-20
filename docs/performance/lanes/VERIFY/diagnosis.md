> Reopened by the orchestrator: attempt timing is now allowed. The closure
> below is historical and superseded for this scope by [first-schedule.md](first-schedule.md).
> Judgment and certificate thresholds remain unchanged.

# LANE-VERIFY Phase 2 diagnosis — no qualifying implementation

Decision: **close this Phase 2 attempt without a new optimization switch**.
The requested expected saving of at least 0.2 s on Si64 MPI8 is not supported
by the relayed profile for any new lever inside this lane's current ownership.
This is a decision about the available evidence and scope, not a proof that
verification cannot be accelerated. The acceptance rule, attempt timing,
tolerance, schedule, continuation and defaults remain unchanged.

## Evidence and provenance

The user relayed six instrumented dev08 Si64 MPI8 runs at
`bad02ee79a02322a5be0969fda200617ac8f53f9` on 2026-11-20. Numbers below are
**relayed**, not reproduced on the Mac. Full inputs, binary hashes and raw
per-rank JSON were not supplied with this relay. The one-core profile is still
queued; no one-core saving is inferred from MPI measurements or the older
4.84 s one-core budget context. No QE comparison is made.

Local source HEAD is `17f93a38221e0db16c0fa4f606ecf7ef500a375b`. A git comparison
against `bad02ee` confirms `session.rs` and `gamma_real.rs` are identical.
The relevant distributed solver differences add lane profiling spans; the
kept-start checks, basis initialization, Ritz-change test and return product
semantics examined here are unchanged. Both revisions default VERIFY_REUSE to
`half`, VERIFY_SCHEDULE to `on`, M2_REUSE to `on`, and REUSE_PSI to `on`.
Effective switches and cache-taken counts in the actual remote run must still
be checked in its raw output; defaults are not proof that every cache hit.

External wall: 8.4884–8.6852 s, median 8.5466 s. Verification pass+fail median:
1.0475 s, or 12.26% of that median wall. All six runs have 11 ordinary fine
maps and two verification attempts:

| Attempt | SCF iteration | Verdict | Coulomb error (Ha) | Density residual (Bohr^-3/2) | Davidson iterations | H columns |
|---|---:|---|---:|---:|---:|---:|
| 0 | 10 | fail | 1.009e-8 | 7.053e-8 | 6 | 380 |
| 1 | 11 | pass | 2.822e-9 | 3.429e-8 | 5 | 388 |

The relay separately quotes fail 0.5439 s, pass 0.4959 s and mean per-attempt
0.5199 s. These do **not** sum to the quoted 1.0475 s median: preserve their
separate aggregation provenance. Schedule: one evaluation, predicted extra
maps zero. Ordinary fine-map reference cost is about 0.52 s (5.77 s / 11).
Delaying an attempt exchanges work of similar cost and also removes a useful
failed-map continuation product. It is not a proposed optimization.

Reported verification map-phase rank maxima:

| Phase | Seconds | Reported asymmetry |
|---|---:|---:|
| Solve | 0.8995 | 0.4% |
| Density | 0.1026 | 21.7% |
| Energy | 0.0499 | 46.1% |
| Potential | 0.0151 | not supplied |

These add to 1.0671 s, **not** 1.0475 s. The reported 84.3% solve share uses
that phase sum. Different maxima/summary statistics cannot be added into a
coherent critical-rank wall or subtracted from the external median. Bucket
asymmetry was reported as zero; that does not make every nested phase balanced.
Use the Phase 0 parser's same-rank closures when the raw JSON is available.

Solve attribution supplied: solver_total 0.8217 s, apply 0.3853 s,
reduced_gram 0.1548 s, rotate 0.1317 s, reduced_eigen 0.0984 s and kept_update
0.0717 s. solver_total is a parent interval, not another additive stage.

## What an attempt computes

ADR 0043 V1–V4 and `session::verify_fixed_point{,_inner}` specify:

1. Install the delivered output density (the owned representation on MPI,
   restoring the ordinary loop input on return), and prepare its warm orbitals.
2. Run the ordinary Kohn–Sham map at that density. Build its Hartree/XC/local
   potential; warm-start Davidson with the **preceding SCF map's output**;
   use `StandardProfile::ratchet(tau, standard_error, electrons)` and the same
   occupied/empty-band tolerance policy. This is a Ritz-change solve, with
   returned eigen residuals for the guard and informational certificate.
3. Build the density of the returned orbitals, respecting occupations and
   symmetry, and perform the map's energy audit. The map already has shared
   and cached pieces; it is not a fresh problem/FFT-plan/Ewald construction.
4. Evaluate the Coulomb metric between this output and the delivered input.
   In the owned route, `stop_error_owned` uses kept Fourier coefficients when
   their **input/output epochs match**, otherwise it materializes/transforms.
   The ordinary stop estimate for a different density pair is not this metric.
5. Evaluate the occupied weighted eigen-residual estimate, density distance,
   occupied eigenvalue change, sampled-gap and applicable clipped-charge checks.
   V2 is the unchanged conjunction of the Coulomb test and enabled residual
   guard, plus the safety checks. Strict residual/eigenvalue criteria are
   computed and reported as V3 information; they are not substituted for V2.
6. Record work and certificate quantities. A pass delivers the preceding SCF
   state, not the verification state. A failure continues via the existing
   HALF reuse path. No result of this analysis changes the judgment or timing.

The residual guard is a separate empirical estimator, not an interchangeable
additive error allowance for the Coulomb test. The small residual term cited
in historical F-VERIFY-COST does not prove that a looser inner solve leaves
F(rho), its Coulomb metric, or the later HALF continuation unchanged. That
historical ticket's budget-tolerance proposal is not authority to alter V1
under the current user instruction. Its ledger closes the cost A/B as
superseded by VERIFY-SCHEDULE; its proposed design file is absent here.

## Arithmetic and bandwidth floors

Let P be the local coefficient-row count, B the returned band count, M a
Davidson search width and N the local real-space grid size. For real f64 rows,
one retained X/H X pair requires at least 16 P B bytes; retaining V/H V at
width M requires 16 P M bytes, before planes, scratch and reduced matrices.
A dense real overlap of two P-by-M blocks costs approximately 2 P M^2 flops;
a rotation from M to B columns costs approximately 2 P M B flops per block.
Distributed overlaps also need their reductions. The changed local-potential
product requires touching the retained planes and new/old potential, plus
transforms where the representation requires them. Keeping H X avoids the
full operator but cannot make this update free.

A fresh density must accumulate the occupied orbital contributions over the
grid (work scaling with N B) and communicate/distribute the resulting density
as required. The Coulomb metric needs the **current pair's** Fourier difference
and weighted norm; retained coefficients avoid transforms only when their
epochs match. These are dependencies, not a claim that current kernels attain
machine bandwidth or arithmetic limits.

The relay supplies total H columns (380 and 388), not per-call widths, local
row/grid sizes, retained-space sizes, measured bandwidth or sustained flop
rates. Consequently it cannot support a numerical hardware floor in seconds,
or a claim that all 768 columns are redundant. Phase 0 size/count records and
raw per-rank profiles are needed for that calculation. The phase costs below
are budget screens, not roofline estimates.

## Product lifetime: what is already reused

Let the ordinary iteration-10 map be `rho10 -> s10 = F(rho10)` with orbitals
`X10`. Attempt 0 computes `v10 = F(s10)` and orbitals `Y10` at `H[s10]`.

| Product | Existing consumer / validity |
|---|---|
| Failed attempt orbitals Y10, including lazy local rows | Both owned and full-density continuation branches select `split_points(verification_orbitals)` and set `previous`/`previous_lazy`. Iteration 11 starts from Y10, not an atomic/free-electron restart. |
| Pair (s10, v10) | HALF adds this pair after the ordinary pair to Pulay history, then uses the midpoint of the two proposals, with the unchanged feasibility safeguard. This already uses the failed map as numerical information. |
| H[s10] Y10 | The distributed keeping return stores X and H X. `prepare_m2_kept_start` updates H X by `(V_new - V_old) X`, using kept real-space planes when available; the solver checks the start/shape/overlap before taking it. It skips the full initial H apply on a hit. |
| Plane-space orbitals | The density build can retain occupied-pair planes for the next local-potential update. `kept_update` is a real update cost, not evidence that reuse is absent. |
| Attempt density Fourier products | Used for the corresponding epoch-matched Coulomb check. The subsequent mixed density is generally a different pair, so blindly reusing the old metric would change the judgment. |
| Final residuals/eigenvalues | Valid for their own operator, eigenpairs and tolerance. They cannot certify a new Hamiltonian without recomputation/update. |
| Full Davidson correction/search history V, H V, reduced matrices | Workspace allocations persist, but the next solve's logical basis is rebuilt from the returned band block X. The entire search space and its Ritz convergence history do not continue across maps. |
| Energy audit of the verification output | Does not replace the delivered ordinary-map energy. STOP_PREDICTOR diagnostics can also read it. Omitting it needs separate consumer/error-path analysis, not substituting an old density's audit. |

After mixing, iteration 11 computes `rho11 -> s11` using Y10 as the warm start.
Its own returned `X11, H[rho11] X11` become the most recent kept block. Attempt 1
then starts from X11 at `H[s11]`; H X is updated for that changed potential.
It does not directly resume attempt 0's full Davidson iteration state, but it
**does inherit its useful orbitals through iteration 11**. Restoring Y10 again
instead would change the selected start and could defeat the most recent cache.

Sources: `session.rs` failed-attempt `verification_reuse`, owned/full `feed`
and `previous` branches; `gamma_real.rs::prepare_m2_kept_start` and its caller;
`distributed_davidson/workspace.rs::solve_core_ws` kept-start branch;
`workspace/lean.rs::run` and `keep_return`. The kept-start branch consumes the
armed product, checks agreement/overlap, copies X/H X into the initial basis,
and suppresses the initial `apply` only on a valid hit. Current raw cache-hit
counts are not in the supplied relay, so this is a verified source contract,
not a fabricated per-attempt cache-hit measurement.

## Why five iterations and 388 H columns remain

A small fixed-point residual is not bitwise equality of potentials. Between
attempts, the failed-map feed and the ordinary map change both the input density
and the orbital block. The new tolerance is freshly ratcheted; the relay lacks
its exact per-attempt values, so no tighter-by-a-specific-factor claim is made.
A corrected H X gives a valid **starting subspace**, not already-converged
Ritz eigenpairs of the new operator.

The distributed workspace starts with the band block X (on a valid kept hit)
and prepares `previous` Ritz values as NaN. Lean convergence checks
`iteration > 1` and successive Ritz changes against root-specific thresholds.
It must therefore establish current-operator changes and add correction
blocks where needed. At return it rotates vectors/H images, forms residual
norms and keeps the returned band block, not all accumulated corrections.
Carrying old eigenvalues as if they were a previous iteration of the new
operator would change the stopping semantics.

The **observed five is not a proven minimum**. Exactly which roots require
which corrections cannot be recovered from two aggregate iteration/H counts.
H columns sum variable-width operator applications across corrections and any
rebuilds; they are not iterations times a fixed band count, nor 380 unique
persistent useful directions. Fewer iterations with 388 rather than 380
columns is therefore consistent with wider corrections. The supplied relay
has no per-iteration selected widths, restarts or cache-decline sequence to
attribute those 388 columns more precisely. Local-potential-only kept updates
also have their own cost/counters and must not be counted as full H columns.

## Floors and candidate screen

The 0.2 s target is 19.09% of the reported verification median, or 2.34% of the
reported external median. It is an absolute **MPI8** request, not a license to
compare or extrapolate between parallel models. The following are deliberately
optimistic work-envelope calculations on the supplied aggregates, not measured
savings or additive critical-path bounds:

| Lever | Evidence / economic screen | Decision |
|---|---|---|
| Carry failed orbitals and H X into the next map | Already default HALF + M2_REUSE. A second implementation supplies no new work saving. Ledger warm/feed variants are not reopened by the same premise. | Do not duplicate. |
| Remove verification energy work | Even deleting the entire reported energy phase offers only 0.0499 s. Some audit products/consumers remain; actual removable cost is lower. | Below 0.2 s alone; no implementation. |
| Eliminate density/energy imbalance | At most illustrative max-to-min phase gaps: 0.1026*0.217 + 0.0499*0.461 = 0.0452681 s, assuming the supplied asymmetry definition. This is not proven removable work or a coherent wall saving. | Far below target; do not claim a latency fix from phase percentages. |
| Reuse all density/potential results | Entire reported density+energy+potential envelope is 0.1676 s; actual fresh verification density and potential remain necessary. No identity witness permits substituting the preceding map's products. | No supported 0.2 s gain. |
| Remove warm copies or all kept-update cost | Lazy matrices carry no full complex payload to copy; the complete reported kept-update component is only 0.0717 s and much of it is the required delta-potential application. | Not a qualifying new lever. |
| Cheaper early failing attempt, or wait one map | Same judgment/time restriction; a partial solve is not the specified F map and cannot supply the identical HALF continuation. A replacement ordinary map costs about 0.52 s. | Excluded; no change to tau, rule, trigger or schedule. |
| Persist the full search subspace / correction directions across changed Hamiltonians | Potentially attacks solve work, but new H V must be computed/updated for **all retained directions**, their validity/overlap checked, and the reduced problem rebuilt. Old return X/H X reuse already exists. Needs solver/bridge ownership, storage/work counters and an independent correctness study. No observed redundancy or contraction trace supports a >=0.2 s prediction. | Unqualified hypothesis, not an implementation or performance claim. |

Solve work is the only demonstrated budget large enough: a solve-only route
would need approximately **22.23%** of the reported 0.8995 s, or an apply-only
route **51.91%** of 0.3853 s. Gram+rotate alone would need about **69.81%** of
0.2865 s. These are required reductions, not expected gains. Linear conversion
of those percentages into H columns is unjustified without widths, operator
mix and per-call costs. General H/FFT/algebra changes belong to their owners;
no single-material constant or generic switch is invented here to meet a
number.

## Conclusion

No new scoped lever clears the requested expected-gain threshold with the
unchanged acceptance/timing contract. Choose the user's authorized closing-note
branch, rather than label a solver-recycling hypothesis a measured saving.
See [closing-note.md](closing-note.md). One-core interpretation remains pending
its own relay; do not scale these MPI fractions to that model. Phase 0
instrumentation and its correctness evidence remain available unchanged.
