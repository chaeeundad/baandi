# P47A switchable Davidson preconditioners

## Frozen rule R (before measurements, 2026-11-08)

Default remains unchanged. Systems: si8, cubic si16, skewed Si16,
si32 Gamma (real, standard and strict), k222, k444s, PseudoDojo MgO and Si,
one low-symmetry/elongated cell, and a synthetic multi-projector fixture.
Each system: five standard starts with relative density noise 1e-10,
seeds 401–405, plus one strict certified run; pair each variant with default.
Single-thread alternating order. Record Davidson steps, H vectors, multiply-adds,
SCF iterations, wall seconds, energies in Ha, certificates, input/build hashes.

Recommend only if median H vectors rise by no more than 3% on every system
and fall by at least 8% in the median over systems; SCF iterations never rise;
no default-converging start fails; strict certifies; delivered energy error is
at most 5e-9 Ha and at most 1.5 times default error (strict also within 1.5 times);
real and complex agree on direction. Stop a variant at the first failed
condition and report the system. No post-measurement tuning of constants or R.
An independent tighter reference is required to establish delivered errors;
variant-minus-default differences alone cannot establish this gate.

## Architect amendment R1 (2026-11-08, before resumed evaluation)

The architect explicitly corrected a specification error in R after seeing
only the Si8 complex strict selection result. Strict accuracy now requires
|E_variant - E_reference| <= max(1.5 |E_default - E_reference|, 1e-12 Ha).
The 1e-12 Ha floor is the documented P17A strict plateau. The standard-profile
accuracy condition and every other condition remain unchanged; no variant
constant is changed. This is an openly declared architect amendment, not an
implementer adjustment to rescue a measured result.

At amendment time, performance selection had evaluated only Si8 complex strict;
no five-start standard ensemble or other system had been evaluated for selection.
The separately labelled identity/determinism validation runs had also occurred
and their raw standard-profile outputs exist; they were not used for selection.
The old stop below is preserved as historical evidence and superseded by R1.
Both Si8 strict accuracy checks pass R1. Evaluation resumes with seeds 401–405,
then proceeds through the required systems only while a variant survives.

Provenance clarification: the stored tight reference was already a separate
run, with 20 SCF iterations versus default strict's 14; its computed energy
happened to equal the default energy exactly. It was not the default run reused
as a reference. Resumed screening reruns that independent tighter calculation.

## Derivation and fixed constants

For normalized plane-wave coefficients c, Parseval gives
<T> = sum_G |c_G|^2 T_G / sum_G |c_G|^2, T_G = |G+k|^2/2 Ha.
In orthonormal real coordinates the paired real/imaginary coefficients have
squared norms summing to the same quantity; each pair has the same T_G.
The per-band scale is s = max(<T>, 0.1 Ha). The floor handles a pure G=0 band
without division by zero and reuses the existing energy scale, not a fitted
constant. Sums follow coefficient order independently of thread budget.
TPA correction is -K(T_G/s) r_G. K is the requested rational polynomial:
(27+18x+12x^2+8x^3)/(27+18x+12x^2+8x^3+16x^4).
For x>1 evaluate after dividing both polynomials by x^4 to avoid overflow.
K(0)=1, K(1)=65/81, K(x)~1/(2x); it is positive for finite nonnegative x.

The diagonal approximation h_diag remains T_G (no potential terms added).
Ritz-shift correction is -r_G / max(h_diag_G-theta, 0.1 Ha).
The denominator is >=0.1 Ha even when h_diag-theta is negative or zero;
therefore it cannot flip the correction's sign relative to -r.
Reject nonfinite differences. The existing default already subtracts theta,
but preserves the sign of the floored denominator; claiming a newly introduced
Ritz shift would be incorrect. The new mechanism is positive clipping.

A14 added mean-local and nonlocal diagonal terms with smoothing and lost;
this experiment changes damping or sign clipping while keeping the kinetic
diagonal. A13 capacity, A19 band counts, A15 stopping and A21 tolerances are
unchanged. Context: diagonalization-literature.md section 3/rank 3 and
qe-study/qe-solvers-and-subspace.md (prose only; no QE source read or copied).
This is not a reproduction of QE's preconditioner.

## Compatibility and evidence

Use `scf.preconditioner = "default" | "tpa" | "ritz_shift"`.
Legacy solver entry points keep default behavior. Nondefault options are
restricted to unpolarized, fixed-occupation Davidson. Checkpoint omission
means default; exact restart compares the selection. Results record it.
Raw evidence lives in `P47A-evidence/`, outside target. Mac timing is observation;
dev08 acceptance must go through the orchestrator. No recommendation yet.

## Historical R screening result: stop at Si8 strict (superseded by R1)

First evaluate an independently rerun tighter SCF (density 1e-10, energy 1e-12,
eigen 1e-11 in the engine's respective units), then default strict and each
variant strict, before spending the five standard starts. No rule or constant
was changed after this measurement. This reference uses the same Hamiltonian,
not an external code, and does not certify its own absolute total-energy error.

| Si8 complex strict | SCF | Davidson steps incl. certificate | H vectors incl. certificate | Logical GEMM multiply-adds | Energy (Ha) |
|---|---:|---:|---:|---:|---:|
| Tight reference | 20 | 659 | 2848 | 8546640594 | -31.3218036010317320 |
| Default | 14 | 318 | 2022 | 4560829350 | -31.3218036010317320 |
| TPA | 14 | 268 | 1944 | 3877466496 | -31.3218036010317817 |
| ritz_shift | 14 | 266 | 1843 | 4095435631 | -31.3218036010316858 |

All four certify. H-vector changes are -3.8576% and -8.8526% on this one
strict trajectory. These are not ensemble medians. The computed default error
is zero; variant errors are 4.973799150320701e-14 and 4.618527782440651e-14 Ha.
Both therefore fail the literal 1.5-times-default comparison. Finite precision
prevents interpreting this as established deterioration in delivered accuracy.
No post-hoc denominator floor was added to rule R. **Screening stopped; neither
variant is recommended.** Standard five-start ensembles, the remaining systems,
and performance-direction agreement between real and complex are not established.
Subsequent tests are implementation/identity validation only, not a resumed
performance study.

The instrumented SCF-loop wall observations (seconds) are reference 3.88394,
default 2.31613, TPA 2.10866, ritz_shift 2.10906. They exclude setup and use
GEMM accounting with clocks; they must not be treated as acceptance timings.
Raw records: `P47A-evidence/si8-screen.log`.
Reproduce the stopped screen (external UPF, no UPF committed):

```sh
P14A_WORK_ACCOUNTING=1 P14B_GEMM_ACCOUNTING=1 \
P47A_SI_UPF=/absolute/Si.pz-vbc.UPF P47A_CASES=si8-gamma \
cargo test --release -p dft-engine --test preconditioner_study sweep \
  -- --ignored --nocapture --test-threads=1
```

No dev08 request to establish speed is justified after this stop. The
orchestrator should retain the negative/inconclusive result and unchanged
default; a future study first needs a resolved independent error oracle.

## Implementation contracts and validation

The complex `DiagonalPreconditioner::with_variant` builder preserves existing
construction call sites. Davidson passes current Ritz columns plus selected-root
indices through an additive trait method; other preconditioners inherit the
old residual-only call. TPA does not approximate the band scale from residuals.
The real `solve_reusing_with_preconditioner` method is additive: existing
`solve_reusing` and all older entry points explicitly choose default. No new
per-step allocation or cross-worker reduction is introduced by either variant.
The real workspace caches storage, not the selection; each call reads the option.
Runtime requests to switch a nondefault preconditioner to the hybrid solver
are rejected before policy mutation. Spin and FD default paths are untouched.

`cargo test --release -p dft-solvers --test preconditioner_variants` passes
2 tests: selected-band scaling, the pure-G=0 floor, stable large-x evaluation,
positive Ritz clipping, dimension rejection, an analytic tridiagonal spectrum,
and complete eigenvalue/vector-bit identity at T=1,2,4,8,12 in both representations.
The analytic eigenvalue tolerance is 1e-10 Ha; this is not a material benchmark.
Checkpoint/config tests cover all three wire values, default omission, invalid
names, and exact-restart rejection when the selection changes.

`scripts/check.sh` exits 0 (build/test/fmt/clippy all PASS).
`python3 -I scripts/tests/test_readme_figures.py` exits 0 (9 tests).
No README result is promoted from this stopped study; README/figures and
artifacts/* remain unchanged, respecting this ticket's ownership.
Initial compile errors (a missing re-export and its malformed insertion) were
fixed before evaluation. Their logs and actual exit codes are retained alongside
successful final checks in `P47A-evidence/build-test-exits.json`.

The identity driver passes absolute binary/UPF paths to `bench_scf.py`.
It pairs the clean-HEAD baseline binary with the final binary on gamma, si8,
k222s and k444s, strict/standard, T=1/8, real off/on. For the k meshes, forced
real `on` is unsupported by design, so the real-enabled configuration is `auto`
and verifies the complex fallback; it must not be described as a real k-mesh run.
The separate determinism driver compares each variant against its own T=1
result on gamma/si8, strict/standard, real off/on, at T=2/4/8/12.
These are software validation runs after the performance stop, not further
selection data. Commands, statuses, summaries and hashes remain under
`P47A-evidence/`; no additional target directory was created.

Final validation outcome: all 32 clean-HEAD/default comparisons have exactly
identical energy bits, density SHA-256, Davidson steps, H counts and SCF counts.
All 40 variant benchmark commands exit 0 (80 system/run conditions), with
64 comparisons against T=1 exactly identical in those same quantities.
See `identity-summary.json`, `determinism-summary.json` and the complete command
manifests under `P47A-evidence/`. This confirms determinism on the stated scope;
it does not remove the rule-R stop or establish broader material performance.

## R1 resumed outcome: SCF-count rejection on Si8 standard

The separately executed tight reference again certified after 20 SCF iterations,
E = -31.3218036010317320 Ha. Default strict again certified after 14 iterations.
Both variants pass the amended strict accuracy bound of 1e-12 Ha; their strict
energies and counts reproduce the table above exactly. The earlier R stop is
superseded and is not used as evidence against either variant.

The standard profile retained its original conditions. Runs used T=1 and
charge-neutral relative density perturbations of magnitude at most 1e-10,
seeds 401 onward. Each pair alternates execution order: default then variant
at seed 401, variant then default at seed 402. Each variant stops independently
at its first failure; no parameters, tolerances or selection constants changed.

| Pair | Default/variant SCF | Default/variant Davidson steps | Default/variant H vectors | Default/variant logical GEMM MACs | Default/variant error (Ha) | Decision |
|---|---|---|---|---|---|---|
| TPA, seed 401 | 28 / 28 | 100 / 89 | 1253 / 1241 | 832166165 / 764930191 | 8.097416070995678e-10 / 8.478551194457395e-10 | Per-start conditions pass |
| TPA, seed 402 | 21 / 28 | 95 / 87 | 1011 / 1276 | 818686905 / 793251482 | 1.454594666938647e-9 / 3.293294525974488e-10 | STOP: SCF count increases by 7 |
| ritz_shift, seed 401 | 28 / 29 | 100 / 90 | 1253 / 1283 | 832166165 / 778379761 | 8.097416070995678e-10 / 1.508126956650813e-10 | STOP: SCF count increases by 1 |

Every measured start converged. Both failed pairs pass the unchanged standard
absolute and relative accuracy conditions: the failure is exclusively increased
SCF iterations. Fewer Davidson steps did not translate into fewer total H
vectors in those pairs. Do not apply the median-H condition to these individual
starts: five-start medians were not completed.

**Neither variant is recommended under R1.** Per the explicit early-stop rule,
TPA seeds 403–405, ritz_shift seeds 402–405, and the remaining systems were not
run. Thus Si16, skewed Si16, Si32 real, k222/k444s, PseudoDojo MgO/Si,
low-symmetry/elongated, synthetic multi-projector selection, and real/complex
performance-direction agreement are not established. These are intentionally
unreached gates after a decisive earlier failure, not claimed passes.

Instrumented SCF-loop times, default/variant seconds: TPA seed 401
1.120539708 / 1.102552500; TPA seed 402 0.917766250 / 1.121339584;
ritz_shift seed 401 1.120815500 / 1.135971458. As before, these exclude setup
and include accounting overhead; no speed claim or dev08 acceptance follows.
Raw evidence: `P47A-evidence/si8-amended-screen.log` and
`P47A-evidence/amended-screen-summary.json`. The release sweep exits 0 because
it successfully records the scientific rejection; exit 0 does not mean R1 passed.
