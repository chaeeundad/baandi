# P15A: frozen-Hamiltonian stopping replay

2026-11-05. Measurement only, derived from `qe-solvers-and-subspace.md`
hypothesis 1 and `qe-scf-cycle.md` E0–E2. Only BAANDI code, existing prose notes
and stored QE outputs were read. No QE source, new default, weakened certificate,
or accepted speedup is involved.

## Answer and its boundary

**The stopping definition can explain the magnitude of the step-count gap on
these fixed BAANDI trajectories; the data do not establish that BAANDI converges
more slowly per expansion at equal delivered accuracy.** Changing only the
stopping predicate, keeping every numerical threshold and strict/buffer prefix
unchanged, changes si8's predicted ordinary stopping work from **283 rows /
2063 H vectors to 39 rows / 749 H vectors**. That is 86.2% fewer recorded
Rayleigh–Ritz rows and 63.7% fewer H vectors. These are prefix counters, **not an
SCF optimization**: the discarded eigenvectors never feed the next iteration.
A residual norm and a Ritz change expressed in Ha are different guarantees.

The threshold schedule does not simply provide another saving. With QE-like
empty-root thresholds, the same numeric native tolerances give 38 rows / 729
vectors; replacing that sequence by the QE ratchet evaluated on the stored
BAANDI Hartree-metric errors gives **54 rows / 1000 vectors**. Thus in this
controlled chain the ratchet adds 16 rows, rather than removing them. Changing
the empty-root policy alone accounts for 39→38. These are specific to the
recorded 15-iteration si8 history, not a universal decomposition.

Across the first **14** native SCF iterations, baseline work is 257 rows.
Applying the **relayed QE threshold sequence** on those same frozen BAANDI
operators gives 89 rows / 1309 vectors, versus the relayed QE 83 expansions.
BAANDI counts the initial projection as row 1; removing one such row per solve
gives **243 native versus 75 replay expansions**. The 75 versus 83 comparison
is still not equal-accuracy or common-Hamiltonian evidence. It shows why the
historical 290–340 versus 83 totals cannot be assigned to a slower convergence
rate by subtraction. The current baseline is 283, not either historical total;
its fifteenth iteration and final certificate must not be silently folded into
QE's 14-iteration number. QE's 98 H calls are blocks, not vectors; their widths
are unavailable. No multiplication by 20 is used to invent QE vector counts.

## Independent accuracy check uncovered a second issue

**Small residuals do not prove that the lowest roots have all been found.** An
independent full-space dense oracle, with true residuals checked at 1e-10 Ha,
was added after the first replay exposed root insertions. On the initial si8
Hamiltonian even the all-root 1e-10 Ha residual replay misses one member of a
low degenerate manifold: maximum occupied sorted-eigenvalue discrepancy is
**0.110723486150125 Ha**. Native first-iteration occupied discrepancy is
**0.133012267588733 Ha**. This is not roundoff or evidence that the looser
criterion is safe. The tight replay is therefore not used as a ground-truth
spectrum at that snapshot.

The phenomenon is not silicon-only. Initial tight replay versus dense maximum
all-root discrepancies are 2.958860388 Ha (20-atom skew), 0.380087693 Ha
(4-atom elongated, here empty roots are affected), 8.225605511 Ha (16-atom
Gamma, also present at its mesh Gamma point), and 0.523220374 Ha (slab).
The converged native skewed run still has an occupied-root discrepancy of
1.514951753 Ha against the dense spectrum of its **own** final ordinary frozen
Hamiltonian. This is an existing spectral-completeness finding, not introduced
by the harness: native energy, density, eigenvalues, history and certificates
match the pre-change executable exactly. No repair or numerical-policy change
is attempted in this measurement ticket.

Dense agreement of the tight replay at si8 iterations 8 and 15 is respectively
9.88e-15 and 2.10e-14 Ha over all returned roots. All five selected k444s dense
checks agree within 7.28e-15 Ha. Dense checks cover every small synthetic solve
and selected large solves (IDs recorded below). All accuracy tables prefer the
dense oracle when available. Elsewhere a tight replay is explicitly only a
comparison spectrum; spectral completeness is not independently established.

At si8 iteration 8 the scheduled Ritz predicate stops at row 3 with occupied
residual as large as 8.00736e-4 Ha and occupied eigenvalue error 1.55524e-6 Ha;
the native terminal occupied error is about 4.04e-11 Ha. At iteration 15 the
scheduled stop has occupied residual 5.97767e-8 Ha but eigenvalue error about
1.37e-14 Ha. The relationship between Ritz change, residual and eigenvalue
error is therefore neither a constant conversion nor a fixed number of steps.

The equal-delivered-accuracy scan illustrates the distinction. At si8 iteration
8, for occupied error <=1e-10 Ha against dense, the best scanned native-trajectory
residual prefix uses 13 rows / 140 vectors (one root allowed), and the Ritz
prefix 10 / 128 (five allowed). These are hindsight selections, not validated
universal tolerances. At iteration 15 the target is already met at row 1 /
20 vectors; a consecutive-Ritz-change rule necessarily waits until row 2 /
37 vectors. At iteration 1 no scanned prefix qualifies. Full thresholds,
allowances, censored entries, and targets 1e-4 through 1e-12 Ha are retained in
the machine-readable table, including the independent tolerance replays.

## Generality, schedules and rate curves

[Generated tables](../../../tests/reference/perf/p15a/tables.md) and
[reproduction contract](../../../tests/reference/perf/p15a/README.md) contain all
seven systems, every ordinary iteration and k point, and snapshots at iteration
1, midpoint and last recorded iteration. si8/k444s are real NC-Si inputs;
20-atom skew, 4-atom elongated/slab, and 16-atom Gamma/2x1x1 are synthetic
soft multi-projector fixtures, not validated material predictions.

- k444s: 1686 native rows / 10290 vectors → 740 / 5631 at the ratcheted Ritz
  predicate, with no censored ordinary solve. This supports a substantial
  stopping-definition effect outside Gamma, not a general adopted schedule.
- The other five systems include censored native-prefix stops. Their partial
  sums are **not** complete predicted costs. Strict-replay schedule hits are
  also retained to show work on a longer, separately specified trajectory.
- Slab and 16-atom mesh reach the unchanged 100-iteration limit in both builds.
  Their late snapshots mean exhausted trajectories, not convergence. Gamma
  16-atom, skew and elongated run for 10, 6 and 9 iterations respectively.
- E_res first falls below 5e-9 Ha at si8 iteration 11 (native final-energy
  difference 9.01e-10 Ha), k444s 7 (8.05e-11 Ha), 16-atom Gamma 7
  (1.25e-10 Ha), skew 4 (1.42e-12 Ha), elongated 5 (5.89e-10 Ha).
  Neither exhausted trajectory reaches that threshold. These differences are
  against each native final state, which does not resolve spectral completeness.
- The native-history first-iteration redo condition fires on skew and the
  16-atom mesh. A redo would change vectors and dr2: only the trigger and
  implied threshold are computed, not an invented new SCF history or redo cost.

Per-root curves record both first crossing and the first step after the last
above-threshold value for residuals 1e-2, 1e-3, 1e-4, 1e-5 and 1e-6 Ha.
At si8 iteration 8 the all-root tight trajectory first reaches 1e-6 at rows
6–36 depending on root. At iteration 15 all roots initially satisfy 1e-6,
but some subsequently exceed it and only remain below after row 42. These
are sorted-root diagnostics, not tracked eigenvector identities. They show
that first crossing is insufficient justification for hard locking and that
empty-root recovery can dominate a strict all-root continuation. They do not
prove a locking algorithm would be unsafe with a margin, cluster rule and
unlocking certificate; none was tested.

The QE relay has average expansion counts but no matched per-root residual or
Ritz-error curves and no transferable frozen operators/starts. Consequently
**QE versus BAANDI convergence rate at equal accuracy remains unidentified**.
Calling the unexplained remainder a convergence-rate penalty would fabricate
evidence. This limitation applies even though the native prefix stopping
experiment itself is complete and deterministic.

## Hypotheses and next experiments

| Hypothesis | Evidence and disposition |
| --- | --- |
| H1: stopping strictness explains many extra steps | Supported for these si8 and k444s trajectories; 283→39 and 1686→740 under the stated different rules. Safe replacement is not established; early root completeness and delivered errors prevent that conclusion. A15's rejected strict-refinement implementation is not retried. |
| H6: a QE-like schedule improves SCF | Not established. Ratcheted native-history scheduling adds work relative to the same-numeric Ritz experiment on si8. E0 exposes different error metrics; A21's fixed-factor rejection stands. |
| Locking removes tail work | Re-crossings support measuring stability first. No hard locking implemented; neither its saving nor its rejection is inferred. A56's zero guaranteed-removable bound for the current unlocked rule stands. |
| Nonorthonormal basis accelerates convergence | Untested. These trajectories retain the existing orthonormal basis; no rate or speedup conclusion follows. |

Ranked next experiments:

1. **Spectral completeness and start subspace**: independently reproduce the
   dense mismatches, log degenerate-cluster multiplicities and initial overlap,
   and test a separate completeness/start safeguard. Keep the old path as oracle
   and preserve all existing certificates; this is a new algorithm ticket.
2. **Matched QE frozen-H accuracy traces**, through the orchestrator on dev08:
   per-root residuals, Ritz values, H block widths, initial-projection accounting,
   and exact input/UPF/output hashes. Use si8 early/mid/late and non-Si fixtures.
   Until then report the 83-expansion relay as context, not a rate oracle.
3. **Controlled SCF coupling after completeness is resolved**: choose tolerances
   by delivered occupied-spectrum and density accuracy, then evaluate the
   Hartree-metric ratchet and a real first-step redo. Keep the strict default;
   this replay does not promise unchanged SCF iteration counts.
4. **Locking stability before structural algebra changes**: measure cluster
   re-crossings, margins and unlocking, then test hard locking alone. Only after
   that, isolate nonorthonormal-basis conditioning/rank diagnostics on the skew
   fixture. Do not combine the mechanisms or repeat rejected capacity/QR trials
   without their required new premise.

## Implementation, provenance and verification

Baseline/current-main source is `ce64b57437e2cb7b14804a4bdd15fb252dab6286`;
P15A changes are uncommitted on `task/P15A-REPLAY-STOPPING`. There is no older
performance-baseline speedup claim. Executable and source hashes, exact commands
and gate exits are in `tests/reference/perf/p15a/verification.json`.

`P15A_REPLAY=1` at build time enables capture and synchronous frozen replays;
`P15A_DENSE_IDS=0,7,14,128,288` selects additional full-space oracles. These IDs
cover si8 1/8/15 and k444s k=0 at 1/5/10, plus two early k444s points. Dimensions
<=128 use dense automatically. The diagnostic dense path is bounded at 4096;
the ordinary DenseEigensolver 1024 guard is unchanged. No new dependency or
manifest edit is required. Normal builds compile out tracing and extra solves.

The real-run harness validates exact energy, density bytes, eigenvalues, entire
scientific convergence history and certificates against the saved pre-change
CLI. Repeated measurement runs reproduce every trace and dr2 exactly. The
unchanged bench_scf identity gate covers gamma, k222s, k222, k444s, si8 at
T=1,2,4,8,12; all energies, density hashes and work/certificate counts also
agree across these budgets. The analytic signed-mode test checks FFT scaling,
Ry/Ha conversion, sphere filtering and omission of G=0. An independent operator
wrapper checks trace H counters; a dense spectrum checks the unit fixture.

Full build/test/fmt/clippy and README figure tests pass. The README evidence
block has no new accepted performance claim; no generated figure/result was
changed. Initial tooling failures (relative UPF path, then a misplaced tracing
insertion caught by the compiler) were corrected and are recorded in the gate
manifest. Mac wall times overlap other work and support no timing claim. MKL
and dev08 were not executed. Dashboard/session uploads remain outside this
ticket's owned files; the canonical delegate registry tracks this run.
