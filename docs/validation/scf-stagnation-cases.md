# Measured synthetic SCF stagnation cases (P11M)

Three distinct structures exhibit persistent P11D stagnation/oscillation:
`skewed`, `slab`, and `slab-offset`. These are actual nonlinear SCF trajectories
on **nonphysical synthetic potentials**, not experimental materials or QE
references. The baseline receives 200 ordinary iterations, not the one-step
budget used in the original P11E demo. All three selected baselines still fail at
1,000 iterations. The symmetric slab remains unrecovered; that negative result
is part of the evidence.

## Reproduce

From the repository root, with the pinned toolchain and existing offline vendor:

```sh
python3 -I scripts/find_stagnation.py --output /tmp/stagnation-all
python3 -I tests/reference/stagnation/validate.py /tmp/stagnation-all
python3 -I tests/reference/stagnation/export_csv.py /tmp/stagnation-all
```

The output directory must be new. No dependencies are installed. A standalone
adapter under `tests/reference/stagnation/harness/` links the existing engine,
P11D diagnostics and P11B runtime. It builds release with `--offline --locked`
and shares this worktree's `target/`; it is outside the root Cargo workspace,
so its separate fmt/clippy checks matter. No solver source/default is modified.
The Python driver imports `scripts/gen_synthetic_upf.py:generate` with the exact
keyword arguments in each model's `generator` object. UPF content lives in a
temporary directory only; the result retains its SHA-256. No UPF is checked in.

One command reproduces each selected case, including its fresh reference,
untreated baseline and diagnostic-triggered intervention:

```sh
python3 -I scripts/find_stagnation.py --case skewed --policy linear-1.0-davidson --output /tmp/stagnation-skewed
python3 -I scripts/find_stagnation.py --case slab --policy pulay-0.3-davidson --output /tmp/stagnation-slab
python3 -I scripts/find_stagnation.py --case slab-offset --policy linear-1.0-davidson --output /tmp/stagnation-offset
python3 -I scripts/find_stagnation.py --case slab-offset --policy linear-0.5-hybrid --output /tmp/stagnation-hybrid
```

Append `--iterations 1000` to increase baseline/intervention budgets without
changing tolerances. The reference always has 1,000 iterations available.
Each child has a 115-second wall limit; a timeout is recorded as `WallTimeout`
with missing history, never inferred stagnation. In the stored main sweep the
largest measured process elapsed time was 0.557 s (rounded upward); this is only
local resource accounting, not a performance comparison.

## Frozen experiment and coverage

[Input matrix](../../tests/reference/stagnation/inputs/sweep.json) freezes 11
models crossed with seven policy tuples (77 baseline runs). This is a sparse
policy sweep, not the full Cartesian product of beta/mixer/solver values.
There are 11 fresh reference attempts and 12 intervention attempts, totaling
100 runs. The full sweep was repeated; every recorded numerical scalar,
certificate field, P11D evidence window, control request and disposition agrees
exactly after excluding timers. Arrays were not compared. The main sweep has
53 converged, 12 exhausted and 12 failed **baselines**; across all arms there
are 66 converged, 18 exhausted and 16 failed runs.

| Axis | Measured coverage |
| --- | --- |
| Cell/structure | cubic; skewed; elongated; slab with vacuum; offset-ion variants |
| Atoms | 1, 2, 4; one synthetic species per model |
| Sampling | legacy Gamma; 2×1×1, 2×2×1, 2×2×2 unshifted Monkhorst-Pack meshes |
| Bands | 4, 6, 16; 1, 2, 4 occupied bands respectively |
| Potential | valence 2; Gaussian radii 1, 3, 5 Bohr; projector channels [0], [0,1], [0,0,1], [0,1,2] |
| Certified reference gaps | 0.004163613860241472–0.4963928002684176 Ha; other candidates fail the integer-occupation gap gate |
| Policy tuples | Pulay β=0.3/1 and Linear β=0.1/1 with Davidson; Kerker-Pulay β=0.3 with Davidson/Hybrid; Linear β=0.5 with Hybrid |

All runs use cutoff 2 Ha, one requested worker, seed 2026, LDA-PZ, uniform cold
starts, history 6 and Kerker q0=1 Bohr^-1 where applicable. The density tolerance
is 1e-7 Bohr^-3/2; energy and eigen tolerances are 1e-8 Ha. Other options are
`ScfOptions::default()` and their complete Debug expansion is stored in every
result. Thus “default” below means **default mixing and Davidson**, not default
cutoff or automatically chosen band count. These tiny low-cutoff fixtures do
not establish converged basis-set physics or robustness at production sizes.

## Before/after results

Residuals are Bohr^-3/2, energies and energy differences are Ha. An iteration
count includes work before the intervention. A failed baseline has no certified
final energy; its last iterate energy remains in its history only.
For a terminal error, the final residual is the last completed boundary record,
not a fabricated record for the failing solve.

| Structure / baseline policy | Baseline result (iterations; residual) | Intervention boundary | After result (iterations; residual) | Certified after energy | Absolute difference from fresh reference |
| --- | --- | ---: | --- | ---: | ---: |
| skewed / Linear 1, Davidson | Exhausted (200; 0.0862306839476192) | 13 | Converged (39; 6.667116061439711e-8) | -0.6441418148101419 | 2.4202861936828413e-14 |
| slab / default Pulay 0.3, Davidson | Exhausted (200; 0.0530841058872408) | 70 | Exhausted (200; 0.01596379424805937) | unavailable | unavailable: reference also exhausted |
| slab-offset / Linear 1, Davidson | Exhausted (200; 0.11227601418830659) | 6 | Converged (31; 2.1375552115903586e-8) | -0.7215584248160892 | 1.8207657603852567e-14 |
| slab-offset / Linear 0.5, Hybrid | Exhausted (200; 0.042297496391354564) | 7 | Converged (33; 3.1922254864187616e-8) | -0.7215584248160788 | 7.771561172376096e-15 |

The single intervention is a typed P11B transaction: switch to Pulay β=0.1,
history=6, and reset history at the first P11D stagnation/oscillation boundary.
It binds exact observed state version/revision and expires at that iteration;
the experiment host supplies approval. The runtime records `Accepted` and
`Applied`, resets history and advances policy revision. No tolerance, model,
budget or eigensolver changes during this transaction. Reference attempts use
fresh Pulay β=0.1/Davidson runs. Energy agreement is required within 1e-10 Ha,
in addition to engine convergence/certification. This is a measured composite
intervention, not an ablation establishing which component causes recovery.

Five interventions meet this energy comparison. Only three recover a baseline
that exhausted its budget: the skewed case and the two slab-offset cases above.
The other two are Kerker-Pulay slab-offset trajectories that converge without
intervention in 45 iterations and with intervention in 37. A diagnostic is
advisory: it can trigger on a transient in an eventually convergent trajectory.
No default-recovery recommendation or speedup claim follows.

## Diagnostic evidence and budget controls

P11D `diagnose` consumes authoritative `IterationRecord` samples, tagged with
policy revision, at each completed boundary. Only triggered findings are stored;
no triggered finding does not imply every early window had sufficient evidence.
Complete six-residual/seven-energy windows are recoverable from numeric history;
`evidence_debug` stores the module's own evidence, including ratios and signs.

- **skewed:** first stagnation covers iterations 8–13, residual range and endpoint
  ratio both 1.1473666682516952. The 200-step baseline triggers stagnation at
  181 boundaries and oscillation at 164. With the budget raised to 1,000 the
  untreated residual is still 0.08623068394761953, while the intervention still
  converges in 39 steps. This is not just an insufficient initial budget.
- **slab:** default Pulay stagnates at iterations 65–70, range ratio
  1.0502702077546342 and endpoint ratio 0.9785361468857847. The 1,000-step fresh
  Pulay β=0.1 reference ends at residual 0.01331797371169278 without convergence.
  Seven baseline policies all exhaust at 200; no reference energy is available.
  Default Pulay also exhausts at 1,000 with residual 0.05528919253889995;
  its intervention arm still exhausts at 1,000 with residual 0.05126231476584983.
- **slab-offset:** Linear β=1 first stagnates at 1–6 (range/endpoint ratio
  1.0143025749675054). It still has residual 0.11227601418830685 after 1,000
  untreated steps; the intervention still converges in 31. Linear β=0.5/Hybrid
  first oscillates at 1–7: all five adjacent energy-difference signs alternate,
  every difference exceeds 1e-8 Ha, and residual ratio is 1.3551901113491622.

[Extended skewed table](../../tests/reference/stagnation/extended-skewed/table.md),
[extended slab table](../../tests/reference/stagnation/extended-slab/table.md)
and [extended slab-offset table](../../tests/reference/stagnation/extended-slab-offset/table.md)
contain the nine additional budget-control runs. These extend the same cold-start
baselines; they are not live budget interventions or checkpoint continuations.

No P11D divergence/nonfinite trigger occurred in this matrix. Linear β=0.1
usually needs about 115–131 steps even for simple convergent models; this is slow
convergence, not proof of charge sloshing. No residual Fourier probe was run:
`charge_sloshing_hint=not_computed` throughout. Neither vacuum nor elongation is
used as a substitute for measuring low-q residual dominance.

## Small-gap and metallic candidates

The converged `soft-mesh` reference has sampled gap 0.004163613860241472 Ha.
The `near-metal` and `skewed-mesh` reference attempts reach small density residuals
but fail the engine's `NonInsulating` gate, with sampled indirect gaps
-0.08669121378531516 Ha and -0.05484133664243912 Ha. The required gap is 1e-6 Ha.
These typed errors are exported separately as `terminal_evidence`; they are not
fabricated P11D diagnoses. No smearing or partial occupation was introduced.

Two interventions on `near-metal` remove persistent residual oscillations but
end at `NonInsulating`, so neither is validated recovery. Failed energies and
reference energy differences remain null. Legacy Gamma does not enforce the
multi-k insulating gate: its reported HOMO/LUMO gap is observational, not an
additional acceptance test. Gaps on exhausted iterates are explicitly marked
as not from a converged state. No ground-state or material-metallicity claim
is made from these synthetic models.

## Evidence, verification and limits

- [Complete main table](../../tests/reference/stagnation/measured/table.md): all
  100 runs, including every non-failing case, final residual/energy, triggered
  rules, gap and recovery comparison. [Aggregate JSON](../../tests/reference/stagnation/measured/results.json)
  contains complete histories, source/binary/input/UPF hashes, compiler/platform,
  expanded numerical options and audit. Per-run input/output files permit direct
  inspection. `traces.csv` contains every recorded scalar history with explicit
  units and is readable by ordinary plotting tools.
- [Verification record](../../tests/reference/stagnation/verification.json):
  100 identical normalized repeats, 91 checked source identities and 5 validated
  interventions. The second full raw bundle was written to a temporary directory;
  its normalized per-run SHA-256 values are retained. Regenerate it with the full
  command above and use `validate.py measured --repeat /tmp/stagnation-all`.
- `validate.py` also checks source/input hashes, individual-versus-aggregate
  equality, certificates on converged runs, P11D predicates against numeric
  windows, revision boundaries, identical pre-intervention histories, audit
  outcomes and energy comparisons. It intentionally fails if source identities
  change; a new build needs a new measurement bundle.
- `scripts/check.sh`: build/test/fmt/clippy all passed (exit 0). Separate harness
  fmt/clippy and evidence checks passed. `scripts/make_readme_figures.py` ran
  successfully and produced no tracked changes: this ticket does not own that
  generator or README, and the existing generator has no stagnation result type.
  README integration is therefore deferred to its owner, rather than inserting
  hand-maintained numerical text. See the [check record](../../tests/reference/stagnation/checks.json).

The audit is an in-memory experimental audit exported after execution, not a
crash-safe journal. No live agent, three-arm P11E integration, cross-build
bitwise identity, array equivalence, real-element validation, exhaustive search,
physical charge-sloshing diagnosis or production-scale band-count claim is made.
A recovered energy agrees with one same-engine reference, not a proof of the
unique/global minimum. The symmetric slab has no recovery under the tested
policy. Defaults remain unchanged.

Repository coordination: this worktree has no P11M running entry in its local
`artifacts/sessions.json` or `.delegate/P11M-STAGNATION-CASES/report.json`.
No replacement agent was spawned. Dashboard/session files and README-generation
code are outside this ticket's owned paths; orchestrator integration must update
them. No dashboard upload or recurring reporting timer was installed here.
Changes remain uncommitted on `task/P11M-STAGNATION-CASES` as requested.
