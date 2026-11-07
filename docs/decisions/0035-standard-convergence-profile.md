# ADR 0035: Additive "standard" convergence profile for QE-comparable runs

- Status: **accepted by the user 2026-11-05** (decision U-perf-1 of docs/performance/exceed-qe-plan.md; the default behaviour is NOT changed).
- Context: at conv_thr 1e-8 QE declares convergence from a Hartree-metric error estimate below 1e-8 Ry, tests eigenvector convergence by the change of the Ritz value, tolerates up to five unconverged roots in SCF and runs no certificate. BAANDI requires density L2 <= 1e-7, dE <= 1e-8 Ha, residual <= 1e-8 and a 1e-9 certificate solve (docs/research/qe-study/qe-scf-cycle.md and qe-solvers-and-subspace.md, tests/reference/perf/p15a). A fair "equal accuracy" speed comparison needs a profile that is comparable to QE's definition and a measurement of delivered accuracy against each code's own tight reference.
- Decision: add an OPTIONAL, explicitly labelled convergence profile `[scf] convergence_profile = "strict" | "standard"`; `strict` is the current behaviour and stays the default; `standard` is defined by a documented set of criteria written in BAANDI's own terms (Hartree-weighted error estimate with its ethr ratchet, Ritz-change stopping rule, bounded count of unconverged buffer roots) and is introduced only after the spectral completeness guard of ticket P15C is merged. Every result produced under `standard` carries `convergence_profile: "standard"`, `certified: false` and the profile parameters in results.json, in the CLI output and in the checkpoint provenance; a `standard` result is never reported as certified, never used in validation claims against other codes, and never mixed with strict results without the label.
- Rules: no tolerance of the strict path is changed; the profile is chosen explicitly by the user; speed comparisons with QE must state which profile and must be accompanied by the delivered accuracy (final energy and eigenvalue and density differences against each code's own tight reference) before any "faster than QE" claim; the numeric parameters are set after experiment E1 (qe-scf-cycle.md) and frozen before the comparison.
- Consequences: new input key and results fields (additive schema), tests that strict results are bit-identical with and without the feature present, and documentation of the equal-accuracy protocol.

## Implementation status (2026-11-05, P17A)

Implemented as an additive unpolarized Davidson profile. E1 freezes the default
at 5e-9 Ha using the architect-approved stored QE error estimates (the rounded
1e-8/1e-12 energy difference is unresolved, not an exact zero) and six independent
BAANDI tight-reference plateaus, including Si/MgO and synthetic geometries.
All occupied and first empty roots are mandatory; higher buffers may remain
unconverged only at the solver limit. No standard result is certified.
The independent Hartree metric remains opt-in: E4 did not meet a uniform
1e-12 Ha energy-equivalence gate on MgO. Three stress fixtures have no tight
reference and retain work-only reporting. No dev08 speed acceptance is claimed.
See [implementation, evidence and limits](../performance/p17a-standard-profile.md).
This status note does not alter the accepted decision.

## Amendment: experimental standard stall guard (P29A, 2026-11-07)

- Status: **Proposed**. The accepted profile and strict default above are unchanged.
- Add `[scf] standard_stall_guard = "off" | "on"`, default `off`; `on` requires
  the explicitly selected standard profile. Results remain uncertified. Persist
  the enabled setting in result/checkpoint provenance; exact continuation rejects
  a different setting. Omitted legacy wire fields mean `off`.
- Inspect the last six completed standard-error records (five intervals). If
  every newer error exceeds half the entry error, cap the next eigen tolerance
  at 0.1 times the **window-entry** requested tolerance, retaining the existing
  ordinary ratchet and 5e-14 Ha floor. Using the entry tolerance avoids repeatedly
  compounding a cut while the same plateau occupies the window. No mixer changes,
  mutable controller state, allocations or I/O are added to the numerical loop.
- Evidence: ten perturbations on each of ten benchmark/material cases reproduce
  the MKL-like pattern in 8/10 Si8 and 2/10 Si16 trajectories. Eight cases have
  two strict references agreeing within 1e-12 Ha. With the guard, all 80 formerly
  convergent trajectories still converge, and the 60 nonstagnating histories are
  identical. Si8's maximum SCF count falls 32→27, but median H work rises 0.99%.
- Negative result: neither tested factor (0.5, 0.1) passes the prespecified
  work-saving selection rule. Held-out Si16 gains no SCF-count improvement;
  elongated/slab remain nonconvergent and median H work rises 128.51%/141.33%.
  Whole-set H work rises 49.61%. The switch is an experimental tail-robustness
  option, **not recommended for general use or default promotion**, and is not
  an accepted speed optimization. No parameters were retuned on held-out cases.
- Accuracy: the largest measured absolute energy difference against the same
  tight references is 3.164040e-9 Ha with the guard, versus 2.245319e-9 Ha without.
  The standard stopping estimate is not a proof of total-energy accuracy. No
  delivered-accuracy claim is possible for the two failed-reference stress cases.
- Full protocol, distributions, two histories, held-out split, reproducible raw
  evidence and limitations: [standard stall study](../performance/standard-stall-study.md).

## Amendment: robust gain evidence and density-forcing experiment (P29B, 2026-11-07)

- Status: **Proposed**. The accepted standard profile and strict default remain
  unchanged. This amendment does not promote a research rule.
- Replace any interpretation of the Si8 historical single-run 2.65x H / 3.0x
  dev08 wall figures as robust guarantees with explicit ensemble statistics.
  Local ten-seed paired strict/standard total-H median ratio is 1.7849x
  (95% paired bootstrap 1.6503–2.0044x); worst observed pair 1.4036x.
  Mandatory strict certificates are included; main-only ratios are separately
  reported. No robust dev08 wall measurement exists in this study.
- One opt-in research key, `[scf] standard_density_forcing = "off" | "on"`,
  defaults off. On requires standard and `standard_stall_guard=off`; unknown
  values and invalid combinations reject. Off serialization stays omitted;
  enabled result/options/checkpoint/inspection provenance is explicit and exact
  restart rejects a changed setting. Standard results stay uncertified.
- On applies `max(5e-14 Ha, min(ordinary_ratchet,
  1 Ha * (previous_density_L2 / (1 Bohr^-3/2))^2))`. The first-step rule, ordinary
  floor, stopping threshold, strict tolerances and certificate are unchanged.
  The inexact-solve motivation and quadratic Rayleigh-error analogy are a
  heuristic, not a theorem for Pulay SCF or an energy bound. No mutable state,
  allocation or I/O is introduced in the numerical loop.
- Training: four cases, seeds 1–5, H 8784→7359 (-16.22%), all 20 converged.
  Freeze coefficient 1 before held-out evaluation; the conditional second
  parameter was not tried. No held-out case/seed retuning.
- Negative generalization: ten cases, held-out seeds 6–10, H 139332→146296
  (+5.00%); whole-case holdouts +6.31%, successful-only +0.75%. Si8/Si16 improve
  but k222s/Si-PBE/MgO-PBE regress; elongated/slab still fail. All 40 previously
  successful held-out trajectories converge, max energy error 8.633592e-10 Ha.
  **No tested rule beats standard on the declared held-out aggregate criterion;
  leave this option experimental, off, and not recommended for general use.**
- Full per-case distributions, accuracy, statistical definitions, failure
  accounting, exact commands, gates and dev08 request:
  [standard robust gain](../performance/standard-robust-gain.md).
