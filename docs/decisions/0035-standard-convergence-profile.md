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
