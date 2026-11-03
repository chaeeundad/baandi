# ESLE roadmap aligned with the repository

Status: plan; follows the directive in [../research/esle-directive.md](../research/esle-directive.md). Dates are not promised; each phase has entry and exit gates. No ESLE code is written before Phase A is complete.

| Phase | Content | Entry gate | Exit gate (evidence) | Notes |
|---|---|---|---|---|
| A (now) | k-point validation against QE, baseline CPU KS-DFT for simple periodic systems | P7A and P7A2 merged (done) | QE k-point comparison recorded in `docs/validation.md` (R6) and `tests/reference`; P7E2 expression wiring merged | no ESLE feature may delay it; QE runs 2x2x2, shifted 2x2x2, shifted 4x4x4 delivered, assembly in progress |
| B | State infrastructure: `run_from`, state/compatibility modules in `dft-engine`, persistence in `dft-cli`, save/load, reproducibility tests | A complete | save then load is bit-identical; restart from a saved density reproduces the stored state within measured noise; noise floors for each metric measured and documented | additive API only, no new crates (ADR 0016) |
| C | Multi-start SCF and state comparison for the non-spin insulating model; registry that records failures | B complete | experiments with uniform, perturbed and restarted starts on Si; distance metrics with measured tolerances; three-way classification tested on constructed cases | UPF `PP_RHOATOM` exposure needed for atomic superposition |
| D | Density transfer and tracking across small distortions | C complete | normalization, electron count and grid transformation tests; continuation versus fresh solve compared (Experiments A to C) | Si lattice deformation |
| E | Spin, smearing, DFT+U, occupation control | explicit user approval and ADRs | each method validated against QE or analytic cases before any ESLE use | out of scope today (AGENTS.md); Al needs smearing, MgO needs a second verified pseudopotential |
| F | Adaptive exploration | E for the interesting systems, or a benchmark with known multiple solutions | beats fixed multi-start on a stated benchmark | any learned policy needs the deterministic baseline |

## Dependencies the directive implies but the code does not have
- Aluminum (metal): needs smearing and partial occupations (now rejected by design) plus k-point convergence study; a new validated feature, so it is a separate ticket with its own QE comparison before it becomes an ESLE test system.
- MgO: needs a verifiable O and Mg norm-conserving LDA UPF (provenance first, see `docs/mvp-review.md`) and a check that the UPF reader handles those files; the Si file exposed three reader defects, so expect more.
- Magnetic and correlated systems: spin-polarized DFT and DFT+U do not exist; Phase E.

## Sequence after the Phase A gate (from the independent review)
(a) model and basis compatibility spec plus an additive-initialization ADR; (b) density-only `run_from` and the evaluated-state envelope with default-behavior regressions; (c) immutable density-state persistence, strict reload, all failures recorded; (d) optional orbital persistence and restart, projector comparison with calibration; (e) deterministic multi-start with conservative grouping; (f) explicit density transfer and branch tracking against independent starts. Atomic superposition joins after its radial-density adapter is validated. Each step is independently reviewable and needs no Hamiltonian redesign.

## Integration milestones
1. M-B1: `ScfProblem::run_from` and `ScfOutcome` extraction with no behavior change (regression: Si Gamma density bit-identical).
2. M-B2: state save/load (manifest committed after payload writes, fail-closed on unknown schema or digest mismatch) and validation tests.
3. M-B3: noise-floor measurements (seeds, thread order, tolerance schedule).
4. M-C1: distance metrics and classification tests on constructed states.
5. M-C2: multi-start experiment report for Si.
6. M-D1: density transfer for small deformations with the Experiment A to C reports.

## Independent review for every ESLE ticket
Checklist for the reviewer (Codex or a separate Claude agent, not the implementer): electron conservation after every transfer or restart; occupation and spin consistency; basis compatibility (cell, cutoff, grid, k set, ordering); energy comparisons only inside one model fingerprint; correct reading of convergence (residual, not energy change alone); no merge of states from similarity of energy alone; reproducibility (repeat run, different seed and thread layout); failures recorded, not dropped. A self-reported success never replaces re-running the gates.
