TASK_ID: PRELAX-PRECON
GOAL: Faster and more reliable ionic relaxation (fewer DFT energy/force evaluations, lower end-to-end wall) through physics-based preconditioning, without compromising accuracy or reliability. User specification of 2026-11-27 (full text: docs/roadmap/precon-relax-spec.md). The algorithm is not the objective; measured end-to-end improvement is.

PHASES (gated in order; each phase ends with a short report and evidence):
- Phase 0: audit of the current relaxation (BFGS, FIRE, force/SCF invocation, warm starts, Hessian init/update, step control, convergence criteria, constraints, PBC, restart/checkpoint, early/late bottlenecks, SCF noise). Reproducible baseline.
- Phase 1: DFT validation through the BAANDI ASE calculator: BAANDI BFGS, BAANDI FIRE, ASE BFGS, ASE LBFGS, ASE PreconLBFGS (Exp), experimental A3 metric + L-BFGS. Systems: Si8-AlP, 63-atom Si-Al-P vacancy, perturbed bulk Si, defects, then surfaces/adsorbates (TiO2 + molecule) only after a measurable benefit. Multiple seeds/perturbations. Fair conditions (section 6 of the spec). Target fmax = max atomic force NORM 0.02 eV/A (3.89e-4 Ha/Bohr) plus intermediate thresholds; common external evaluator; count every attempted evaluation.
- Phase 2 (only if Phase 1 shows a meaningful benefit): native Rust preconditioned optimizer (sparse P = mu(L + eps I) style, SPD on the active subspace, sparse solves, warm starts preserved, constraints/PBC, restart/checkpoint, safeguards of spec section 8), tests of spec section 12, benchmarks of section 9/10, QE comparison later.

OWNED_FILES (Phase 0/1): docs/roadmap/precon-relax-*.md, tests/reference/relax-precon/** (benchmark harness, inputs, raw results), python/ or scripts/ harness files for the ASE comparison (no runtime dependency of the engine on Python/ASE). Phase 2 ownership is set when Phase 2 is approved by evidence.

REFERENCES: Packwood et al., J. Chem. Phys. 144, 164109 (2016), doi:10.1063/1.4947024; ASE PreconLBFGS (LGPL-2.1: study and run as an external benchmark only, never copy code into the MIT tree); experimental files (outside the repo) /Volumes/990-2TB/Develop/rust-dft-wt/exp-tabpfn/tools/icl/{preconditioner.py,run_precon.py}, exp-tabpfn/cl/{stage2_emt.json,stage2_preconlbfgs.json}, exp-tabpfn/relax/ and exp-tabpfn/cl/{si8_AlP,...}.

RULES: AGENTS.md (no unrelated electronic-structure changes; fairness; same parallel model; Mac timings are observations, acceptance on dev08); evidence with commands, versions, inputs, UPF SHA-256; negative results documented.
GIT: dedicated branch per phase; leave uncommitted; the orchestrator commits.
