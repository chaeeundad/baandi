# baandi-py

Purpose: optional PyO3 binding to dft-engine and Python/ASE integration (ADR 0030).
Read root AGENTS.md, ownership.toml and ADR 0030 before edits. No core crate depends
on Python. Unit conversion uses dft-units; no duplicated physical constants.
Python receives owned results and drives the existing session; do not duplicate
SCF or policy validation. ASE energy requires certified convergence. Free energy remains unsupported, never synthesized as zero. P22A adds certified frozen-Miller stress for the existing Gamma/unpolarized session. P20A adds
post-convergence forces for the existing Gamma, unpolarized NC model; the worker
uses retained model inputs and requires the engine certificate. ASE converts
Ha/Bohr to eV/Angstrom using dft-units conversions.

Verify: `scripts/check.sh baandi-py`, then `scripts/check.sh`. Native wheel and ASE
commands/evidence are in `VALIDATION.md`, `native-evidence.json` and
`ase-evidence.json`. Tests that skip when the extension is absent do not establish
native acceptance. No MLIP loop or low-level external coupling is implemented.

ASE stress is tensile-positive in eV/Angstrom^3, Voigt xx,yy,zz,yz,xz,xy.
The native method returns Ha/Bohr^3 and warns on the separate Pulay assessment.
`test_ase_stress_sign_voigt_and_energy_derivative` uses ASE's own numerical
stress oracle. Do not claim installed-wheel acceptance from skipped tests.

Stage R routing ([ADR 0047](../../docs/decisions/0047-reference-driver-single-path.md)):
`Options.driver = "reference"` is the only accepted driver spelling. Default
`reference_strict = False` uses the reference facade/runtime; `True` selects
strict-oracle compatibility (Gamma only in this adapter). `"legacy"` was removed
in R10/R11, see git history. Reference defaults are Broyden beta 0.7/history 8,
atomic start, threshold 5e-9 Ha and verification on. Strict defaults remain
Pulay beta 0.3/history 6. Forces/stress use retained delivery and deferred
observables, not repeated SCF. Native tests must assert the selected route.
Historical pre-R5 evidence retains its original driver.

R8d: initialize/step/status/iteration drive `refscf::runtime::Runtime` around
`api::stepper`, on the existing owning worker thread. `run()` initializes a new
runtime or completes its live continuation. Status uses reference RuntimeState
names (Created/Initialized/Converged/Exhausted/Failed); strict oracle retains StepStatus.
Reference step returns an owned, immutable StepReport object with every engine
field, including map kind and safe_point. Verification/retry maps can repeat an
iteration number. The terminal delivery is retained for result and derivatives;
no SCF replay or legacy fallback. Explicit strict-oracle step still returns its status
string for compatibility. Stress still requires threshold <= 5e-11 Ha selected before the run;
the separate Pulay assessment warning is preserved. Keep all pre-R5 evidence as
historical legacy evidence. Native acceptance lives in `tests/test_reference.py`
(pytest-compatible unittest tests); ordinary Rust test builds still omit PyO3.

`Options.kgrid` defaults to `[1, 1, 1]` (Gamma); other positive meshes use the
engine's unshifted Monkhorst-Pack/time-reversal helpers through the approved
`dft-engine::input::kpoints` re-export. Only Gamma uses GammaReal. Explicit strict-oracle
sessions reject non-Gamma meshes. New Si2 Gamma/k222 native/CLI bit-identical
energy, force and stress evidence is in `tests/r5-reference-evidence.json`.
