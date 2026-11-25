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

Historical ADR 0043 ONEPATH-4 defaults (2026-11-16; legacy after R5): `Options()` defaults to the standard path with
the fixed-point verification (`standard_threshold_ha` 5e-9, `verify` True,
`mixing_beta` 0.5); `reference_strict` (default mirrors `BAANDI_REFERENCE_STRICT`)
selects the strict reference. Forces need a passed verification or the strict
certificate; `stress()` refuses a threshold looser than
`StandardProfile::STRESS_THRESHOLD_HA` (5e-11 Ha). The native/CLI bit-equality
tests set `cutoff_ladder = "off"` on the CLI side (the session has no ladder).


R5 reference batch routing (2026-11-25): `Options.driver = "reference"` is
now the default. Only the explicit `"legacy"` driver constructs `ScfSession`.
Reference `run()` retains the reference facade delivery (R8d uses its runtime/stepper); result/density access reads
that output, and forces/stress use its deferred `Observables` in prepare,
forces, stress order. Repeated access does not run SCF again. No fallback on
engine errors. Defaults are reference Broyden beta 0.7/history 8, atomic density
and orbitals, standard threshold 5e-9 Ha, verification on. Legacy defaults remain
Pulay beta 0.5/history 6 (strict beta 0.3); explicit mixing overrides survive a
driver change. `reference_strict` requires the explicit legacy driver.

R8d: initialize/step/status/iteration drive `refscf::runtime::Runtime` around
`api::stepper`, on the existing owning worker thread. `run()` initializes a new
runtime or completes its live continuation. Status uses reference RuntimeState
names (Created/Initialized/Converged/Exhausted/Failed); legacy retains StepStatus.
Reference step returns an owned, immutable StepReport object with every engine
field, including map kind and safe_point. Verification/retry maps can repeat an
iteration number. The terminal delivery is retained for result and derivatives;
no SCF replay or legacy fallback. Explicit legacy step still returns its status
string for compatibility. Stress still requires threshold <= 5e-11 Ha selected before the run;
the separate Pulay assessment warning is preserved. Keep all pre-R5 evidence as
historical legacy evidence. Native acceptance lives in `tests/test_reference.py`
(pytest-compatible unittest tests); ordinary Rust test builds still omit PyO3.

`Options.kgrid` defaults to `[1, 1, 1]` (Gamma); other positive meshes use the
engine's unshifted Monkhorst-Pack/time-reversal helpers through the approved
`dft-engine::input::kpoints` re-export. Only Gamma uses GammaReal. Explicit legacy
sessions reject non-Gamma meshes. New Si2 Gamma/k222 native/CLI bit-identical
energy, force and stress evidence is in `tests/r5-reference-evidence.json`.
