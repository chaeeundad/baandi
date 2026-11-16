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

ADR 0043 ONEPATH-4 (2026-11-16): `Options()` defaults to the standard path with
the fixed-point verification (`standard_threshold_ha` 5e-9, `verify` True,
`mixing_beta` 0.5); `reference_strict` (default mirrors `BAANDI_REFERENCE_STRICT`)
selects the strict reference. Forces need a passed verification or the strict
certificate; `stress()` refuses a threshold looser than
`StandardProfile::STRESS_THRESHOLD_HA` (5e-11 Ha). The native/CLI bit-equality
tests set `cutoff_ladder = "off"` on the CLI side (the session has no ladder).
