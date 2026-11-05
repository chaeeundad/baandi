# baandi-py

Purpose: optional PyO3 binding to dft-engine and Python/ASE integration (ADR 0030).
Read root AGENTS.md, ownership.toml and ADR 0030 before edits. No core crate depends
on Python. Unit conversion uses dft-units; no duplicated physical constants.
Python receives owned results and drives the existing session; do not duplicate
SCF or policy validation. ASE energy requires certified convergence. Forces,
stress and free energy remain unsupported, never synthesized as zeros.

Verify: `scripts/check.sh baandi-py`, then `scripts/check.sh`. Native wheel and ASE
commands/evidence are in `VALIDATION.md`, `native-evidence.json` and
`ase-evidence.json`. Tests that skip when the extension is absent do not establish
native acceptance. No MLIP loop or low-level external coupling is implemented.
