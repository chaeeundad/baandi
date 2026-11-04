TASK_ID: P11G-PYTHON-ASE
GOAL: PyO3 binding and ASE calculator (directive D1) exposing energy, density summary, eigenvalues and stepwise control; separate crate/package outside the core; maturin build; ASE is an optional Python dependency. Forces/stress are NOT available yet and must raise a clear NotImplemented error (never zero).
OWNED_FILES: crates/baandi-py/ (new), python/baandi/ (new), docs/architecture/scientific-integration.md (usage section)
DEPENDENCIES: P11A merged; ADR for PyO3 (license, Apple Silicon check, unsafe policy in FFI adapter); network needed for crates/pip: Claude fetches.
ACCEPTANCE_CRITERIA: ASE calculator energy equals CLI result bit-for-bit for the same input; unsupported properties raise; wheel builds on Apple Silicon.
GIT: branch task/P11G-PYTHON-ASE, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
ARCHITECT DECISION 2026-11-04: PyO3 0.30.1 approved (ADR 0030, docs/dependencies.md). Crate skeleton `crates/baandi-py` is registered in Cargo.toml/ownership.toml; the local virtual environment is /Users/mino/Develop/rust-dft-wt/pyenv (maturin 1.16.0, ASE 3.26.0, numpy); use PYO3_PYTHON=/Users/mino/Develop/rust-dft-wt/pyenv/bin/python for builds. The sandbox may be unable to run maturin develop: write Python tests that skip cleanly when the extension is not importable and report what ran. Use at least two different systems in tests.

