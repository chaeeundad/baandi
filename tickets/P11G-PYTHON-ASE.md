TASK_ID: P11G-PYTHON-ASE
GOAL: PyO3 binding and ASE calculator (directive D1) exposing energy, density summary, eigenvalues and stepwise control; separate crate/package outside the core; maturin build; ASE is an optional Python dependency. Forces/stress are NOT available yet and must raise a clear NotImplemented error (never zero).
OWNED_FILES: crates/baandi-py/ (new), python/baandi/ (new), docs/architecture/scientific-integration.md (usage section)
DEPENDENCIES: P11A merged; ADR for PyO3 (license, Apple Silicon check, unsafe policy in FFI adapter); network needed for crates/pip: Claude fetches.
ACCEPTANCE_CRITERIA: ASE calculator energy equals CLI result bit-for-bit for the same input; unsupported properties raise; wheel builds on Apple Silicon.
GIT: branch task/P11G-PYTHON-ASE, uncommitted
RETURN_FORMAT: docs/schemas/report.schema.json
