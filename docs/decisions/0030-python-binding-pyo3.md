# ADR 0030: Python binding through PyO3, outside the core

Status: accepted (architect, 2026-11-04)

## Context
Directive D1 asks for an ASE calculator and a Python workflow API. Python 3.9 or later is available on the reference Mac; ASE 3.26 and maturin 1.16 install into a virtual environment.

## Decision
- A separate crate `baandi-py` (module `_baandi`) built with PyO3 0.30.1 and packaged with maturin; the pure-Python package lives in `python/baandi/` and imports `_baandi`.
- Core crates (`dft-*`) never depend on `baandi-py`; the engine API is used as published (stepwise session, events, policy, checkpoint).
- ASE is an optional Python dependency. Forces and stress raise `NotImplementedError` (never zero) until implemented.
- The GIL is released during long numerical calls; the numerical determinism contract is unchanged (same build, threads, seed, control log).
- Python never runs inside the SCF inner loop; control requests are applied at safe points through the same validated policy path as the MCP server and CLI.

## Consequences
New dependency `pyo3` (MIT OR Apache-2.0, review in `docs/dependencies.md`). Wheel distribution and CI for other platforms are not covered here.
