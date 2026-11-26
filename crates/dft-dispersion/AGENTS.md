# dft-dispersion

**Purpose:** Grimme DFT-D3 dispersion (rational and zero damping) through the official simple-dftd3 C API, opened at run time with libloading ([ADR 0048](../../docs/decisions/0048-dftd3-dispersion.md)). Energy (Ha), gradient (Ha/Bohr) and virial (Ha) of one periodic structure in Bohr. No physics is re-implemented here.

**Invariants to preserve and test:** without feature `dftd3` nothing is loaded and `compute` returns `FeatureDisabled`; `unsafe` lives only in `src/ffi.rs`, every block with a `SAFETY:` comment; only the C functions actually used are declared; library errors are surfaced verbatim; lattice rows = lattice vectors at this boundary (callers with column matrices transpose); stress = virial / volume (tensile positive) and force = -gradient, both established by `tests/finite_differences.rs`; partial parameter overrides are rejected (the C API cannot read back a loaded set). Never call this inside an SCF loop: once per geometry.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-dispersion`, then with the library (`python3 -I scripts/fetch_dftd3.py`): `BAANDI_DFTD3_LIB=<dest>/lib/libs-dftd3.dylib cargo test -p dft-dispersion --features dftd3 -- --include-ignored --nocapture` (oracle against the `s-dftd3` program of the same package; finite differences, about 40 s, ignored by default). Tests skip with a message when the variable is unset.
**Rules:** keep public API changes out of unrelated tickets; never depend on the LGPL `dftd3` Rust crate or link the library; never commit library binaries; every behavior change adds or updates a test with stated tolerance and units; record design decisions in an ADR.
