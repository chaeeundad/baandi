# dft-basis

**Purpose:** Plane-wave basis: G-vectors with 0.5*|G|^2 <= Ecut, FFT grid choice, G-index <-> FFT-index mapping.

**Invariants to preserve and test:** Negative-frequency wrap; cutoff boundary; round trip; Parseval; aliasing policy documented.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-basis`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Frozen-Miller deformation (P22A)

`PlaneWaveBasis::in_cell` preserves exact Miller order and fractional k, updating
only the reciprocal metric. Cutoff metadata is a reference value, not a new
sphere assertion. The engine freezes the density support separately; calling
`density_basis` still builds its documented production cutoff sphere. Existing
constructors remain unchanged. Strain oracles must freeze FFT dimensions too.

P23A real-file follow-up: `in_cell_with_frozen_density` explicitly carries the
reference cell for later `density_basis` requests, including Hamiltonian setup.
The engine selects this opt-in path for `ScfOptions.strain_reference`; plain
`in_cell` and all ordinary basis constructors retain their previous semantics.
The shell-crossing regression rejects the old rebuilt-density/grid combination
and checks exact frozen Miller ordering and changed reciprocal metrics.

## P26A real Gamma foundation

`half.rs` owns the library-only Gamma mirror table and ORC conversion. Zero is
first; positive Miller representatives preserve full-list order. Construction
checks inversion symmetry and the FFT map. Pack rejects nonfinite values and
absolute Hermitian defects above `ORC_TOLERANCE` (1e-14); unpack completes the
negative representatives. Neither operation allocates. No engine uses this API.

R1 local validation (2026-11-06): `cargo test -p dft-basis --lib --
--nocapture` exercises cubic/skewed/elongated/slab bases of 27/23/19/75
plane waves. Maximum dimensionless ORC norm/dot defect: 3.3306690738754696e-16
(gate 1e-15). Broken inversion symmetry, non-Gamma, wrong FFT map, nonfinite
coefficients, broken conjugacy and incorrect lengths return errors.
