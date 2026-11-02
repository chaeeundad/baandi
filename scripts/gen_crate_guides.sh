#!/usr/bin/env bash
# Regenerates crates/<name>/AGENTS.md (purpose, invariants, verify command). Edit the table below, then run.
set -euo pipefail
cd "$(dirname "$0")/.."
gen() {  # name | purpose | invariants
  local f="crates/$1/AGENTS.md"; [[ -d "crates/$1" ]] || return 0
  cat > "$f" <<EOF
# $1

**Purpose:** $2

**Invariants to preserve and test:** $3

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** \`scripts/check.sh $1\`, then \`scripts/check.sh\` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.
EOF
}
gen dft-units    "Hartree atomic units and constants. The only place for unit conversions (1 Ha = 2 Ry); newtypes over bare f64." "Ry<->Ha and Bohr<->Angstrom round trips; no mixed-unit APIs."
gen dft-geometry "Cell (lattice columns a1..a3), atoms, fractional/Cartesian conversion, reciprocal lattice B = 2*pi*inverse(A)^T." "a_i . b_j = 2*pi*delta_ij; volume = |det A|; fractional<->Cartesian round trip."
gen dft-numerics "Adapters to general-purpose numerical libraries (linear algebra, FFT). Thin traits so backends stay swappable; the only crate allowed to contain unsafe." "Column-major complex layout stated; FFT normalization stated; plans and scratch reused; unsafe preconditions documented and tested."
gen dft-basis    "Plane-wave basis: G-vectors with 0.5*|G|^2 <= Ecut, FFT grid choice, G-index <-> FFT-index mapping." "Negative-frequency wrap; cutoff boundary; round trip; Parseval; aliasing policy documented."
gen dft-upf      "UPF v2 norm-conserving reading and validation layer, radial-to-reciprocal transforms. Never copy QE (GPL) code." "Unsupported files fail with explicit errors; no silent field skipping; Hermitian D_ij; finite, monotonic radial mesh."
gen dft-physics  "Kinetic, local ionic, nonlocal KB, Hartree, LDA XC, density; matrix-free H*psi." "Hartree G=0 convention; XC v = d(n*eps)/dn finite-difference test; H Hermitian; density integrates to N_valence."
gen dft-solvers  "Eigensolvers (dense baseline, then Davidson/LOBPCG), mixers (linear, then Pulay), convergence checks." "Residual and orthonormality reported; convergence needs density, energy and eigen residual, never energy alone."
gen dft-engine   "Orchestrates SCF, occupations, total energy and its audit (ion-ion/Ewald, double counting)." "Every energy term has a unit test and appears in results; no hidden constants."
gen dft-cli      "CLI, TOML input, output files. No numerical logic." "run/validate/inspect-upf subcommands; nonzero exit with cause on failure; outputs carry reproducibility metadata."
echo "generated: $(ls crates/*/AGENTS.md | wc -l) guides"
