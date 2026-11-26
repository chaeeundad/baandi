# ADR 0048: Grimme DFT-D3 dispersion through simple-dftd3 loaded at run time

- Status: **Accepted** (orchestrator review 2026-11-26: ISC licence of libloading accepted; three-body term on by default as in s-dftd3 and QE) (implementation of ticket D3-DISPERSION; user decision
  2026-11-26, option A: official library behind a Cargo feature).
- Scope: D3 with rational (Becke-Johnson) and zero damping, added to energy,
  forces, stress and relaxation. Not in scope: D3M(BJ), D4, many-body
  dispersion, nonlocal vdW-DF/rVV10, D2, gCP, LDA parameter sets; the QE
  `vdw_corr = 'dft-d3'` comparison is run afterwards on dev08.

## Decision

1. **Library, not a re-implementation.** D3 is a geometry-only post-SCF term
   whose reference data (C6 tables, radii, parameter sets per functional) are
   large and easy to get subtly wrong. BAANDI uses the official simple-dftd3
   implementation (Grimme group, LGPL-3.0-or-later), consistent with ADR 0003
   (prefer maintained general-purpose libraries for numerics).
2. **Run-time loading, own adapter.** The library is opened with `dlopen`
   (`libloading`, ISC) by the new MIT crate `dft-dispersion`, which declares
   only the C functions it uses from the public header `dftd3.h` (version,
   error handle, structure, D3 model, damping load/new, get_dispersion,
   deletes). It does not link the library and does not depend on the LGPL
   `dftd3` Rust crate, so BAANDI's own code and binaries stay MIT; the LGPL
   library is a separate user-supplied file (fetched by
   `scripts/fetch_dftd3.py` from conda-forge with SHA-256 pins, never committed
   or redistributed by this repository). `unsafe` is confined to
   `crates/dft-dispersion/src/ffi.rs`, each block with its preconditions.
3. **Feature flag.** Everything is behind the Cargo feature `dftd3`
   (`dft-cli/dftd3` forwards to `dft-dispersion/dftd3`). The default build
   compiles neither the adapter nor libloading and rejects `[dispersion]` at
   preflight; without the table no new code path runs in any build (the helper
   returns before touching anything), so default results are unchanged.
4. **Once per geometry, after the SCF.** The term does not enter the Kohn-Sham
   Hamiltonian; it is evaluated on the final structure of a run and of every
   relaxation frame, never inside the SCF loop. E_disp is added to
   `total_energy`, `energies.total` and, on spin/Fermi-Dirac runs, to the
   internal and free energies; forces and stress add the library's analytic
   derivatives; relaxation frames carry the corrected energy, forces and stress
   so the optimizer follows the corrected surface.
5. **Parameters.** No damping parameter given: the library's set for the run's
   functional (PBE), with the three-body ATM term on by default (`s9 = 1`, as
   QE's `dftd3_threebody` default and the library's Python API); `s9 = 0` turns
   it off. All core parameters given: an explicit damping (`alpha6` default 14,
   zero-damping `rs8` default 1). Partial sets are rejected because the C API
   cannot return a loaded set to be modified. LDA has no D3 parameter set, so
   LDA requires explicit parameters.
6. **Units and signs (measured, not assumed).** The library takes Bohr and
   returns Hartree, gradient dE/dR (Ha/Bohr) and the virial dE/d eps (Ha) for
   the homogeneous deformation r -> (1 + eps) r. Central finite differences
   (h = 1e-4 Bohr, strain 1e-5) fix the conventions: force = -gradient and the
   engine's tensile-positive stress = **+virial / Omega** (the ticket's
   provisional "-sigma/Omega" was not used; the FD test would fail by
   2|stress| with the opposite sign). Our lattice matrix has the vectors as
   columns; the adapter boundary passes them as rows.

## Validation (Mac, 2026-11-26, simple-dftd3 1.7.0 conda-forge osx-arm64)

- Adapter against the independent `s-dftd3` program of the same package
  (Turbomole coord input in Bohr, JSON output), graphite AB, a triclinic
  Si/O/B/N cell and H2O + CH4 in a 16 Bohr box; PBE-D3(BJ) and PBE-D3(0) with
  and without ATM, plus explicit BJ and zero parameter sets (18 cases):
  max |dE| 5.2e-18 Ha, max |d gradient| 8.7e-18 Ha/Bohr, max |d virial|
  9.0e-17 Ha (acceptance 1e-10, 1e-9, 1e-9).
  `crates/dft-dispersion/tests/oracle_sdftd3.rs`.
- Finite differences of E_disp, same structures, BJ and zero damping:
  max force error 8.4e-11 Ha/Bohr, max stress error 2.3e-12 Ha/Bohr^3
  (acceptance 1e-6 both). `crates/dft-dispersion/tests/finite_differences.rs`
  (ignored in the gate, about 40 s).
- CLI end to end (`crates/dft-cli/tests/dispersion_cli.rs`), each against the
  same input without `[dispersion]`: the SCF energy is bit-identical, every
  total equals SCF + E_disp bit for bit, forces and stress differ from the
  plain run by exactly the D3 term (max difference 0). Synthetic LDA carbon
  cell (explicit BJ parameters, scalar path) and real PBE graphite AB with the
  PseudoDojo nc-sr-04 C file, 40 Ry, 2x2x1 k, Fermi-Dirac 0.01 Ha (reference
  spin/FD path): E_disp(BJ) = -2.1608568683e-2 Ha, E_disp(0) =
  -1.2722007644e-2 Ha. One relaxation step: frame 0 differs from the plain
  relaxation by E_disp (5.6e-17 Ha) and F_disp (7e-18 Ha/Bohr), and the step
  moves differently (2.7e-4 Bohr), so the optimizer used the corrected forces.

## Consequences and limits

- Users must obtain the library themselves (`scripts/fetch_dftd3.py`); only
  osx-arm64 is pinned so far. dev08 needs a linux-64 pin table.
- The library's default real-space cutoffs are used (no Ewald summation option
  exposed); this matches the `s-dftd3` program defaults, not necessarily QE's
  internal D3 cutoffs, which the dev08 comparison must state.
- The strict-oracle spin/Fermi-Dirac drivers (`BAANDI_REFERENCE_STRICT=on` with
  spin or smearing) reject `[dispersion]`; the default reference path and the
  scalar strict path carry it.
- Element data are accepted for H to Pu (Z <= 94); the atomic number comes from
  the UPF header element.
