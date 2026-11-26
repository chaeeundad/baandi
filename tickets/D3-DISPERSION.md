TASK_ID: D3-DISPERSION
GOAL: Grimme DFT-D3 dispersion correction (BJ and zero damping) added to energy, forces, stress and relaxation, computed by the official simple-dftd3 library loaded at run time behind the Cargo feature `dftd3` (user decision 2026-11-26, option A). Default builds stay unchanged and bit-identical.

OWNED_FILES:
- crates/dft-dispersion/** (new crate: the D3 adapter; the only place with `unsafe`/FFI for this task)
- crates/dft-cli/Cargo.toml (feature `dftd3` forwarding to the new crate, optional dependency)
- crates/dft-cli/src/config.rs, crates/dft-cli/src/preflight.rs (new `[dispersion]` table and checks)
- crates/dft-cli/src/run.rs (add the dispersion term after SCF to energy, forces, stress and the relax loop; JSON and table output)
- crates/dft-cli/tests/dispersion_*.rs
- Cargo.toml (workspace members), ownership.toml (new crate entry)
- docs/decisions/0048-dftd3-dispersion.md, docs/dependencies.md, docs/cli.md (`[dispersion]` section), docs/feature-status.json (dispersion entry)
- scripts/fetch_dftd3.py (fetches the conda-forge simple-dftd3 binary and its shared-library dependencies into a user directory; verifies SHA-256)

READ_ONLY_CONTEXT: AGENTS.md, ownership.toml, crates/dft-cli/AGENTS.md, docs/physics-conventions.md, docs/dependencies.md, ADR 0003/0004/0017, crates/dft-cli/src/run.rs (requested_forces, requested_stress, relax loop), the public simple-dftd3 C header and documentation (https://github.com/dftd3/simple-dftd3, `include/dftd3.h`, docs on the C API).

INPUT_CONTRACT:
```toml
[dispersion]
method = "d3bj"        # "d3bj" | "d3zero"  (later: "d3mbj", "d4" through another adapter)
# optional overrides; when absent the library's parameters for the run's functional are used
# s6, s8, a1, a2, rs6, alpha6, s9 (three-body, library default)
# library = "/path/to/libs-dftd3.dylib"   # else env BAANDI_DFTD3_LIB, else the platform default name
```
Functional name passed to the library: "pbe" for PBE, "lda" is not parameterized by D3, so `method` with LDA must be refused by preflight unless all damping parameters are given explicitly.

PHYSICS_CONVENTIONS: Hartree atomic units. The library takes Bohr positions and lattice (rows = lattice vectors in its API; our A has lattice vectors as columns, transpose at the boundary), returns energy in Hartree, gradient dE/dR in Ha/Bohr (force = minus gradient), and sigma = virial in Ha (stress = -sigma/Omega with our sign convention; check against docs/physics-conventions.md and against finite differences, do not assume). Periodic in all three directions.

IMPLEMENTATION_STEPS:
1. Thin own FFI adapter (libloading, dlopen at run time; MIT/Apache crates only, do NOT depend on the LGPL `dftd3` Rust crate) declaring only the C API functions used: version, error handle, structure create/update, damping load (rational/zero), param load/new, model create, get_dispersion, deletes. Preconditions documented next to each unsafe block; error strings surfaced.
2. Safe API: `compute(species_z, positions_bohr, lattice_bohr, method, functional, overrides) -> {energy_ha, gradient_ha_per_bohr[n][3], sigma_ha[3][3]}`.
3. CLI: parse `[dispersion]`, preflight (feature present, library loadable, functional parameterized), add E_disp to the total energy (reported separately as `energy.dispersion_ha` and included in total and free energy), forces and stress (reported as their own component too), relaxation uses the corrected forces/stress.
4. scripts/fetch_dftd3.py: download conda-forge simple-dftd3 1.7.0 for the host platform plus its runtime shared-library dependencies into a directory given on the command line (default under /Volumes/990-2TB/Develop/rust-dft-wt/vendor/dftd3), record URLs and SHA-256 in a manifest. Never commit binaries.
5. Docs: ADR 0048 (why library, licence reasoning: LGPL library loaded dynamically, own MIT adapter; feature flag; what is validated), dependencies.md entry (licence, maintenance, Apple Silicon), cli.md section, feature-status entry.

TEST_COMMANDS:
- scripts/check.sh (default features; must stay green and bit-identical: no change without `[dispersion]`)
- cargo test -p dft-dispersion --features dftd3 and cargo test -p dft-cli --features dftd3 --test 'dispersion_*' with BAANDI_DFTD3_LIB set; tests that need the library skip with a clear message when it is absent
- independent check against the `s-dftd3` command-line program from the same conda-forge package on the same structures (energy, gradient, virial)

ACCEPTANCE_CRITERIA:
- Adapter vs `s-dftd3` CLI on at least 3 structures (graphite AB, a skewed cell, a molecule-in-box with >= 2 elements), PBE-D3(BJ) and D3(0): |dE| <= 1e-10 Ha, max |d gradient| <= 1e-9 Ha/Bohr, max |d sigma| <= 1e-9 Ha.
- Central finite differences of E_disp (h = 1e-4 Bohr, strain 1e-5): forces within 1e-6 Ha/Bohr, stress within 1e-6 Ha/Bohr^3 equivalent.
- One end-to-end CLI run (small graphite cell, PBE) where total energy = SCF energy + E_disp exactly and forces/stress include the term; a relax step uses it.
- Default build: RG-ID style identity unchanged (no code path touched without `[dispersion]`); scripts/check.sh PASS.

NON_GOALS: D4, many-body MBD, nonlocal vdW-DF/rVV10, D2, gCP, LDA parameter sets, QE comparison (done afterwards by the orchestrator on dev08 with `vdw_corr='dft-d3'`).

DEPENDENCIES: libloading (MIT/Apache, needs a docs/dependencies.md entry). Runtime: simple-dftd3 1.7.0 shared library (LGPL-3.0-or-later, conda-forge build).

GIT: leave changes uncommitted in the worktree; the orchestrator commits after review.
RETURN_FORMAT: docs/schemas/report.schema.json
