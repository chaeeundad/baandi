# P21A fixed-cell relaxation evidence

The driver and MF07 density-transfer minimum are implemented in CLI/engine
setup code. There are no changes in the SCF iteration arithmetic, no new
numerical dependency, and no force-tolerance relaxation. Usage and restart/output
contracts are in [the user guide](../user/relaxation.md). Evidence is stored in
[`tests/reference/relaxation`](../../tests/reference/relaxation/). All numbers below
are BAANDI results or independent analytic/finite-difference checks, not QE data.

## Measured oracles

- Analytic positive quadratic and Rosenbrock-like objectives: both BFGS and FIRE
  recover the analytic minimum to coordinate error below 1e-6 in the test's
  dimensionless coordinates, with force stopping 1e-8 and step cap 0.1.
- Independent scalar energy scan: the synthetic periodic two-atom Gaussian
  local-potential fixture has a fitted [111] minimum distance
  6.928203229686338 Bohr. Three fine energy samples at fractional
  t=0.4999, 0.5, 0.5001 determine the parabolic vertex, without using forces.
  Starting at t=0.4 and t=0.6, capped BFGS converges in 12 steps each; distance
  errors are 6.0442e-10 and 5.7377e-10 Bohr, below the ticket's 1e-4 Bohr cap.
  This is a **periodic synthetic minimum**, not an isolated molecule or a
  hydrogen bond-length prediction. The two unwrapped starts are related by
  periodic inversion; they test opposite approach directions to the same basin.
- Along the first four nonzero trajectory directions from each start, central
  energy differences at 0.01, 0.005, 0.0025 Bohr show approximately fourfold
  error reduction. The largest finest-step error is 3.466e-9 Ha/Bohr,
  below the unchanged 1e-6 Ha/Bohr synthetic force cap.
- FIRE on the same displaced periodic fixture reaches maximum force
  3.923755e-8 Ha/Bohr in 190 steps with requested tolerance 1e-7 Ha/Bohr.
- Three displaced structures with tight SCF: warm/fresh energy differences
  2.220446e-16, 5.551115e-16 and 2.220446e-16 Ha (cap 1e-10 Ha).
  Invalid density values, charge, shape/layout, cell, species, PP content, XC, cutoff,
  occupation layout, k mesh and density-floor changes reject before initialization.
- Restart at a serialized accepted boundary produces byte-identical restart,
  trajectory, POSCAR, final TOML, log and results files to the uninterrupted run.
  The test also verifies an exactly fixed atom and exactly fixed Cartesian
  component. A one-iteration initial SCF budget returns exit 4 and no invented
  forces. Raw forces remain distinct from projected stopping forces.
- Auto-symmetry trajectory: all five structures retain 16 detected operations;
  forbidden Cartesian displacements remain below 1e-14 Bohr (arithmetic
  averaging residuals are not asserted to be exactly zero). Symmetry-breaking
  constraints reject. Standard-profile and invalid constraints/tolerances reject.
- ASE BFGS and CLI BFGS on an independently generated synthetic nonlocal fixture:
  maximum coordinate difference 3.619402e-7 Bohr; energy difference
  7.771562e-16 Ha; 7 versus 12 ionic steps. ASE independently reads the extxyz
  positions/cell/energy/forces and POSCAR coordinates. The full Python adapter
  suite passes with a current native wheel, including existing force/unit tests.

## External PseudoDojo runs and the MgO limitation

The readable external directory was
`/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex`; the environment variable
was set explicitly for the ignored real-file test. UPFs are not committed.
Hashes, fixture source hashes, commands and logs are in the evidence manifest.
Both tests use fixed atom 0, PBE with density floor 1e-10, strict SCF tolerances
(density 1e-9, energy 1e-11 Ha, eigen 1e-10 Ha), Gamma sampling and a 0.03 Bohr
Cartesian displacement on each axis of atom 1.

Si uses primitive diamond a=10.2612 Bohr, cutoff 30 Ha. It returns to fractional
(0.25,0.25,0.25) within 1.8e-10 Bohr in each tested axis, converges in 5 steps
and has final maximum free force 5.570571e-11 Ha/Bohr. Its energy decreases
monotonically within 1e-9 Ha numerical noise. Repeating with a 17-iteration SCF
budget yields two rejected trial SCFs and exit 4, retaining good ionic step 3.
The records include both failed candidate geometries and their causes.

MgO uses a=7.96 Bohr, cutoff 60 Ha. It converges monotonically in 57 steps to
maximum free force 3.967598e-9 Ha/Bohr and oxygen fractional coordinates near
(0.75,0.75,0.75), **not** the input rocksalt (0.5,0.5,0.5). The initial test's
extra ideal-position assertion therefore failed; that failed run is retained.
The ticket explicitly requires ideal-bond recovery for Si, but only a displaced
MgO oracle. MgO acceptance here means a certified, energy-lowering stationary
state of this Gamma-only model; it does not qualify bulk MgO structure.

An independent fresh-SCF scalar energy-curvature check at the ideal MgO
position resolves negative curvature, approximately -0.08178 Ha/Bohr^2 along
[111], across three displacement sizes. This confirms that demanding return to
the ideal position is inappropriate for this finite Gamma model. No force or
SCF tolerance was loosened. K-mesh-converged MgO structural validation remains
open; neither this minimum nor the Gamma-only Si test establishes broad
pseudopotential transferability or global-minimum behavior.

## Reproduction

```sh
scripts/check.sh
cargo test -p dft-cli --lib relax -- --nocapture
cargo test --release -p dft-engine --test density_transfer -- --include-ignored --nocapture
cargo test --release -p dft-cli --test relaxation -- --include-ignored --skip real_si_and_mgo --nocapture
DFT_NC_DIR=/path/to/nc cargo test --release -p dft-cli --test relaxation real_si_and_mgo -- --ignored --nocapture
python3 -I tests/reference/relaxation/check_mgo_curvature.py --binary target/release/dft-cli --nc-dir /path/to/nc
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /absolute/path/Si.pz-vbc.UPF --cases gamma,k222,si8,k444s --repeats 1 --baseline tests/reference/forces/followup-legacy.json --require-identical
python3 -I scripts/tests/test_readme_figures.py
```

The legacy benchmark reports delta E=0 Ha and identical density bytes for all
four cases. Wall times on this Mac are observations, not performance acceptance.
The existing README freshness suite passes unchanged. The long SCF trajectory
and real-file tests are ignored in debug to preserve the fast-gate contract.

Python: build a release native wheel with `maturin build --release --manifest-path
crates/baandi-py/Cargo.toml --features python-extension --interpreter PYTHON`.
Use that extension with `python/baandi` on `PYTHONPATH`, ASE installed, and
`BAANDI_CLI=/absolute/path/target/release/dft-cli python -m unittest discover
-s python/baandi/tests -v`. This run extracted the extension to `/tmp/p21-python`
and used the already available `../pyenv/bin/python`; no shared environment
installation was performed. The initial ASE exception failure is retained;
raising ASE's `PropertyNotImplementedError` enables its documented fallback to
ordinary energy without adding a free-energy alias.

## Scope and orchestrator handoff

- Spin/finite-temperature ionic relaxation, changed-cell/grid transfer, orbital
  transfer and cross-build/cross-platform exact continuation reject or are not
  claimed. Initial density transfer is not exact SCF checkpoint restoration.
- The original `stage-oracles.md` lacked MF15; its new row records the criteria
  already explicitly supplied by this ticket, without changing ADRs 0033/0034.
- The README result generator, generated README/figures, dashboard registry and
  upload/dataset actions are outside this ticket's `OWNED_FILES`. The orchestrator
  must add this result type to `scripts/make_readme_figures.py`, regenerate the
  README block, and update/upload the dashboard artifacts. No new README claim
  was hand-entered, and no task registry status was fabricated.
- `docs/schemas/` is frozen: the orchestrator should add the tagged
  `calculation_type="relaxation"` results profile and optimizer restart schema
  if formal external JSON-schema validation is required. Existing SCF schemas
  and their consumers are unchanged.
- Changes remain uncommitted on `task/P21A-RELAXATION` for review.
