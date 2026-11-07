# P27A: fixed-density PW92 prefactor diagnostic

2026-11-07. Diagnostic only; no engine behavior, default, tolerance, dependency,
or pseudopotential changed. No QE source was accessed. The implementation is
`crates/dft-physics/tests/xc_prefactor_diag.rs`.

## Result

**Yes for the dominant contribution; partly for the complete offset.** Changing
only the PW92 prefactor from 0.0310907 to 0.031091 lowers the fixed-density
PBE XC energy by 2.1002508994e-6 Ry in si-20 and 6.8705824638e-6 Ry in mgo-40.
These are 94.61% and 99.00% of the relayed XC-term offsets, and 96.34% and
84.30% of the rounded total offsets. The predicted remaining total offsets
are +7.97491006e-8 Ry and +1.27941754e-6 Ry. Thus the hypothesis explains the
sign and most of the magnitude, but does not establish a complete explanation
of MgO. Neither convention is claimed to be more correct.

The comparison uses the existing [dev08 relayed evidence](../../tests/reference/perf/dev08_materials_qe_relayed.md),
including its provenance and rounding limits, and the architect's
[study note](../research/qe-study/qe-pseudo-treatment-differences.md).
No new QE run was made. Printed XC terms belong to different self-consistent
densities; their difference need not equal the explicit fixed-density shift.
The total-energy comparison is the appropriate first-order stationary test.

## Evaluator and conventions

The scalar uniform-gas expression follows PW92, PRB **45**, 13244 (1992),
Eq. (10) and Table I ([paper](https://doi.org/10.1103/PhysRevB.45.13244));
the existing independently transcribed unpolarized coefficients in
`spin_oracles.rs` agree with that expression. PBE follows PRL **77**, 3865
(1996), Eqs. (7), (8), (10), (14)
([author-hosted paper](https://dft.uci.edu/pubs/PBE96.pdf)).
Only the unpolarized PW92 A changes. Other coefficients are
0.21370, 7.5957, 3.5876, 1.6382, 0.49294. PBE retains
kappa=0.804, mu=0.2195149727645171, beta=0.06672455060314922
(full precision as recorded in ADR 0018), and gamma=(1-ln(2))/pi^2.
In particular, gamma is not replaced by either PW92 A.

The independent scalar evaluator does not call xcx. The engine comparison
uses `gga::evaluate`. The core uses the engine's spherical Fourier transform,
density cutoff, volume normalization and ionic phases. The combined density is
valence plus core; no core-only XC energy is subtracted. The explicit material
floor, max(n,1e-10), is applied before FFT differentiation; energy is zero at or
below the GGA threshold. No points encounter either floor in these runs.
Each even-axis Nyquist Miller component is zeroed before multiplication by the
reciprocal lattice, exactly as in the engine, including for skew cells.
The integral weight is volume/grid-point-count. All internal quantities are
Hartree atomic units, with 1 Ha = 2 Ry.

The test invokes the existing material smoke runner, which checks UPF hashes,
SCF convergence, a fresh certificate and electron count. It reads the CLI's
little-endian C-order f64 valence density, verifies its SHA-256, lattice, shape,
units and layout, and reconstructs core density. Agreement with the stored CLI
XC energy additionally checks reconstruction. No Rust dependency was added.

The acceptance cap was fixed at 1e-12 Ha before measurement. On two small
synthetic skew-cell grids, evaluator-minus-engine was 0 and
-3.5527136788005009e-15 Ha. The real-grid comparisons below also pass; the test
asserts this before evaluating the alternative prefactor for each case.

## Primary measurements

| Quantity | si-20 | mgo-40 |
| --- | ---: | ---: |
| Cutoff (Ha) | 20 | 40 |
| FFT grid | 30 x 30 x 30 | 36 x 36 x 36 |
| SCF iterations | 15 | 18 |
| Total energy (Ha) | -8.462317836906042 | -76.05088085171062 |
| Engine XC on reconstructed density (Ha) | -3.0924726204846902 | -12.249880409165776 |
| Modified-A evaluator minus engine (Ha) | 0 | 0 |
| Reconstructed minus stored CLI XC (Ha) | -1.3322676295501878e-15 | -1.7763568394002505e-15 |
| Original-A minus modified-A XC (Ha) | -1.0501254497086165e-6 | -3.4352912319235429e-6 |
| Original-A minus modified-A XC (Ry) | -2.1002508994172331e-6 | -6.8705824638470858e-6 |
| Relayed BAANDI minus QE XC (Ry) | +2.220e-6 | +6.940e-6 |
| Relayed BAANDI minus QE total, rounded (Ry) | +2.18e-6 | +8.15e-6 |
| Predicted remaining total offset (Ry) | +7.9749100583e-8 | +1.2794175362e-6 |

At the SCF minimum, stationarity at fixed electron number makes the
first-order total-energy change equal to the explicit change of the XC
functional on the unchanged density. **The negative fixed-density shifts in
the table are therefore also the first-order predictions including
self-consistency.** A changed-functional SCF has not been run; density response
and higher-order terms have not been measured. The remaining MgO offset must
not be silently attributed to another candidate without a separate test.

## Density and threshold diagnostics

Counts refer to raw valence-plus-core density. Sigma is the squared Cartesian
spectral gradient, in Bohr^-8 (electron counts treated as dimensionless).

| Quantity | si-20 | mgo-40 |
| --- | ---: | ---: |
| Grid points | 27000 | 46656 |
| Minimum density (e/Bohr^3) | 3.2526580185412552e-3 | 1.0005323458352609e-2 |
| Density below 1e-6 e/Bohr^3 | 0 | 0 |
| Density below 1e-10 e/Bohr^3 | 0 | 0 |
| Negative density | 0 | 0 |
| Sigma below 1e-10 Bohr^-8 | 0 | 8 |
| Removing gradient correction at density <=1e-6 or sigma <=1e-10: shift (Ha) | 0 | 0 |

The eight smallest sigma values are:

- si-20: 1.0969409133340512e-8, 3.6000062404330004e-8,
  8.801950422041407e-8, 1.1910446852381593e-7,
  1.191044685298544e-7, 1.1910446853553954e-7,
  2.1155045497017178e-7, 2.1155045501559866e-7 Bohr^-8.
- mgo-40: 1.3297310859718006e-30, 4.3727488225202565e-30,
  5.66805646670577e-30, 1.2359213054653374e-29,
  1.7209657346214353e-29, 4.6426472961733645e-29,
  5.853661997291347e-29, 4.57597366086612e-28 Bohr^-8.

Density-floor and negative-density explanations are excluded for these grids.
MgO has eight nearly zero gradients, but suppressing their gradient correction
changes the energy by zero at f64 resolution. These threshold effects cannot
explain the measured micro-Ry offsets on these densities. This is not a claim
about surfaces, vacuum regions or different materials.

## Separated follow-up: cutoff dependence

All four additional listed cases were run (the ticket calls them “three” but
lists four). These use their own converged densities at their respective
cutoffs, with the A comparison always on the same density within a case.

| Case | Iterations | Shift (Ha) | Shift (Ry) | Evaluator minus engine (Ha) |
| --- | ---: | ---: | ---: | ---: |
| si-30 | 15 | -1.0501249168015647e-6 | -2.1002498336031294e-6 | -4.4408920985006262e-16 |
| si-40 | 15 | -1.0501240939042589e-6 | -2.1002481878085177e-6 | 0 |
| mgo-50 | 19 | -3.4353075353266149e-6 | -6.8706150706532298e-6 | 0 |
| mgo-60 | 19 | -3.4352892441802396e-6 | -6.8705784883604792e-6 | 0 |

Peak-to-peak shifts across all three cutoffs: Si 2.71161e-12 Ry,
MgO 3.65823e-11 Ry. This confirms cutoff independence at the scale of the
observed offsets. All additional density counts below 1e-6, below 1e-10 and
negative are zero. Sigma below 1e-10 counts are 0 for both Si cases and 8 for
both MgO cases; suppressing the gradient correction gives 0 Ha in each.

## Reproduction and evidence

Run from the worktree root. The only build directory is the existing `target`.
All outputs and temporary absolute-path input TOMLs are outside the repository;
UPFs are read in place, never copied. The environment variable must point to
the extracted directory containing `Si.upf`, `Mg.upf`, `O.upf`.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
XC_DIAG_WORK=/tmp/p27-xc-final \
XC_DIAG_CASES=si-20,mgo-40,si-30,si-40,mgo-50,mgo-60 \
cargo test --offline --release -p dft-physics --test xc_prefactor_diag \
  -- --include-ignored --nocapture
scripts/check.sh dft-physics
scripts/check.sh
```

Omitting `XC_DIAG_CASES` selects si-20,mgo-40. Omitting `XC_DIAG_WORK` selects a
PID-specific temporary directory, printed by the test. The ignored test
returns successfully with an explicit SKIP when `DFT_NC_DIR` is absent or any
required UPF is absent. A configured but incorrect UPF fails the runner's hash
check. A missing release CLI, failed SCF, malformed density, or invalid probe
is a failure, not a skip. The default test only evaluates two tiny grids;
material SCFs remain ignored in debug.

The script invoked by the test is `scripts/run_materials_smoke.py`; it reads
`tests/reference/materials/<case>.toml`, checks the manifest/UPF SHA-256 values,
and changes only the runtime pseudopotential paths. Si uses shifted 4x4x4 and
MgO shifted 3x3x3 k meshes with time reversal, four threads, original strict
convergence settings, and density cutoff four times the wavefunction cutoff.
The UPF identities are those in the relayed reference and source TOMLs:

- Si: `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`
- Mg: `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`
- O: `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`

Primary valence density binary SHA-256:

- si-20: `3c48ca15ad7672e3f3a1265a16b6b3ebb31d21f483feb0d4b25c55742a031de6`
- mgo-40: `f5f7e17893f6d85c03402e3d6658822e6109999bd1d584a132a52a18725d6a4a`

Integrated primary valence counts are 8.000000000000059 and
16.00000000000016 electrons. Fresh-certificate density residuals are
2.8361100426239544e-9 and 1.5246186177956903e-9 in the engine's normalized
L2 metric; occupied eigenvalue changes are 1.3106684348951347e-8 and
1.0625336499536786e-8 Ha, with maximum strict residuals
9.928556539536035e-10 and 9.691087789815237e-10 Ha.

Raw CLI results, densities, metadata, convergence certificates and source
input hashes are retained under `/tmp/p27-xc-final/` (`smoke.json` and per-case
subdirectories); these temporary files are not durable repository artifacts.
Measured summaries are preserved in this document. No performance claim is
made from Mac execution times. No engine switch, new SCF under the changed
functional, radial-transform attribution, or beyond-Si/MgO validation is made.

## Verification record

Mac arm64, rustc 1.99.1 (58cf8f9c3 2026-10-13), base commit
`ea7e55b3e2e9b40d1d46fc2f56e8ec753267db07`, branch
`task/P27A-XC-PREFACTOR-DIAG`. Changes are uncommitted as requested.

| Command/check | Exit | Result |
| --- | ---: | --- |
| `cargo build --offline --release -p dft-cli` | 0 | Release CLI built |
| Release diagnostic, all six cases, `--include-ignored --nocapture` with environment above | 0 | 2 passed; six SCFs converged; 27.34 s test runtime |
| `env -u DFT_NC_DIR cargo test --offline --release -p dft-physics --test xc_prefactor_diag materials_prefactor_diagnostic -- --ignored --nocapture` | 0 | Explicit missing-variable SKIP |
| `DFT_NC_DIR=/tmp/p27-no-upfs cargo test --offline --release -p dft-physics --test xc_prefactor_diag materials_prefactor_diagnostic -- --ignored --nocapture` | 0 | Explicit missing-files SKIP |
| Initial `scripts/check.sh dft-physics` | 1 | Build/test/fmt passed; clippy requested `as_chunks` instead of `chunks_exact` |
| Final `scripts/check.sh dft-physics` | 0 | Build/test/fmt/clippy all PASS after byte-reader idiom fix |
| `scripts/check.sh` | 0 | Workspace build/test/fmt/clippy all PASS |
| `git diff --check` | 0 | No whitespace errors |

The byte-reader lint fix changed no formulas or tolerances; the six-case
release diagnostic was rerun afterward. Only the new test and this document
are changed. README, shared dashboard/registry artifacts and production files
are outside this ticket's owned scope and remain untouched.

## P27B: 2026-11-07 explicit convention and self-consistent results

The optional `[xc] pw92_prefactor = "original"` now selects unpolarized PW92
A=0.031091 inside PBE. Omission and `"modified"` retain A=0.0310907 through the
unchanged xcx path. LDA, collinear PBE and unknown values reject the option.
[ADR 0036](../decisions/0036-pw92-prefactor-convention.md) is Proposed; the
recommendation is to retain modified as default. Neither convention is claimed
more correct. No QE source was accessed; no new QE run was performed.

### Self-consistent energies

Same six inputs, UPF hashes, cutoffs, grids, NLCC, density floors and strict
convergence settings as P27A above. The existing smoke runner supplies external
UPF paths and checks hashes, convergence, certificates and electron counts.
Only the explicit convention changes in each paired calculation.

| Case | Modified total (Ha) | Original total (Ha) | Change (Ry) | P27A prediction requested (Ry) | Original iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| si-20 | -8.462317836906042 | -8.462318887031506 | -2.1002509278389425e-6 | -2.1002e-6 | 15 |
| si-30 | -8.462329006518452 | -8.46233005664336 | -2.100249815839561e-6 | -2.1002e-6 | 15 |
| si-40 | -8.462336822426911 | -8.462337872551002 | -2.100248181591269e-6 | -2.1002e-6 | 15 |
| mgo-40 | -76.05088085171062 | -76.05088428700304 | -6.870584854823392e-6 | -6.8706e-6 | 18 |
| mgo-50 | -76.0509692404303 | -76.05097267573814 | -6.870615663956414e-6 | -6.8706e-6 | 19 |
| mgo-60 | -76.05097717077302 | -76.05098060606232 | -6.870578602047317e-6 | -6.8706e-6 | 19 |

Against the *unrounded*, same-cutoff P27A fixed-density predictions, si-20
changes by another -2.8422e-14 Ry and mgo-40 by -2.3910e-12 Ry. These are
consistent with first-order stationarity within the SCF/roundoff resolution;
we do not claim those tiny remainders as resolved density-response energies.
Using the rounded relayed BAANDI-minus-QE totals from P27A gives residuals
+7.9749072e-8 Ry (si-20) and +1.2794151e-6 Ry (mgo-40). The MgO residual remains
unexplained by this ticket. No radial-mesh change or attribution is made.

### Fixed gates and observations

- 20,000 deterministic random points, rs=0.1..20 and exchange reduced
  gradient s=0..3, including sigma=0: maximum relative errors against xcx
  correlation are 3.7601554e-14 (energy/electron), 1.6216984e-14 (density
  derivative) and 1.3035022e-14 (sigma derivative), each below 1e-12.
- Original analytic derivatives versus independent central energy differences:
  maximum relative errors 6.4058202e-10 and 2.0810671e-9 for density and sigma.
  The fixed cap is 1e-7; steps of 1e-5 times natural density/sigma scales balance
  second-order truncation and f64 subtraction. The analytic continuation near
  sigma=0 permits a centered test without taking square roots of sigma.
- Two skew grids: modified total XC error 0 Ha; maximum potential error
  2.2204460e-16 Ha. Original directional energy derivative errors are
  1.2427036e-10 and 9.9432684e-11 Ha (cap 1e-7 times max(1,|derivative|)).
- All six reconstructed real combined-density grids: modified total XC error
  0 Ha, potential error at most 6.6613381e-16 Ha (both caps 1e-12).
- Synthetic two-ion nonlocal PBE fixture, original convention: three Cartesian
  force differences at h=0.0025 Bohr have maximum error 3.3044003e-8 Ha/Bohr
  (existing cap 1e-6); all six stress components at strain h=0.0005 have maximum
  error 2.8668675e-10 Ha/Bohr^3 (existing cap 1e-7). Three step sizes show
  second-order convergence. This is a synthetic fixed-basis derivative oracle,
  not a real-material force/stress or cutoff-convergence claim.
- Five default LDA benchmark cases (gamma, k222s, k222, k444s, si8), each at
  thread budgets 1,2,4,8,12: energy difference 0 Ha and density byte-identical
  to the pre-plumbing executable, with unchanged iteration counts. Timings
  overlapped builds/other tests and support no performance claim.
- Exact checkpoint restoration rejects both original-to-modified and
  modified-to-original with `Rejection::Options`; same-convention restoration
  succeeds. Explicit WarmStart with a new calculation ID succeeds. Wire tests
  verify original round-trip, omitted modified, legacy default and unknown-value
  rejection. Existing NLCC/floor checkpoint exclusions are unchanged.

### Reproduction

All builds use this worktree's single `target/`. Build the baseline CLI from
base commit `6606cc8c8a14b5f8bf0620fd020bf4895f06c842` before implementation
and retain its executable outside the worktree (here `/tmp/p27b-baseline-cli`).
The baseline executable was retained before option plumbing; the new unused
physics module may have been compiled during that build, but no existing
numerical path or caller had been changed. No git metadata was written.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
cargo test --offline -p dft-physics --test pbe_pw92_oracles -- --nocapture
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
PW92_WORK=/tmp/p27b-materials-final PW92_BASELINE_CLI=/tmp/p27b-baseline-cli \
cargo test --offline --release -p dft-physics --test pbe_pw92_materials \
  -- --ignored --nocapture
cargo test --offline --release -p dft-cli --test pw92_prefactor_config \
  -- --include-ignored --nocapture
```

The real-file test skips explicitly when DFT_NC_DIR or required files are
absent; malformed configured inputs and failed SCFs fail. With
PW92_BASELINE_CLI it additionally compares each default PBE energy and density
against the pre-change executable. Explicit modified versus omission is always
checked. Raw inputs, outputs, density binaries and comparison JSON are retained
in PW92_WORK, outside the repository; they are temporary evidence. This section
preserves the measured summary. UPFs were read in place, never copied.

Benchmark: `python3 -I scripts/bench_scf.py --binary target/release/dft-cli
--upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF
--cases gamma,k222s,k222,k444s,si8 --repeats 1 --threads T
--baseline /tmp/p27b-baseline.json --require-identical`, for T=1,2,4,8,12.
The baseline JSON was produced with the retained baseline executable, T=1,
using `--record /tmp/p27b-baseline.json`.

Initial integration tests exposed missing original-path scratch allocation
(exit 101), corrected before material results. The restart test fixture then
needed its actual PBE/cutoff model identity and a distinct warm-start calculation
ID (exit 101); production restart checks correctly rejected the inconsistent
fixture. The initial workspace check exited 1 for that fixture. No tolerance or
production restart check was relaxed. The final verification record follows.

### Final verification record

All following final commands exited 0 on this Mac. All four check.sh runs
reported PASS for build, test, fmt and clippy; gates were not relaxed.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo test --offline -p dft-physics --test pbe_pw92_oracles -- --nocapture` | 0 | Two point/grid oracle tests passed |
| Six-material release test with PW92_BASELINE_CLI, command above | 0 | All six original SCFs converged; real-grid gates passed; baseline default and explicit-modified energy/density identical; 114.90 s including extra SCFs |
| Release CLI convention suite with `--include-ignored --nocapture` | 0 | Config, exact/warm restart and force/stress tests passed |
| `scripts/check.sh` | 0 | All workspace gates passed, including convention wire tests |
| `scripts/check.sh dft-physics` | 0 | All gates passed |
| `scripts/check.sh dft-engine` | 0 | All gates passed |
| `scripts/check.sh dft-cli` | 0 | All gates passed |
| `bench_scf.py --require-identical`, command above, T=1,2,4,8,12 | 0 each | 25 energy/density comparisons passed |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | Nine tests passed; existing unclosed-file ResourceWarnings only |
| Real-file release test with DFT_NC_DIR unset | 0 | Explicit missing-variable SKIP |
| Real-file release test with DFT_NC_DIR=/tmp/p27b-missing-upfs | 0 | Explicit missing-files SKIP |
| `git diff --check` | 0 | No whitespace errors |

A separate JSON comparison of the six default versus baseline outputs also
confirmed exact equality of all `energies`, `total_energy`, `eigenvalues`,
`occupations`, `kpoints`, `electrons`, `bands` and `band_gap` records.
Execution timing is intentionally excluded from numerical identity.

Executable SHA-256:

- Retained baseline: `95639033c78e5eac20f1f7adcaacb1eb6d0b260eda7f5eb65c641205aeca22dd`.
- Implemented release CLI: `a7fb625ef1dd23a989b739dc86debe6966c9b5b7098197ab2cd4d81fab8624f4`.

Original-convention density SHA-256 in case order:

- si-20: `8a676fde37d52ea96d44635f985742c1542890cb011a4c70736dd31f418e3b4d`.
- si-30: `3075c2f1428bf571aed7af0a011cfc46a762968ecb8387977409bf3e78888223`.
- si-40: `009cc0100343de79f649d6a0cf2d7b9f32a806820c09829f65713ec812d77f92`.
- mgo-40: `4d775867de0b829f987447e01eb39308add84a4b7def646a35cf0b8edbe4631c`.
- mgo-50: `30ad9bca77fd6aec82e1fc08b5b7988c58bf942e0d3d65e352ade3845e0b937d`.
- mgo-60: `ae962911d612041b6a059a3bfcae0b928dec6474d7865a84777231d325a04f91`.

Additional engine/CLI files contain option plumbing only: engine `lib.rs`,
`force_identity.rs`, `forces.rs`, `stress.rs`, `spectral_scf.rs`, and CLI
`checkpoint_io.rs` (including wire tests). The derivative implementations are
unchanged apart from selecting the new XC path. No existing physics numerical
path was edited. README, real-gamma-r5-results.md, pseudopotentials and shared
orchestrator dashboard/session artifacts remain untouched per ticket scope.
Changes are uncommitted on `task/P27B-PW92-ORIGINAL-OPTION` for review.
