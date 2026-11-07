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
