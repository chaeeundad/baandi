# P27C: radial convention residual diagnostic

2026-11-07. Test-only diagnostic; no engine, default, tolerance, dependency or
UPF change. No QE source was accessed. Implementation:
`crates/dft-physics/tests/radial_convention_diag.rs`.

## Conclusion

**No measured candidate establishes an explanation of the positive residual.**
Dropping the final O point raises the fixed-density local energy by
7.6682999881e-7 Ry, so BAANDI minus that alternative is **negative**, equivalent
to -59.96% of the +1.279e-6 Ry MgO target. It explains none of the positive
residual; applying it would increase the predicted residual to +2.04583e-6 Ry.
Its zero-mode-only contribution is +5.2494331965e-7 Ry in the alternative-minus-
BAANDI direction (-41.04% in the comparison direction). Si is unchanged.
The nonzero-G local contribution is significant and must not be omitted.

Trapezoidal local integration accounts for only +0.2914% of the MgO target
in the BAANDI-minus-alternative direction, and -1.168% of the Si target.
The tested local interpolation accounts for +0.2645% MgO and -0.02620% Si.
Core XC shifts are at most 3.86e-12 Ry in these alternatives. Projector
transform sensitivities are measured below, but their nonlocal energy
contribution is **not measured**: the retained density files contain no occupied
orbital coefficients. These are not bounds on the complete interpolation error.

Moving the radial upper limit produces effects much larger than the residual,
with species-dependent signs. This demonstrates convention sensitivity, not a
reference mismatch. No reference endpoint, quadrature or interpolation choice
is inferred from these results, and neither convention is called more correct.
There is no new QE run or newly established QE-side convention. The targets
+7.97e-8 Ry (Si) and +1.279e-6 Ry (MgO), BAANDI minus QE, are the rounded
P27A predictions, not newly measured post-change SCF residuals. See
[xc-prefactor-diagnostic.md](xc-prefactor-diagnostic.md) and its original
[relayed reference provenance](../../tests/reference/perf/dev08_materials_qe_relayed.md).

## Actual engine rules and energy derivation

`dft-upf/src/radial.rs` integrates in mesh-index coordinates: each sampled
integrand is multiplied by PP_RAB. Composite Simpson 1/3 panels use a final
3/8 panel when the point count is even. The independent scalar diagnostic
assembles weights then uses compensated summation, and asserts agreement with
production transforms within 1e-11 absolute in each transform's stated units.
Polynomial tests check odd/even Simpson exactness through degree three and
trapezoid exactness for a linear function. No acceptance cap was relaxed.

`LocalTransform` selects the last point **at or below** min(10 Bohr, mesh end),
not the first point beyond it. Si/Mg have 1510 points, O 936; the default local
prefix has 1001 points for Si/Mg and 936 for O. Odd truncation therefore affects
only O's local transform. The retained O endpoint changes 9.35 -> 9.34 Bohr.
Projectors integrate their own stored prefix, and core density integrates the
full mesh; neither has the local 10 Bohr limit. `data.rs` only reexports types
and erf/erfc; it has no additional integration rule. `hamiltonian.rs` computes
direct local/projector transforms at the needed wave numbers; `density.rs`
computes direct core transforms. No production radial interpolation table is used.

Let S_s(q) = integral r^2 [V_s(r)+Z_s erf(r)/r] j0(qr) dr in Hartree units.
Define A_s = 4 pi S_s(0) + pi Z_s (Ha Bohr^3). Production `alpha_z()` returns
A_s/Omega, already including Z and the volume factor. Nonzero-q local
coefficients before volume division are
L_s(q) = 4 pi [S_s(q)-Z_s exp(-q^2/4)/q^2].
`ionic_local` sums these divided by Omega with ionic phases; at G=0 each phase
is one. `dft-engine/src/energy.rs::local_energy` contracts
E_local = Omega sum_G Re[n_G* V_G], and the audit records
`local_zero_mode = charge * ionic_g[0].re`. Therefore

    delta E_zero (Ha) = N_el / Omega * sum_s N_atoms,s delta A_s.
    delta E_local (Ha) = sum_G Re[n_G* sum_I exp(-i G.R_I) delta L_s(G)].

If alpha_s = A_s/Z_s, the first formula is N_el sum_s N_atoms,s Z_s delta
alpha_s/Omega. A single alpha times Z-sum is valid only for a common alpha;
it must not double-count Z or use one species' alpha for a mixture.
Ry values multiply the Hartree result by two. Volumes are 2*(5.13)^3 =
270.011394 Bohr^3 (Si, 2 atoms, 8 electrons) and 2*(3.975)^3 =
125.61496875 Bohr^3 (MgO, 16 electrons).

Fixed-density local shifts are exact explicit linear energy differences on the
retained converged density, and first-order stationary total-energy predictions
when only the local potential changes. No perturbed SCF or orbital response is
computed. Core transforms are reconstructed with the same density sphere,
phases and volume, then evaluated with production PBE on valence plus core;
these are finite fixed-valence XC differences, including nonlinear terms.
The baseline local energy reproduces stored CLI values to 2.38e-14 Ha (Si)
and 5.05e-13 Ha (MgO); core/XC reconstruction agrees within 1e-12 Ha.

## Experiments

Alternatives are applied one at a time. `odd` drops the final mesh sample when
the selected local/core mesh count is even. Projector prefixes are only
shortened if they extend past the resulting full mesh endpoint; their own
prefix parity is not changed. `r9`, `r9.5`, `r10.5` select samples
at or below that radius, clamped to available data: O cannot be extended beyond
9.35 Bohr, and no tails are extrapolated. For core/projectors this additionally
tests imposing the named cutoff on their full arrays. `trap` integrates the
same PP_RAB-weighted samples with composite trapezoid in index coordinates.
These real files have uniform 0.01 Bohr meshes.

`table` is a diagnostic four-point Lagrange table at spacing 0.01 Bohr^-1,
using floor(q/dq)-1 through +2, clamped at zero. Only the smooth erf-subtracted
local integral is interpolated; the Coulomb pole stays analytic and alpha-Z
stays directly evaluated. Core and projector radial factors are interpolated
directly. This is a stated alternative, not an assertion that it exactly
reproduces a reference implementation's table or endpoint handling.

All deltas below are **alternative minus BAANDI**, opposite to the target's
BAANDI-minus-reference sign. Fractions are -delta / target, not fitted values.

| Alternative | Si zero mode (Ry) | Si full local (Ry) | MgO zero mode (Ry) | MgO full local (Ry) |
| --- | ---: | ---: | ---: | ---: |
| odd | 0 | 0 | +5.24943320e-7 | +7.66829999e-7 |
| r9 | +2.99879570e-6 | +2.82161563e-6 | -9.34743415e-5 | -1.24719269e-4 |
| r9.5 | +1.23341948e-6 | +1.11737232e-6 | -5.56733506e-5 | -7.87606107e-5 |
| r10.5 | -8.21013031e-7 | -7.06559728e-7 | +5.31150223e-5 | +4.73216195e-5 |
| trap | +1.28687817e-9 | +9.30930652e-10 | -2.75066601e-9 | -3.72652000e-9 |
| table | 0 | +2.08825111e-11 | 0 | -3.38339728e-9 |

| Alternative | -delta local / Si target | -delta local / MgO target | Si core XC shift (Ry) | MgO core XC shift (Ry) |
| --- | ---: | ---: | ---: | ---: |
| odd | 0 | -0.599554 | 0 | 0 |
| r9 | -35.402956 | +97.513111 | 0 | 0 |
| r9.5 | -14.019728 | +61.579836 | 0 | 0 |
| r10.5 | +8.865241 | -36.998921 | 0 | 0 |
| trap | -0.011680 | +0.002914 | +1.44062540e-12 | +3.85114163e-12 |
| table | -0.000262 | +0.002645 | -7.10542736e-13 | +2.09610107e-13 |

Baseline A_s (Ha Bohr^3): Si 9.34294065587965, Mg 9.23041038857398,
O 6.19432687126154. The following table gives their changes.

| Alternative | Si delta A (Ha Bohr^3) | Mg delta A (Ha Bohr^3) | O delta A (Ha Bohr^3) |
| --- | ---: | ---: | ---: |
| odd | 0 | 0 | +2.0606480842e-6 |
| r9 | +2.5303406478e-5 | -4.4298789793e-4 | +7.6057382650e-5 |
| r9.5 | +1.0407416056e-5 | -2.1854394364e-4 | 0 |
| r10.5 | -6.9275897783e-6 | +2.0850130815e-4 | 0 |
| trap | +1.0858492772e-8 | -1.8442811722e-9 | -8.9533696013e-9 |

## Radial transform probes

The test prints every projector separately, all signed deltas, and baselines
are checked against production at q=0, 0.137, 0.513, 1.017 Bohr^-1. The compact
table below gives maximum absolute changes over the three **positive** probe
wave numbers (and all projector channels). It is a probe summary, not a bound
across the full plane-wave sphere. Units: local Ha Bohr^3 (volume-free),
projector Bohr^(3/2), core electrons. Mg has no core density. Projector energy
requires both occupied orbital overlaps and D couplings; squaring a radial
factor alone is not a total-energy estimator. Such missing energies are not
silently taken as zero.

| Species | Alternative | max local | max projector | max core |
| --- | --- | ---: | ---: | ---: |
| Si | odd | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 |
| Si | r9 | 1.881218350e-05 | 0.000000000e+00 | 0.000000000e+00 |
| Si | r9.5 | 7.585202184e-06 | 0.000000000e+00 | 0.000000000e+00 |
| Si | r10.5 | 4.871449619e-06 | 0.000000000e+00 | 0.000000000e+00 |
| Si | trap | 1.019361662e-08 | 1.544660644e-06 | 1.015743045e-12 |
| Si | table | 4.897628969e-10 | 1.338307243e-09 | 3.950950678e-11 |
| Mg | odd | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 |
| Mg | r9 | 3.281351837e-04 | 0.000000000e+00 | 0.000000000e+00 |
| Mg | r9.5 | 1.591228738e-04 | 0.000000000e+00 | 0.000000000e+00 |
| Mg | r10.5 | 1.464430825e-04 | 0.000000000e+00 | 0.000000000e+00 |
| Mg | trap | 1.810803951e-09 | 1.364751675e-06 | 0.000000000e+00 |
| Mg | table | 1.312514541e-09 | 6.526601481e-11 | 0.000000000e+00 |
| O | odd | 1.542095106e-06 | 0.000000000e+00 | 0.000000000e+00 |
| O | r9 | 5.756338169e-05 | 0.000000000e+00 | 0.000000000e+00 |
| O | r9.5 | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 |
| O | r10.5 | 0.000000000e+00 | 0.000000000e+00 | 0.000000000e+00 |
| O | trap | 1.013489737e-08 | 1.875433090e-06 | 1.391886606e-12 |
| O | table | 9.940777090e-10 | 3.665601156e-11 | 5.916378498e-13 |

## Reproduction and evidence identity

Run from the worktree root (absolute UPF path is needed because Cargo launches
tests from the crate directory). No UPF file is copied or committed. To obtain
the densities again, build the existing CLI and use the existing smoke runner:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex
cargo build --release -p dft-cli
python3 -I scripts/run_materials_smoke.py --cases si-20,mgo-40 \
  --work-dir /tmp/p27c-materials --output /tmp/p27c-materials/smoke.json
RADIAL_DIAG_DENSITIES=/tmp/p27c-materials \
  cargo test --release -p dft-physics --test radial_convention_diag \
  -- --include-ignored --nocapture
```

The CLI rebuild/SCF commands above are reproduction instructions, **not new
runs in P27C**. This run reused P27A's retained `/tmp/p27-xc-final` densities.
Their metadata, SHA-256, layout, units, lattice, grid, convergence flag and
UPF hashes are checked by the test before contraction. Stored local and XC
energies additionally validate the reconstruction. Si uses a 30^3 grid and
MgO 36^3. The density directory is externally retained runtime evidence;
if unavailable, regenerate it using the commands above.

| File | SHA-256 |
| --- | --- |
| Si.upf | `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282` |
| Mg.upf | `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14` |
| O.upf | `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11` |
| si-20/density.bin | `3c48ca15ad7672e3f3a1265a16b6b3ebb31d21f483feb0d4b25c55742a031de6` |
| mgo-40/density.bin | `f5f7e17893f6d85c03402e3d6658822e6109999bd1d584a132a52a18725d6a4a` |

Measured command:

```sh
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
RADIAL_DIAG_DENSITIES=/tmp/p27-xc-final \
cargo test --release -p dft-physics --test radial_convention_diag \
  -- --include-ignored --nocapture
```

Exit 0: three tests passed, including both external diagnostics (not skipped).
The final run took 1.61 s after compilation. Stdout contains the signed
per-projector probes and complete energy tables; captured as
`/tmp/p27c-radial.log`. The two external tests are ignored in the default
fast gate. An unset UPF directory or missing Si/Mg/O files prints an explicit
SKIP. A supplied invalid density directory fails, rather than claiming a
measurement. The pure radial test can run without any density directory.

## Scope and remaining questions

- Nonlocal energy sensitivity is uncomputed because the available retained
  artifacts contain density, not occupied orbitals. The radial probe measures
  all projector channels; no orbital-based bound or complete total-energy
  attribution is claimed. Endpoint and radius alternatives leave the measured
  projector probes unchanged, while trapezoid and interpolation do not.
- Reference-side endpoint, quadrature, smooth-part tabulation and actual table
  interpolation conventions must be established independently before any
  BAANDI-minus-alternative number can be assigned to BAANDI-minus-QE. No
  reference behavior is guessed or tuned to fit the residual.
- No reoptimized orbitals, changed-convention SCFs or new QE evidence are
  included. First-order shifts do not establish higher-order response.
- The narrow ticket owns only this document and the new test. README figure
  generation, artifact registry/dashboard publication and changes outside
  these paths are left to the orchestrator. No new `tests/reference` result
  schema or engine implementation was introduced. Git changes stay uncommitted.

## Verification

| Executed command | Exit | Result |
| --- | ---: | --- |
| Measured external command above | 0 | 3 passed; radial and fixed-density diagnostics executed |
| `env -u DFT_NC_DIR -u RADIAL_DIAG_DENSITIES cargo test --release -p dft-physics --test radial_convention_diag -- --include-ignored --nocapture` | 0 | Polynomial passed; two explicit external-data SKIPs |
| `DFT_NC_DIR=/tmp/p27c-absent cargo test --release -p dft-physics --test radial_convention_diag radial_conventions -- --ignored --nocapture` | 0 | Explicit missing-UPF SKIP |
| `scripts/check.sh` (final source tree) | 0 | PASS build, test, fmt, clippy; `/tmp/p27c-check-final.log` |
| `git diff --check` | 0 | No tracked whitespace errors; new owned files separately checked for trailing whitespace |

The worktree's checked-in session registry has no P27C entry and its local
`.delegate/P27C-MGO-RESIDUAL-DIAG/report.json` does not exist. The orchestrator
must record this run and publish dashboard artifacts from its authoritative
registry; no session status was invented here.
