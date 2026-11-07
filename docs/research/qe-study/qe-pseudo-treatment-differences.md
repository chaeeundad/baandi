# QE 7.6 versus BAANDI: where the same pseudopotential and functional can give different energies

Status: research note written by Claude Code from reading QE 7.6 (read-only clone outside the repository), in its own words, per the "QE source study policy" in AGENTS.md. It contains no QE code. Routine and file names are references only. Implementers must not be given QE source; the diagnostics below are specified against BAANDI's own design and public papers.

## 0. Problem and the one finding that matters

Six dev08 bulk cases (PseudoDojo NC-SR v0.4 PBE standard, NLCC for Si and O, 4x4x4 shifted mesh for Si, 3x3x3 for MgO, ecutrho = 4 ecutwfc, labels are in Ha: si-20 = 40 Ry) give BAANDI minus QE of about +2.15e-06 Ry (Si, per cell of 2 atoms) and +8.14e-06 Ry (MgO), nearly cutoff independent. The relayed breakdown (`tests/reference/perf/dev08_materials_qe_relayed.md`, addendum) puts 102 percent (Si) and 85 percent (MgO) of the difference in the printed XC term, and shows that Ewald, kinetic, local and nonlocal agree.

**Primary finding (read from source on both sides, the numerical estimate is mine):** the two codes use two different published constants for the uniform-gas PW92 correlation inside PBE.

- QE 7.6, XClib `qe_funct_corr_lda_lsda.f90`, routine `pw`: the prefactor is the originally published PW92 value a = 0.031091. This routine supplies the uniform-gas correlation both to the LDA part of PBE (driver `qe_drivers_lda_lsda.f90`, correlation case 4) and inside the gradient part (`pbc` in `qe_funct_corr_gga.f90`). The PBE gamma inside `pbc` is the exact (1 - ln 2)/pi^2 = 0.0310906909, kept separately.
- BAANDI, through `xcx` 0.3.1, `gga_c_pbe.rs`: the uniform-gas limit is called with the "modified" PW92 set `A_MOD` (a = 0.0310907, plus the exact second derivative of f at zeta = 0), as libxc does. The xcx source states that this differs from the standard set by about 1e-5. The other PW92 coefficients (a1, b1 to b4) are equal in both.

Relative change of the prefactor: (0.031091 - 0.0310907)/0.031091 = 9.6e-06. QE has the larger a, which makes its uniform-gas correlation energy slightly more negative, so QE total energy is lower and BAANDI minus QE is positive. That is the observed sign.

Magnitude. I evaluated the standard PW92 formula (unpolarized) with both constants (a plain scalar script, no density involved). QE minus modified, in Ha per electron: rs = 0.5: -2.93e-07; 1: -1.78e-07; 1.5: -1.24e-07; 2: -9.3e-08; 3: -5.7e-08; 5: -2.8e-08. Relative to ec this is 1 to 4e-06.

- Si: 8 valence electrons at rs mostly 1.5 to 2.5 (mean density 8/270 = 0.030, rs = 2.0), so about 8 x 1.0e-07 = 8e-07 Ha = 1.6e-06 Ry. The NLCC model core adds a little (high density, up to 3e-07 per electron). Estimate 1.7 to 2.3e-06 Ry versus observed 2.2e-06 Ry in the XC term and 2.15e-06 in the total. (Estimate, not a measurement.)
- MgO: 16 valence electrons at mean density 0.127 (rs = 1.2, about 1.5e-07 Ha each) gives about 2.4e-06 Ha = 4.8e-06 Ry, plus the O NLCC core and the denser O 2p region. Estimate 5 to 7e-06 Ry versus observed 6.9e-06 (XC term) and 8.14e-06 (total). The estimate is rough because I did not integrate the actual density.
- Cutoff independence: the difference is a smooth functional of the density, so it cannot depend on the cutoff beyond the density's own convergence. This is what is observed (spread 3e-08 Si, 1e-08 MgO).
- Scaling: the effect is proportional to the electron number times a slowly varying function of rs. Per electron: Si 2.7e-07 Ry, MgO 5.1e-07 Ry; the ratio 1.9 is the rs dependence (denser MgO) times core charge, as expected (rs 1.2 versus 2.0 gives about 1.6 to 2 in the table above).
- Term-by-term breakdown: at a self-consistent minimum a change of the functional changes the total energy at first order only through the explicit functional (second order through the density relaxation, here about 1e-14 Ry given the relative L2 density difference of 1e-7). Hartree (-2.6e-07 Si, +1.0e-06 MgO) and one-electron terms absorb the density response and do not point to a separate cause. Only the total and the XC term are meaningful, and the XC term accounting for 85 to 102 percent is exactly what this candidate predicts.

Ranking: this is candidate X1 below; everything else is at least one to two orders of magnitude smaller or excluded by evidence.

## 1. Ranked candidate list

Rank by likelihood and predicted size, with the pattern test (constant, about proportional to electron number, cutoff independent).

### X1. PW92 prefactor 0.031091 (QE) versus 0.0310907 (BAANDI/libxc modified)

Already described in section 0. Predicted sign: BAANDI higher, matches. Predicted size: about 2e-06 Ry Si, 5 to 7e-06 Ry MgO, matches. Pattern: yes on all counts. Confidence: high (about 90 percent) that this is the dominant part. Possible residual is below.

BAANDI-side diagnostic (no new SCF needed, very fast): add a test/diagnostic-only function in `dft-physics` (not in the hot path, no change of the default) that evaluates the unpolarized PBE XC energy on a given density grid from the published PBE and PW92 formulas, with a selectable prefactor a (0.031091 or 0.0310907) and the exact gamma (1 - ln 2)/pi^2 in the gradient part. Use the retained converged `density.json/bin` of each case with the core. Steps and expected outcomes:
1. With a = 0.0310907, the diagnostic must reproduce the BAANDI XC energy from xcx to about 1e-12 Ha (this validates the diagnostic against the library).
2. With a = 0.031091 the energy must move by about -1.0e-06 Ha (Si, -2.2e-06 Ry) and about -3.5e-06 Ha (MgO, -7e-06 Ry); compare with the relayed XC-term differences +2.220e-06 and +6.940e-06 Ry. Agreement to better than 10 percent confirms; agreement of the total energy difference after also re-running the SCF with the same switch should reduce 2.15e-06 and 8.14e-06 to the 1e-07 level.
3. Optional: make the same switch available as an explicit, default-off XC parameter set (name it for what it is, "pw92 standard prefactor") so the SCF itself can be rerun; the cost is one SCF per case. Keep the default unchanged until the user decides (see section 3).

The diagnostic must be written from the PW92 and PBE papers (Perdew and Wang, PRB 45, 13244 (1992); Perdew, Burke, Ernzerhof, PRL 77, 3865 (1996)), not from QE source.

### X2. GGA density and gradient thresholds in QE (rho < 1e-6, |grad rho|^2 < 1e-10)

QE: in the GGA driver the whole gradient correction (energy and both potential pieces) is set to zero where the total density (valence plus core) is at most 1e-6 or the squared gradient is at most 1e-10 (parameters in `dft_setting_params.f90`, defaults in `dft_setting_routines.f90`; the LDA part has its own cut 1e-10 on the absolute density in `v_of_rho.f90`/the LDA driver). BAANDI: the full PBE (library) is evaluated down to a floor of 1e-10 e/Bohr^3, with a zero below it (`gga.rs`, `functional.rs`, `density_floor` in the input).

Size: the gradient correction scales like rho x (small s) in the tails. In bulk Si and MgO the minimum total density is orders above 1e-6 (valence minima around 1e-3 to 1e-2 in the bonding region, higher with the NLCC core near atoms), so the number of affected points is zero, and the energy effect is below 1e-12 Ry. Points with |grad|^2 < 1e-10 sit at isolated symmetry points; the skipped correction there is of order rho x sigma-small, below 1e-12. Pattern: predicts nothing at the observed level. Excluded by estimate; confirm with a count.

Diagnostic: report from the BAANDI result the number of grid points with total density below 1e-6 and below 1e-10 and with squared gradient below 1e-10, and the energy change when forcing QE's cuts (threshold switch in `gga::evaluate`, default off). Expected: 0 points, change 0.0 to 1e-12.

### X3. Negative or tiny densities and the NLCC core in the XC energy and potential

QE: the XC energy uses the signed total (valence + core); the LDA part evaluates the energy per electron at the absolute value of the density but multiplies by the signed value; the model core is not clipped (a comment in `set_rhoc.f90` explains why the positivity clipping was removed; the core charge is built in reciprocal space inside the density sphere, then transformed, so it may ring slightly negative). The XC potential's "double counting" integral uses the valence density only. The gradient is taken from the G-space total (valence plus core) inside the density sphere. The core-only XC energy (E_xc of the core alone) is not subtracted in 7.6 (the call is commented out).

BAANDI: `core_density` (`density.rs`) builds the core in reciprocal space inside the density basis with the same normalization; `evaluate_with_core` (`functional.rs`) forms valence + core, evaluates XC, and the potential integral uses the valence density only; the GGA rejects negative densities and zeroes below 1e-10 (`gga.rs`). No core-only energy is subtracted. The two assemblies are equivalent by reading.

Size: zero when the total is positive everywhere (the audit reports negative charge 0.0 for both input and output in `results.json`). Core ringing is a high-G truncation effect and would fall with cutoff, which is not observed. Predicts the wrong pattern. Excluded unless the diagnostic below finds negative points.

Diagnostic: confirm the minimum of the total density and the number of negative points in the retained density plus core grid (expected: none), and the model-core number of electrons from the grid versus from the radial integral (expected equal to 1e-8).

### X4. BAANDI density_floor policy (1e-10 e/Bohr^3)

The floor is applied only where the total density is below 1e-10 (setting it to the floor value for the energy and zeroing the potential there, `evaluate_with_floor`). With no such points it does nothing. Same pattern argument as X2: excluded by estimate; the count in X2's diagnostic confirms it. Do not change the floor.

### X5. Hartree, and the G-space density comparison

The relayed numbers put Hartree at -2.6e-07 (Si) and +1.0e-06 (MgO). These are density-response terms of the changed XC (X1), by the first-order argument above: the SCF densities differ at relative L2 1e-7 because the potentials differ by the XC constant change, and the Hartree energy is quadratic in the density, so a relative density change of 1e-7 on E_H of 1.1 Ry (Si) and 41 Ry (MgO) gives 1e-7 to 4e-6 Ry depending on the sign structure. The magnitudes are consistent with the density response (about 2e-07 and 1e-06). Not an independent cause. Both codes sum the Hartree energy over the density sphere with G = 0 omitted, so no lattice-sum difference is expected.

Diagnostic: evaluate the Hartree energy of BAANDI's own density and of QE's density (as read from the QE density file by the existing comparison script) with BAANDI's Hartree routine. The difference should reproduce the relayed one; if it also shows up with identical density, the Hartree path is the cause (not expected).

### X6. Radial interpolation tables for beta(q), the local potential and the core (QE `dq` = 0.01)

QE (upflib `beta_mod.f90`, `vloc_mod.f90`, `rhoc_mod.f90`): each radial transform is tabulated on a uniform grid in |q| with spacing 0.01 Bohr^-1 and interpolated by the four-point Lagrange formula in the table index; the table is built with the Simpson rule on the mesh. BAANDI computes the transform directly at every |G| (`radial.rs`).

Size: the Lagrange error is about 0.56 x dq^4 x f''''/24 with f'''' of the order r^4 f, r the radial extent of the function (2 to 3 Bohr for projectors), so the relative error is about 1e-8 to 3e-8 and alternates in sign with q, so its first-order effect on the energy averages to below 1e-7 Ry. It does not scale with electron number and would carry a cutoff dependence (more G shells average it out). Predicts at most 1e-8 to 1e-7 Ry, below the observed offsets; likely source of the small cutoff dependent part (si-20 versus si-30 differs by 3e-08 Ry; MgO by 8e-09 Ry).

Diagnostic: in a BAANDI test-only path, round the projector and local transforms through a 0.01 table with four-point Lagrange (written from the formula, which is a textbook interpolation rule) and report the energy change on the fixed density. Expected below 1e-7 Ry.

### X7. Radial integration upper limit, mesh handling, and the G = 0 "alpha Z" term

QE (`Modules/read_pseudo.f90`): all local, core and projector-independent integrals stop at the first mesh point beyond 10 Bohr, then forced to an odd count (Simpson); the Simpson variant for an even count uses a special end formula. The local potential is split with an error-function term (Z erf(r)/r, Ry units with e2 = 2), integrated, then the analytic transform is added back; the G = 0 value is the integral of r(rV + Z e2) over the same limited mesh. BAANDI (`LocalTransform`) uses the same erf split, the same 10 Bohr limit, the standard Simpson rule for an odd count and a three-eighths end panel for an even count, and the analytic pi Z / Omega term for the G = 0 constant.

I computed both choices from the Si, Mg and O PseudoDojo files (scalar script, no engine). Results:
- Si and Mg meshes run to 15.09 Bohr: the limit selects 1001 points in both codes, odd, same rule; alpha agrees to 1e-14 Ry Bohr^3.
- O mesh ends at 9.35 Bohr with 936 points: QE drops to 935 (odd), BAANDI uses 936 (even, three-eighths panel). alpha differs by -4.1e-06 Ry Bohr^3 (BAANDI lower). The energy shift is N_el x delta alpha / Omega = 16 x 4.1e-06/125.6 = -5e-07 Ry for MgO: opposite sign to the observed offset and about 6 percent of its size. This is a real, small, predicted residual to look for after X1 (it would make MgO after X1 appear about 5e-07 Ry lower than QE).
- The local potentials in these files have a non-decaying tail of r V + 2 Z of order 1e-5 (noise in the stored digits) out to the end of the mesh. alpha therefore depends on the integration limit: Mg changes by about 1.5e-03 Ry Bohr^3 per Bohr near 10 Bohr (alpha(10) = 18.4608, alpha(15) = 18.4638), Si by about -5e-05 per Bohr. A limit mismatch of 0.01 Bohr would be 2e-06 Ry for MgO. Both codes use 10 Bohr, so there is no difference today, but this is a strong convention dependence: neither choice is "more correct", and any future change of the BAANDI limit must be treated as a model change.

Pattern: not a candidate for the main offset (sign and species sensitivity do not fit). Diagnostic: add an option to the local transform to use the QE-style odd count (drop the last point of an even mesh) and report alpha and the energy; expected O alpha shift +4.1e-06, MgO energy +5e-07 Ry.

### X8. Ewald parameters

QE (`PW/src/ewald.f90`): alpha starts at 2.9, is reduced in steps of 0.1 until a bound for the truncated G sum is below 1e-7 Ry, the real-space sum stops at 4/sqrt(alpha), the G sum covers the density sphere. For ecutrho above about 200 Ry alpha is 2.8 and both sums are converged to far below 1e-7. BAANDI (`ewald.rs`) uses a truncation budget of 1e-13 relative. The relayed differences are 3e-08 (Si) and 0 (MgO). Excluded by the evidence; the small si-20 versus si-30 change (3e-08) is plausibly QE's alpha stepping at ecutrho = 160 Ry (alpha about 2.2).

### X9. Other items checked and excluded

- Slater exchange constant, PBE exchange kappa and mu, PBE gamma and beta: equal to 1e-12 relative in both (QE uses full-precision constants; its truncated (3/(4 pi))^(1/3) = 0.6203504908994 is accurate to 3e-14).
- Local potential tail term -Z exp(-G^2/4)/G^2, projector normalization, Y_lm, units (UPF in Ry, BAANDI in Ha): equal by the one-electron agreement of 2.1e-07 Ry.
- Valence density sphere: for any k the product of two wave functions has components only inside |G|^2 <= 4 ecutwfc, so QE's sphere truncation of the density changes nothing; likewise the sphere-truncated divergence of the GGA flux only drops components that do not couple to the Hamiltonian.
- Kinetic cutoff test |k+G|^2 <= ecutwfc, k-point weights, fixed occupations: part of the one-electron agreement.

## 2. Which candidates predict the observed pattern

| Candidate | Sign (BAANDI - QE) | Si, Ry | MgO, Ry | Cutoff independent | Scales with electrons | Verdict |
| --- | --- | --- | --- | --- | --- | --- |
| X1 PW92 prefactor | positive | about 2e-06 | 5 to 7e-06 | yes | yes (times rs factor) | matches |
| X2 GGA thresholds | positive or zero | below 1e-12 | below 1e-12 | n/a | n/a | excluded by estimate |
| X3 core and negative density | any | 0 | 0 | no (ringing falls with cutoff) | no | excluded |
| X4 density floor | zero | 0 | 0 | n/a | n/a | excluded |
| X5 Hartree | follows X1 | -2.6e-07 | +1.0e-06 | yes | n/a | response, not a cause |
| X6 radial tables | random | at most 1e-8 to 1e-7 | at most 1e-7 | weakly dependent | no | too small, may explain the 3e-08 drift |
| X7 alpha and mesh end | negative for O | 0 | about -5e-07 | yes | yes | small opposite residual |
| X8 Ewald | none | 3e-08 | 0 | alpha stepping | no | too small |

Only X1 produces the whole pattern. The numerical size of X1 is an estimate within about a factor of 1.5; step 2 of the X1 diagnostic turns it into a measurement.

## 3. Recommended order of experiments

1. X1 diagnostic on the retained densities (minutes of work, no SCF): evaluate the XC energy of the BAANDI density with both prefactors, compare with +2.220e-06 and +6.940e-06 Ry. Decisive.
2. Count test X2 to X4 (points below 1e-6 and 1e-10, negative points, gradient below 1e-5): expected all zero. Cheap, closes three candidates with data.
3. If X1 confirms, rerun the six SCFs with the prefactor switch (default-off) and re-compare against the stored QE outputs: the total difference should fall to the level of the remaining items (X7 about -5e-07 for MgO, X6 up to 1e-7, X8 3e-08). Residuals above 1e-6 after that mean something else is open and the list must be re-ranked.
4. X7 odd-count option for the O mesh, to see the predicted -5e-07 in MgO and remove it from the residual.
5. X6 table-rounding diagnostic, to attribute the remaining 1e-8 to 1e-7.

No tolerance in any existing comparison should be loosened at any step; the explained difference is documented as a named model difference with its own measurement, not absorbed into a cap.

## 4. "More correct" versus "different convention"

For X1 the decision is a convention, not a defect, unless an oracle says otherwise:

- Original PW92 (the paper's Table I) gives a = 0.031091. The "modified" value 0.0310907 is the same coefficient with more digits, chosen to agree with the exact high-density-limit coefficient (1 - ln 2)/pi^2 = 0.03109069 (about 0.0310907). PBE itself uses the exact gamma in its gradient term and takes the uniform-gas correlation from PW92 without prescribing the digits. So neither value is wrong, and the two differ by 1e-5 relative in a (about 1 to 4e-06 relative in the uniform-gas correlation energy).
- BAANDI could only be called "more correct" against an independent oracle that defines the target: for example the exact high-density limit of the correlation energy (a property of the physics, favors 0.0310907 asymptotically at rs -> 0 only, and the fit's other coefficients were optimized with the original digits), or a stated reference implementation. A system-level oracle (all-electron or other code) cannot resolve a 1e-06 Ry difference because their own functional constants differ by the same amount.
- What BAANDI should do about it is the user's decision: keep the modified default (libxc convention, the high-density-limit consistent one) and offer an explicit selectable "standard PW92 prefactor" for QE-equivalence tests, or switch the default. The note recommends the first: keep the default, add the named switch, and report the comparison with it enabled as the QE-equivalent number. Record the choice in an ADR; the evidence for it is the X1 diagnostic and the post-switch SCF comparison.
- For X7 the question is similar: the tail of the local potential in the file is noise, the true infinite-range integral is not defined by the file, and a cut at 10 Bohr is a convention; QE-equivalence is achieved by reproducing the odd-count rule, "more correct" has no oracle.

The numerical claims in section 0 (the table of QE minus modified per electron, the Si and MgO estimates) are my own estimates; the alpha and mesh numbers in X7 are computed from the PseudoDojo files with a scalar script and are not engine output. Neither is recorded as a QE reference number. The three source facts used (QE prefactor 0.031091 in `pw`, xcx modified prefactor 0.0310907, same coefficients otherwise) were read directly in the two source trees.
