# QE 7.6 radial endpoint convention, and the endpoint variants to score against the MgO residual

Status: research note by Claude Code from reading QE 7.6 (read-only clone outside the repository), in its own words, per the QE source study policy in AGENTS.md. No QE code, excerpts or tables are reproduced; routine and file names are references only. Implementers must not be given QE source. UPF files were read read-only from `../rust-dft-wt/pseudo/pseudodojo` (extracted to a scratch directory); nothing from them is copied or committed.

Context: after `pw92_prefactor=original` the BAANDI minus QE total energy is +4.5e-08 to +7.6e-08 Ry (Si) and +1.27e-06 Ry (MgO). P27C (`docs/validation/mgo-residual-diagnostic.md`) showed the local-potential G=0 term is extremely sensitive to the radial upper limit.

## 1. How QE chooses the radial limit and what it does with the count

Fact base: `Modules/read_pseudo.f90` (the limit), `upflib/vloc_mod.f90`, `upflib/rhoc_mod.f90`, `upflib/beta_mod.f90`, `upflib/simpsn.f90`, `upflib/read_upf_new.f90`.

Two separate counts exist per species.

**A. The "msh" count (local potential and core charge).** A fixed limit of 10.0 Bohr is a parameter of the reader. For each species the mesh is scanned from the first point and `msh` is set to the index (1-based) of the FIRST point whose radius is strictly greater than 10.0 (the comparison is "r greater than limit", so a point exactly at 10.0 is not the trigger). If no point exceeds 10.0 (mesh ends at or below 10.0), `msh` is the full mesh size. Afterwards, for every species and in all cases, the count is forced odd by an integer-division rule: msh -> 2*floor((msh+1)/2) - 1. For an even count this subtracts one; an odd count is unchanged. Consequences, with exact decimal meshes r_i = 0.01 (i-1):
- first point beyond 10.0 is index 1002 (r = 10.01) -> even -> msh = 1001 (r = 10.00, included, exactly 10.0 as parsed);
- mesh that ends at or below 10 with an even size (O: 936 points, r_max = 9.35) -> msh = 935 (r = 9.34): the last file point is dropped;
- the comment in `atom.f90` says "last point below the limit", which is not what the code does (the code starts from the first point beyond, then forces odd); for a mesh where the first-beyond index is odd the result would include one point past 10 Bohr. Not the case for these files.

**B. The "kkbeta" count (nonlocal projectors).** The beta transform does NOT use msh and has no 10 Bohr limit and no odd forcing. Each species integrates all of its beta functions over kkbeta points, where kkbeta is the largest of the per-projector cutoff indices stored in the UPF (`cutoff_radius_index`, read in `read_upf_new.f90`). Samples beyond a projector's own index are zero in the file, so the common count only matters for projectors with a shorter index (not the case in these files, all projectors of a species share one index).

What the count is used for: msh bounds (i) the erf-split local-potential transform for every |q| of the table (and the G=0 constant), (ii) the core-charge transform, and (iii) the atomic-charge routines (not relevant to the converged energy). kkbeta bounds the beta transform. All tabulated at spacing 0.01 Bohr^-1 and interpolated (see `qe-pseudo-treatment-differences.md`, X6); the table build is where the count enters.

## 2. Quadrature rule and end weights

The one routine `simpson` (`upflib/simpsn.f90`) is used for all of them, as a sum over mesh index with the supplied `rab` weights (integrand times rab, composite rule in the index coordinate).
- Odd count n: the ordinary composite Simpson 1/3 rule over all n points (end weights 1/3).
- Even count n: composite Simpson on the first n-1 points (odd count), plus a three-point quadratic-interpolation integral over the LAST interval only (through points n-2, n-1, n; equivalent to weights -1/12, 8/12, 5/12 on those three points; as a result the net weights on points n-2, n-1, n are 15/12, 1 and 5/12 in the units of the 1/3 prefactor convention). The source comment says the formula follows the DFTK convention.
- BAANDI (`dft-upf/src/radial.rs`, `integrate`): odd n identical (Kahan summation only). Even n: Simpson 1/3 up to point n-3, then a Simpson 3/8 rule on the last four points (cubic interpolation over the last three intervals).
The two even-n rules differ in the last three intervals. They agree to the level of the interpolation error of the integrand there, so the difference is negligible when the integrand is smooth and about zero at the end (local potential, core), and visible when the integrand ends at a hard cut with a nonzero sample (projectors, see section 4).

## 3. Same logic for all species and all integrals?

- Same for all species: yes. The loop over species applies the same limit, parity rule and Simpson routine; there is no species switch.
- Local potential (both q > 0 and the G=0 constant) and the core charge: same count msh, same rule. The G=0 constant uses the integrand r(rV + Z e^2) without the error-function split, over msh; the q > 0 values use the erf-split integrand and add back the analytic Gaussian transform. I verified numerically from the three files that BAANDI's split form plus the analytic pi Z term equals QE's unsplit form at the same count to 1e-14 Ry Bohr^3, so only the count and end rule matter, not the algebra.
- Nonlocal: separate count (kkbeta), no limit, no odd forcing, same Simpson routine.
- BAANDI today: local uses the last point at or below min(10 Bohr, mesh end) (`LocalTransform`, via `partition_point(r <= limit)`), no parity forcing; `ionic_local` in `hamiltonian.rs` re-derives the same n; projectors use their own stored prefix of `cutoff_radius_index` samples; the core transform uses the FULL mesh (`spherical_density_transform`, 1510 and 936 points, even, 3/8 end).

## 4. UPF mesh data and the distinct endpoint variants

Data (read from the files; all three meshes are uniform: r_i = 0.01 (i-1) to 2e-15, rab = 0.01, first point r = 0; r is stored with four decimals so r(1001) parses to exactly 10.0 in both codes):
- Si: 1510 points (r_max 15.09); projector cutoff index 196 (r = 1.95) for all 6 projectors; core density is zero beyond index 462 (r = 4.61).
- Mg: 1510 points (15.09); cutoff index 160 (r = 1.59) for all 4 projectors; no NLCC.
- O: 936 points (r_max = 9.35, below 10); cutoff index 152 (r = 1.51) for all 5 projectors; core zero beyond index 176 (r = 1.75).
Indices below are 1-based mesh points (0-based = minus one); "n" is the point count used.

Local potential and G=0 term (variants L):

| Variant | Si n / last r | Mg n / last r | O n / last r | Rule | Reading it represents |
| --- | --- | --- | --- | --- | --- |
| L0 BAANDI today | 1001 / 10.00 | 1001 / 10.00 | 936 / 9.35 | odd: Simpson; O even: 3/8 end | last point at or below the limit, no parity forcing |
| L1 QE as I read it | 1001 / 10.00 | 1001 / 10.00 | 935 / 9.34 | all odd: Simpson | first point beyond, then forced odd |
| L2 strict "below" | 1000 / 9.99 | 1000 / 9.99 | 936 / 9.35 | even: either end rule | r < 10, no parity forcing |
| L3 strict below then odd | 999 / 9.98 | 999 / 9.98 | 935 / 9.34 | Simpson | r < 10, then forced odd |
| L4 first beyond, no forcing | 1002 / 10.01 | 1002 / 10.01 | 936 / 9.35 | even rule | first point beyond the limit kept |
| L5 first beyond, forced odd upward | 1003 / 10.02 | 1003 / 10.02 | 936 / 9.35 (no point to add) | Simpson | odd forcing rounding up |

Even-count end rule in L0, L2, L4 (BAANDI 3/8 versus the QE-style quadratic end interval) changes the local alpha only at the 1e-9 Ry Bohr^3 level; it is not a distinct scoring variant for the local term (I treat it as one variant per n, with either rule).

Scalar zero-mode sensitivity (my script, plain Simpson on the file data, not engine output; E = N_el/Omega x delta alpha in Ry for the alpha-Z term only; positive means the alternative energy is higher than L0). Nonzero-G local terms add to this: P27C measured the full local change for the O drop as 1.46 times the zero-mode part (+7.668e-7 against +5.249e-7), the factor for Mg is not known and must come from the harness.
- Mg, per mesh point at the 10 Bohr limit: one point fewer -1.091e-06, one more +1.090e-06, two fewer -2.183e-06, two more +2.179e-06.
- O, near the end: drop last point (935) +5.250e-07 in MgO; 934 points (r = 9.33) +1.05e-06. Alpha for O cannot be extended past 9.35 because the file ends.
- Si per point: about 1e-08 (negligible, consistent with Si residual being small).
Sign requirement: the target is BAANDI minus QE = +1.279e-06 Ry, so a variant that matches QE must have alt minus L0 of about -1.28e-06 Ry (QE lower). Single-point shifts of Mg are the only individual changes of the right size and sign (one point fewer: -1.09e-06 zero-mode, probably about -1.6e-06 in full local if the O factor 1.46 applied); the O drop (L1 for O) goes the other way (+0.5 to +0.8e-06). I do not claim which variant is right: my reading of QE predicts L1 (Mg, Si unchanged; O dropped), which on its own would widen the residual to about +2.0e-06 Ry, so L1 does not explain the residual and some other term must carry about -2e-06. The harness should score L0 to L5 in all species combinations (each species may take a different row), not just the rows as written.

Nonlocal projector transform (variants P):

| Variant | Si n / last r | Mg n / last r | O n / last r | Rule |
| --- | --- | --- | --- | --- |
| P0 BAANDI today | 196 / 1.95 | 160 / 1.59 | 152 / 1.51 | even: Simpson to n-3, then 3/8 |
| P1 QE as I read it | 196 / 1.95 | 160 / 1.59 | 152 / 1.51 | even: Simpson on n-1 points + quadratic last interval |
| P2 drop the last sample (odd forcing) | 195 / 1.94 | 159 / 1.58 | 151 / 1.50 | Simpson |
| P3 include the first zero sample (odd, upward) | 197 / 1.96 | 161 / 1.60 | 153 / 1.52 | Simpson |

These differ in the transform because the last stored sample is nonzero for several projectors (a hard cut): Si l=2 (last values 8.5e-05 and -2.7e-04), Mg l=1 (2.1e-05 and -7.3e-05), O l=1 (3.6e-08 and -6.9e-08); the others end at exactly zero and are unaffected. My scalar comparison of P1 minus P0 over q = 0 to 8 Bohr^-1 gives maximum differences (units of the transform, Bohr^-1/2-type radial integral) of 1.2e-07 and 3.6e-07 on Si l=2 (relative 7e-07 and 1.5e-06), 5.7e-08 and 2.1e-07 on Mg l=1 (relative 4e-07 and 1e-06), 4e-10 and 1.8e-09 on O l=1, and below 1e-15 for all other projectors. A relative change of 1e-06 on a semicore projector whose nonlocal energy is of order 1 Ry gives a first-order energy change of order 1e-07 to 1e-06 Ry (estimate; sign not determined here; P27C did not measure the nonlocal energy because it had no orbital coefficients). So the projector end rule is a second candidate of the right order of magnitude for MgO (Mg p semicore) and a smaller one for Si, consistent with a tiny Si residual and a 1e-06 MgO residual. This variant is cheap to score because it needs only the fixed density plus the occupied orbitals, or one extra SCF.

Core-charge transform (variants C):

| Variant | Si n | Mg | O n | Rule |
| --- | --- | --- | --- | --- |
| C0 BAANDI today (full mesh) | 1510 (even, 3/8) | none | 936 (even, 3/8) | |
| C1 QE as I read it (msh) | 1001 / 10.00 | none | 935 / 9.34 | Simpson |

C0 and C1 differ by 4e-16 (Si) and 3e-16 (O) in the transform at q = 0 and 3 because the core density is zero beyond r = 4.61 and 1.75; not a scoring variant. The P27C core XC shifts of at most 4e-12 Ry agree.

## 5. Scoring plan and caveats

1. Score the local variants per species independently (the O row choices L0, L1, L2 and the Mg row choices L0 to L5; Si changes are negligible), with the full local energy, then add the projector variants P0 to P3 (Si, Mg, O); the residual target is +1.279e-06 (MgO) and +7.9e-08 (Si) for the prefactor-original runs. The matrix is small (about 6 x 3 x 4 combinations if species are combined independently).
2. A combination that reproduces both residuals simultaneously, with one coherent rule (a single comparison, a single odd/even rule, the same Simpson), is evidence for that convention; a combination found by fitting different rules per species or per integral is not.
3. Interaction with interpolation: QE tabulates at 0.01 Bohr^-1 and interpolates, which P27C found to be at the 3e-09 Ry level for the local part; the projector and core tables were not scored there and remain a small, separate item.
4. "More correct" cannot be decided from QE: the local-potential file has a noisy tail of order 1e-05 in r V + 2Z (alpha changes about 1.5e-03 Ry Bohr^3 per Bohr for Mg), so the finite-radius integral has no true infinite-range value. Whichever convention reproduces QE is "QE-equivalent". Any change of the default BAANDI radius or rule is a model change needing the user's decision and an ADR; the residual can be reported as a named, measured convention difference without loosening any tolerance.
5. Limits of this note: my scalar script uses plain double Simpson on the stored samples, not the engine's compensated sum; the zero-mode energies convert alpha with N_el/Omega only; QE's actual endpoint was read from source, not measured by a QE run (an instrumented or verbose QE run printing the mesh counts would settle it, if the user wants that).
