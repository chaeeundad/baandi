# QE 7.6 density mixing, specified for an independent implementation (issue #25)

Policy: AGENTS.md, "QE source study policy". Written by Claude Code from a read of QE 7.6 (read-only clone outside the repo), in its own words. Prose and mathematics only. Routine and file names are references. Every statement not confirmed in the source text is marked INFERRED. Scope: Si128-type run (Gamma-only real wavefunctions, nspin = 1, norm-conserving, LDA, mixing_mode plain, defaults otherwise). PAW, DFT+U, meta-GGA, GC-SCF, SIC, electric field terms are skipped on purpose (they add terms to the same bilinear form and the same update, all absent here).

Published method to cite (not QE): D. D. Johnson, PRB 38, 12807 (1988), modified Broyden. Related: Anderson mixing, Pulay/DIIS (Kresse and Furthmueller). Section 3 states exactly in which sense QE equals them.

Primary sources: PW/src/mix_rho.f90 (mix_rho, rho_ddot, high_frequency_mixing, open_mix_file), PW/src/electrons.f90 (electrons_scf), PW/src/sum_band.f90, PW/src/potinit.f90, PW/src/setup.f90 (first ethr), PW/src/input.f90 and Modules/read_namelists.f90 (defaults), Modules/io_base.f90 and PW/src/io_rho_xml.f90 (file), Modules/recvec_subs.f90 (G ordering), PW/src/v_of_rho.f90 (Hartree, charge).

## 0. Notation and units

- Omega: cell volume in bohr^3. rho in electrons per bohr^3 (no unit conversion between Ry and Hartree: a density is the same in both). G in 1/bohr; QE stores |G|^2 in units of tpiba2 = (2 pi / alat)^2 (array gg), so |G|^2 [bohr^-2] = gg * tpiba2.
- QE energies are Rydberg: e2 = 2. Cutoffs in Ry satisfy ecut = |G|^2 [bohr^-2] for a plane wave (hbar^2/2m = 1 Ry bohr^2). ecutrho = 120 Ry means |G|^2 <= 120 bohr^-2 (60 Ha in Hartree convention of a Hartree code: |G|^2/2 <= 60).
- Iteration index k = 1, 2, ... is the SCF counter `iter` of electrons_scf (it is also the value passed to mix_rho and stored as mixrho_iter). x^k = rho_in^k is the density that generated the potential of iteration k, y^k = rho_out^k the density from sum_band at iteration k, F^k = y^k - x^k the residual (sign: output minus input).
- Mixed vectors hold the coefficients rho(G) only (complex numbers), nothing in real space.

## 1. The mixed quantity and its representation

1.1 What is mixed. For nspin = 1 and no PAW/+U/meta-GGA the mixed object (type mix_type) contains exactly one complex array of length ngms: the valence charge density total rho(G). Created by create_mix_type, filled from the scf_type density by assign_scf_to_mix_type for entries 1 to ngms.
- Valence only. The NLCC core charge (rho_core, rhog_core) is never part of rho%of_g; it is added only inside the xc evaluation (v_xc, via rho + rho_core in real space). So the mixed vector, the residual and dr2 do not contain the core charge.
- The ionic/pseudo local potential and Ewald are not density-dependent: not mixed.
- The electric dipole field has a slot in mix_type but is identically zero without tefield/dipfield (its term in the dot product vanishes).

1.2 Representation, normalization. rho%of_g(ig, 1) is the coefficient in rho(r) = sum_G rho(G) exp(i G.r), and rho(G) = (1/N_r) sum_r rho(r) exp(-i G.r) on the dense FFT grid (N_r = 72*72*144 here). Confirmed consequence: the total charge is Omega * Re rho(G=0) (potinit and v_h both compute the charge as omega times the real part of the first coefficient), so rho(G=0) = N_el / Omega. rho(G) is NOT divided by Omega beyond that and carries no sqrt(Omega) factor. The real-space density is the inverse transform of exactly these coefficients (rho_g2r), so rho(r) is real and equals the sum above with no extra factor.

1.3 Gamma-only half sphere. With gamma_only = true only half of the G sphere is stored; the partner coefficient is conj(rho(G)) at -G and is implicit. Kept vectors (Miller indices (i,j,k) with G = i b1 + j b2 + k b3, b = reciprocal basis incl. the 2 pi factor): i > 0, or i = 0 and j > 0, or i = j = 0 and k >= 0 (generation loop of ggen). G = 0 is therefore stored once, as the first entry.

1.4 Index ordering. "G index 1" (ig = 1) is G = 0 (the generator forces it to the front; checked by an assertion that gg(1) is below eps8). The remaining vectors are sorted by increasing |G|^2 with an equality tolerance of 1e-8 in |G|^2 units of tpiba2 (hpsort_eps). Ties inside a shell are resolved by the sorter using the original generation order (loop order i outer, j middle, k inner within the half-sphere filter); that precise tie rule is INFERRED from the sorter's name and use and must be checked against the Miller indices in the density file (section 7). The membership test of the sphere is |G|^2 <= gcutm (inclusive comparison in the generator); for ecutrho = 4*ecutwfc the shell at the boundary may contain vectors whose |G|^2 differs from the cutoff by rounding only: INFERRED risk for an independent generator, so take the vector list from the file, not from a re-generation.
- Under plane-wave MPI the vectors are distributed over ranks by FFT sticks; the local ordering differs per rank, the global ordering (ig_l2g) is the one above. Every global quantity below is rank-count independent up to summation order.

1.5 Dimension of the mixed vector and the high-G cutoff. The mixed vector has ngms entries (also named ngm0 in mix_rho). ngms is the number of G in the smooth FFT grid sphere (gcutms). The smooth grid is distinct from the dense grid only if doublegrid, i.e. dual = ecutrho/ecutwfc > 4 (setup: dual > 4 + eps8). Here ecutrho = 120 = 4 * 30, so doublegrid is false and ngms equals ngm: every stored G is in the mixed vector and the "high-frequency" part is empty. (The equality ngms = ngm in the non-double-grid case is INFERRED from the descriptor reuse; the output of a real run shows it, section 7.)
- General rule for the case ngms < ngm: the entries ngms+1..ngm of the new input density are NOT part of the history scheme. They are set by simple mixing from the OLD input and the NEW output: rho_in_new(G) = x(G) + beta * (y(G) - x(G)) for G beyond the smooth set (high_frequency_mixing). The low part is zeroed and then overwritten by the Broyden result. When ngms = ngm that routine reduces to zeroing everything and then the Broyden result is copied back, so the new input is exactly the mixed vector.
- G = 0 is inside the mixed vector: it is carried through the linear combinations like any other entry. Its residual is zero up to roundoff (y(0) = x(0) = N_el/Omega), the dot product ignores it (below), so the new input has G = 0 equal to N_el/Omega up to roundoff. Nothing renormalizes the mixed density afterwards (INFERRED: no call to a renormalization in mix_rho; the starting density is renormalized once in potinit to N_el if it deviates by more than 1e-7 relative).

1.6 Symmetry. sum_band symmetrizes rho_out(G) with sym_rho before mixing. If the 128-atom cell has symmetry operations beyond identity found by QE (displaced geometry usually has none), y is symmetrized; this is outside mix_rho but changes y. Check the "Sym. Ops." line of the QE output, or run nosym (section 7).

## 2. Residual, startup update, metric dr2

2.1 Residual. At entry mix_rho forms r = y - x on the first ngms coefficients (rhout_m minus rhoin_m). Sign: output minus input. Both vectors are the complete G-space coefficients (full precision complex double).

2.2 The bilinear form (rho_ddot). For two coefficient vectors a, b (same length, half sphere):

<a, b> = (Omega / 2) * e2 * 4 pi * S,   S = sum over G != 0 of Re( conj(a(G)) b(G) ) / |G|^2   (full sphere),

and with the half sphere the sum over stored vectors G != 0 is counted twice: S = 2 * sum_{half sphere, G != 0}. G = 0 is excluded in the default path (the loop starts at the first non-zero vector; the G = 0 term appears only in the GC-SCF branch). In terms of stored data with e2 = 2:

<a, b> = 4 pi Omega * 2 * sum_{ig = 2..ngms} Re( conj(a_ig) b_ig ) / |G_ig|^2   [Ry; a, b in e/bohr^3].

So the factor is 8 pi Omega over |G|^2 (|G|^2 in bohr^-2, which is gg * tpiba2). It is the Hartree energy of the density a (as a charge distribution) with b, in Rydberg: E_H[rho] = (Omega/2) e2 4 pi sum_{G != 0} |rho(G)|^2/|G|^2 (full sphere). In Hartree units E_H(Ha) = 2 pi Omega sum_{G!=0, full sphere}|rho(G)|^2/|G|^2, so <a,a>_QE [Ry] = 2 * <a,a> [Ha-formula with e2 = 1]. A Hartree-based implementation compares dr2 as dr2_QE = 2 * E_H[F] with E_H in Hartree.
- Extra terms present in the function but zero here: spin (nspin >= 2), meta-GGA, +U, dipole. PAW is explicitly commented out.
- The sum runs over ngms vectors (ngm0) for dr2 and for all history dot products; same range everywhere (the line comment says it used to run over all of ngm). With ngms = ngm that is the whole sphere.
- Parallel: partial sums per rank over the local vectors; the G = 0 rank is recognized by the flag gstart = 2 (only affects the GC-SCF branch here); result summed with mp_sum over the plane-wave communicator.

2.3 dr2 and stopping. dr2 = <r, r>, computed before any history is touched, from the residual of THIS iteration (input and output of the same k). Negative dr2 is a fatal error. Convergence: conv = (dr2 < tr2) with tr2 = conv_thr, in Ry, default 1e-6 (read_namelists). If conv, mix_rho returns at once leaving rhoin (input) and rho (output) unchanged (no history update, history buffers released); electrons_scf then takes the "converged" branch (section 4) and the density that is written is the OUTPUT density.

2.4 Early exit by the diagonalization-error rule: also returns unchanged if dr2 < tr2_min. tr2_min is zero except in the very first iteration of the run, where tr2_min = ethr * max(1, N_el) (ethr in Ry, N_el = number of electrons; for Si128 N_el = 512). Section 4.2 gives the follow-up ("redo with lowered threshold").

2.5 Next inner threshold ethr (electrons_scf, at the top of every iteration k > 1): at k = 2 ethr is set to 1e-2 (Ry); at every k > 1 ethr = min(ethr, 0.1 * dr2 / max(1, N_el)) where dr2 is the value returned by mix_rho of the previous iteration; then ethr = max(ethr, 1e-13). Orbital-resolved +U uses 0.05; irrelevant. The threshold is thus monotonically nonincreasing after k = 2. In iteration 1 ethr comes from setup: 1e-2 for a start from atomic charges; 1e-5 if starting_pot = 'file' (unless the input diago_thr_init is non-zero, which is then used unchanged). Adaptive threshold (adaptive_thr) is off by default (default false) and is not used here. The existing note davidson-and-scf-thresholds.md covers the threshold rule from the diagonalization side; the only new fact here is that dr2 is the Hartree-metric number above.
- Unit remark for a Hartree code: ethr is a Ry eigenvalue threshold, ethr_Ha = ethr_Ry / 2, and dr2_Ry = 2 * dr2_Ha-equivalent; the product 0.1 * dr2 / N_el is dimensionally consistent in either unit system if both are converted.

## 3. History-bearing update (modified Broyden, as implemented)

3.1 Parameters. mixing_beta (beta) default 0.7 (internal sentinel -1 replaced by 0.7 for plain runs); mixing_ndim (n_iter, here called m) default 8; hard maximum 25. mixing_mode plain means no preconditioning (imix = 0). `TF` and `local-TF` are the alternatives (section 4.5).

3.2 Counters. With k = mixrho_iter (the SCF iteration):
- iter_used = min(k - 1, m). At k = 1 it is 0.
- ipos (slot receiving the newest difference) = k - 1 for k <= m + 1, then wraps 1, 2, ..., m (k = m + 2 gives slot 1; the oldest difference is overwritten first). At k = 1 no slot is written (ipos evaluates to 0, not used).
- So at iteration k the history is the last iter_used differences, each formed at a distinct earlier iteration; slot order is irrelevant for the result because the result is invariant under permutation of the stored differences.

3.3 What is stored. mix_rho keeps a file/buffer (iunmix, opened by open_mix_file) with: record 1 = residual F of the previous mix_rho call, record 2 = its input x; and for slot i the pair (dF_i, dx_i) in records 2i+1, 2i+2. Each record is the vector of ngms complex coefficients (plus zero-length extras here). At iteration k > 1:

dF_new = F^{k-1} - F^k,   dx_new = x^{k-1} - x^k   (previous minus current, so the sign is minus the forward difference)

stored in slot ipos; then the pair (F^k, x^k) overwrites records 1 and 2. Older slots are re-read unchanged (they hold differences formed in earlier iterations, never recomputed). The pair (F^k, x^k) is the state BEFORE the update of this iteration (x^k = the input density of iteration k, F^k = y^k - x^k). No normalization of the difference vectors in the default build (an optional compile-time normalization macro exists; INFERRED not enabled by default because it is a preprocessor switch that is not set by the standard build).
- Both vectors have the G = 0 entry (tiny) and the same ngms entries. No explicit "startup" memory beyond this: iteration 1 stores only (F^1, x^1).
- When mix_rho returns early (conv or dr2 < tr2_min) nothing is stored; the redo of iteration 1 then calls it again with the same k = 1 and finds an empty history (the stored pair is simply written on the second call).

3.4 Algebra. Let j = 1..p (p = iter_used) index the stored slots, dF_j, dx_j as above, and r = F^k.

A_ij = <dF_i, dF_j>,   c_i = <dF_i, r>   (same bilinear form of section 2.2, same range ngms).

A is built by the upper triangle and mirrored. The inverse is computed with a symmetric indefinite (Bunch-Kaufman) factorization and its inverse (LAPACK DSYTRF then DSYTRI, uplo = U), then the lower triangle is filled from the upper; no regularization, no pseudo-inverse, no check on conditioning beyond the library's info flag (nonzero info is a fatal error). The coefficients are

gamma_i = sum_j (A^{-1})_{ji} c_j

(the matrix is symmetric so the index order is immaterial). gamma minimizes ||r - sum_i gamma_i dF_i||^2 in the metric <.,.> (normal equations), so this is an Anderson/DIIS-type (Pulay-type) least-squares step on the difference vectors, equal in structure to Johnson's modified Broyden with equal weights and a zero-weight first-Jacobian term (the second half of that sentence is INFERRED from the formulas; a term-by-term comparison to Johnson's eq. is left to the implementer).

Then, in place:

x~ = x^k - sum_i gamma_i dx_i,
F~ = F^k - sum_i gamma_i dF_i,

and finally the new input (no preconditioning for plain)

x^{k+1} = x~ + beta * F~ = x^k + beta F^k - sum_i gamma_i (dx_i + beta dF_i).

Equivalent statement: because dx_i and dF_i have the "previous minus current" sign, the combination moves x^k toward the previous points weighted by gamma; with p = 0 the formula is x^{k+1} = x^k + beta F^k.

3.5 beta enters in two places only: the last step above, and the high-G simple mixing when ngms < ngm. It does not enter A, c or gamma.

3.6 Startup. At k = 1 (and at a redo of k = 1): p = 0; x^2 = x^1 + beta (y^1 - x^1). With beta = 1 this is y^1 up to roundoff of the expression x + 1.0 * (y - x). No other damping exists. The `mixing_fraction` of BAANDI has no QE counterpart.

3.7 After update: the new input is placed in rhoin (scf_type): of_g from the mixed vector (entries 1..ngms), of_r by inverse FFT of those coefficients on the dense grid (no clipping of negative values, INFERRED: no clipping statement in the routine).

## 4. Order inside one SCF iteration (electrons_scf)

4.1 Sequence (iteration k):
1. ethr set (section 2.5). Compute deband_hwf (with the old potential, for the Harris-Foulkes energy).
2. rhoin := copy of rho. At this moment rho holds the current input density x^k (the mixed density of the previous iteration, or the starting density), so rhoin = x^k in both G and r space.
3. c_bands: Davidson diagonalization with V built from x^k. In the very first iteration tr2_min = ethr*max(1, N_el).
4. sum_band: new valence density y^k in real space on the smooth grid (gamma-only real-pair path), forward FFT to G space, symmetrization, back to real space. It overwrites `rho` (rho now = y^k). eband and the HWF terms are computed.
5. deband = delta_e (uses y^k).
6. mix_rho(rho = output, rhoin = input, beta, dr2, tr2_min, k, m, ...) on the root pool: computes dr2, tests conv / tr2_min, and, if neither, mixes and writes the new input into rhoin; broadcast of rhoin, dr2 and conv_elec to all pools.
7. First iteration only: if dr2 < tr2_min (the Davidson threshold was too loose compared to the self-consistency error), print the "Threshold (ethr) on eigenvalues was too large" message, set ethr = max(0.1 * dr2 / max(1, N_el), 1e-13), and redo steps 3 to 6 (cycle of scf_step; state: rhoin still x^1 because mix_rho returned unchanged). The redo is done at most once.
8. If not converged: potential from the NEW input: v_of_rho(rhoin): Hartree (G = 0 excluded) and LDA xc on rho_valence + rho_core (real space). descf = delta_escf = - sum_r (rhoin(r) - rho_out(r)) V_new(r) * Omega / N_r, a first-order correction to make the energy variational; then rho := rhoin (so rho = x^{k+1} from here on); set_vrs (V_loc + V_H + V_xc), newd (nothing for norm conserving).
   If converged: V is recomputed from the OUTPUT density y^k for the final energies and forces, descf = 0, and rho stays = y^k.
9. Energy: E = eband + (etxc - etxcc) + ewld + ehart + deband + demet + descf; the "estimated scf accuracy" printed per iteration is dr2 from step 6.

4.2 Therefore the double-counting correction uses (input minus output) with the potential of the new (mixed) density; there is no separate "mixing fraction".

4.3 Stopping order: conv is tested BEFORE mixing and with input/output of the same iteration; there is no extra post-mixing check.

4.4 Loop end: after convergence or after electron_maxstep iterations the density in `rho` is written (write_scf). Consequence of 4.1 step 8: after a converged iteration, rho = y^k (output); after an unconverged last iteration (electron_maxstep reached), rho = x^{k+1} (mixed new input, i.e. the density that would have been used next). Exception: if scf_must_converge = .false. the last iteration is declared converged, so the output is written (INFERRED: follows from line order, not from a run).

4.5 mixing_mode alternatives. `TF` (imix = 1): after the history subtraction F~ is multiplied by the Thomas-Fermi factor |G|^2 / (|G|^2 + a), with a = (12/pi)^(2/3) / r_s and r_s = (3 Omega / (4 pi N_el))^(1/3) [from the formula in approx_screening; |G|^2 here in tpiba2 units, dimensionally a also divided by tpiba2]; only F~, not dx. `local-TF` (imix = 2): a local-density-dependent screening, longer routine approx_screening2 (not read; INFERRED irrelevant for this task). `plain` does nothing between 3.4 and the final step. mixing_mode `potential` is gone.

## 5. Parallel facts

- Mixing is executed by the root pool only (my_pool_id == root_pool); the result is broadcast to other pools (bcast_scf_type for rhoin, mp_bcast for dr2 and conv).
- Inside a pool with plane-wave (G-vector) MPI the mixing vectors are distributed: each rank owns the entries of its G vectors; the vector ops (axpy, copy, scale) are rank-local. Reductions: each bilinear form <a, b> is a local sum plus one mp_sum of a scalar over the band-group communicator; per iteration there are 1 (dr2) + p(p+1)/2 (matrix) + p (rhs) such reductions, a few scalars each, with p <= 8 default. Replicated: the small p x p factorization/inverse and gamma on every rank (bit-identical inputs after the reductions, up to reduction order); loop counters, beta, ethr.
- The history lives per rank (in the mix buffer/file iunmix, record length 2*ngms*nspin reals plus extras); whether it stays in memory or goes to disk follows io_level (INFERRED, not traced).
- Parallel FFT: one inverse transform (rho_g2r) per iteration after mixing for the new real-space density (dense grid). G-space to real-space and back for sum_band belong to the density build, not to the mixing.
- The rank holding G = 0 is flagged by gstart = 2; for the plain-mode dot product the G = 0 entry is simply never summed, so no special reduction is needed.
- Band-group and pool parallelism only affect who computes the density; the mixed vector is identical on all ranks.

## 6. Reference quantities a QE run can provide

6.1 Printed per iteration (stdout): "estimated scf accuracy < dr2 Ry" (format with 8 decimals if dr2 > 1e-8, else one digit in exponent form, so the print has limited digits: use it as a cross-check of orders of magnitude and the first digits only); "ethr = ..., avg # of iterations = ..." from c_bands (ethr used for that iteration); the total energy and its components including the double-counting term (printed when |descf| > 1e-8); "iteration # k ecut = ... beta = ..." lines. The total charge is not printed unless the starting charge needed renormalization ("starting charge X, renormalised to Y") or v_h's integrated charge differs from N_el (warning). With verbosity high more detail is printed (INFERRED, not traced).

6.2 Density file. At the end of every run (unless disk_io is none) write_scf writes <outdir>/<prefix>.save/charge-density.dat (formatted as HDF5 only if the code was built with HDF5; default is Fortran unformatted sequential). Content for nspin = 1 (INFERRED nothing else for nspin = 1, since nspin_ is 1 here): record 1 = logical gamma_only, integer ngm_g, integer nspin; record 2 = the three reciprocal basis vectors b1, b2, b3 (real double each, in 1/bohr, with the 2 pi/alat factor included: bg * tpiba); record 3 = Miller indices (3, ngm_g) integer in the global G order of section 1.4; record 4 = rho_tot(G) as ngm_g complex double in that same order (for nspin = 2 a second record holds the magnetization). The coefficients are the same rho%of_g of section 1.2 (units e/bohr^3, normalization with rho(0) = N_el/Omega), the half sphere is stored (gamma_only true, the minus-G partner is implicit). The file contains the whole ngm_g vectors (dense sphere), not just ngms. Reading back (starting_pot = 'file') reads exactly these coefficients into rho%of_g with no renormalization unless |charge - N_el| exceeds 1e-7 relative.

6.3 What the file holds (from the code order of section 4): converged run -> y^k; otherwise -> x^{k+1}. See question 7.1 for the check that decides the matter independently of reading.

6.4 One-step chain argument. A run with starting_pot = 'file', electron_maxstep = 1, mixing_beta = 1.0:
- iteration 1 uses x^1 = file density, ethr as in 2.5 (set diago_thr_init explicitly to a tight value, for instance 1e-10 or smaller, so the Davidson error does not contaminate y; also startingwfc atomic and a fixed seed or nonrandom start to remove run-to-run noise: INFERRED necessary for a bitwise comparison, not verified);
- p = 0, so x^2 = x^1 + 1.0 * (y^1 - x^1) = y^1 up to one roundoff on the sum, and ngms = ngm so the high-G branch is empty;
- not converged (unless dr2 < conv_thr, in which case the output y^1 is written directly) so rho = x^2 is written; if the first-iteration redo is triggered (dr2 < ethr*N_el) then the redo changes y only at the Davidson tolerance and the same relations hold. Either way charge-density.dat = y^1 (+ Davidson error + roundoff).
- Therefore a chain: x^1 (arbitrary density in the file format of 6.2) -> run -> file = y(x^1). Repeating from any set of densities {x_j} gives y(x_j) for each, without ever touching mixing: the map x -> y is recovered for arbitrary inputs. A history-bearing update can then be tested with an implementation fed pairs (x_j, y_j) generated this way, compared against the QE trajectory with electron_maxstep = k (file after k steps = x^{k+1}); the pairs of the QE trajectory itself need x^1..x^k (runs with electron_maxstep = 1..k - 1, files = x^2..x^k, with deterministic start) and y^1..y^k (one chain step from each x^j). The pre-requisite for this is that the QE path to x^{j+1} is reproducible between runs (question 7.2).

6.5 Quantities for term-by-term checks (all derivable from the file and the printed dr2): rho(G=0) times Omega (= N_el); dr2 recomputed from x and y with the form of 2.2 (compare to the printed first digits); the first non-zero shells (lowest |G| entries) of x^2; gamma_1 for p = 1 from A_11 = <dF,dF>, c_1 = <dF, F^2>, gamma_1 = c_1/A_11 (from two QE steps), etc.

## 7. Questions that need a dev08 run (not settled by reading)

7.1 What does charge-density.dat contain after electron_maxstep = k with a non-converged result: the mixed input x^{k+1} (code reading says yes) or y^k? Test: a run with beta = 1, electron_maxstep = 1 versus a run with beta = 0.5 from the same x^1: file difference must be exactly 0.5*(y - x) if x^{k+1} is stored, zero if y is stored. Report also what is stored when conv is reached at iteration 1.
7.2 Run-to-run reproducibility of QE: is the density file after electron_maxstep = k bit-identical between two identical 8-rank runs? Between 8 and 1 ranks (expected: roundoff-level differences from summation order and the random/atomic wavefunction start)? This sets the usable tolerance floor.
7.3 Numerical tolerance to expect between an independent implementation and QE for the update itself (given identical x, y, history as input arrays): the arithmetic is linear algebra in double; for p = 0 the difference should be roundoff (about 1e-15 relative per coefficient); for p >= 1 the error in gamma scales with the condition number of A (grows quickly near convergence since residual differences become nearly parallel; INFERRED, not measured). The tolerance is to be measured by comparing QE against QE with a different rank count and against a second solver of the p x p system, not guessed.
7.4 Whether ngms equals ngm for Si128 (look for the absence of a separate "Smooth grid" line and a "Dense grid" line only in the output; alternatively compare ngm_g in the file with the printed number of G vectors). This decides whether the high-G simple mixing branch matters at all.
7.5 Symmetry: "Sym. Ops." in the output; if more than identity, sym_rho changes y inside sum_band, so the QE reference densities contain a symmetrized y and an independent implementation must use its own y (or QE's y from the chain) before mixing. A run with nosym = .true. and noinv = .true. removes the issue.
7.6 Ordering ties: compare the Miller index list in the file with the independent generator, shell by shell.
7.7 Does the first-iteration redo trigger in the Si128 reference runs (look for the message about lowered threshold) and what is the printed ethr of the first Davidson call? Both change what the first dr2 and y^1 mean.
7.8 Printed dr2 precision limits a number comparison of the convergence metric to 8 decimals in Ry (or one significant digit below 1e-8): a higher-precision dr2 may need the instrumented route (a fork of QE on dev08 outside this repo), which is a policy decision, not done here.
7.9 Whether the density file really keeps the whole sphere (ngm_g) with the same ordering for different MPI rank counts (ig_l2g gives a rank-independent global order, INFERRED from the merge routine).

## 8. Short specification summary (for the ticket)

1. Mixed vector: valence rho(G), half-sphere coefficients including G = 0, e/bohr^3, rho(0) = N_el/Omega, no core charge, ngms entries (= all stored G when ecutrho <= 4 ecutwfc).
2. Residual F = y - x.
3. Metric <a, b> = 8 pi Omega sum_{half sphere G != 0} Re(conj(a) b)/|G|^2, in Ry (= twice the Hartree-unit expression), dr2 = <F, F>, conv when below conv_thr (Ry).
4. Next ethr = min(previous, 0.1 dr2/max(1, N_el)), floor 1e-13, 1e-2 at k = 2, first-iteration redo if dr2 < ethr*max(1, N_el).
5. History: m = 8 differences dF = F^{k-1} - F^k, dx = x^{k-1} - x^k, ring buffer (newest overwrites oldest after m), no normalization.
6. gamma = A^{-1} c with A_ij = <dF_i, dF_j>, c_i = <dF_i, F^k>.
7. x^{k+1} = x^k + beta F^k - sum_i gamma_i (dx_i + beta dF_i), beta = 0.7. First step: p = 0, x^{k+1} = x^k + beta F^k.
8. Density written at the end: y^k if converged, x^{k+1} otherwise (code order; question 7.1 confirms).
