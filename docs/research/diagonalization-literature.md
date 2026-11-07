# Literature survey: faster iterative diagonalization in plane-wave Kohn-Sham SCF

Written 2026-11-07 (research note, no code changed, QE source not read). Audience: the architect planning the next eigensolver work. Context: si8 Gamma, 20 bands, one core, standard profile, real-Gamma path: wall 0.615 s, of which diagonalization 0.339 s (55%), energy audit 0.054, density 0.047, potential 0.017, mixing 0.014, engine setup 0.127 (`docs/performance/profile-real-standard-si8.md`, Mac observation). Davidson already warm-starts from the previous SCF step's states (`scf.rs`, `initial.map(...)`).

## 0. Source status (read this first)

- **OPENED** = the text or abstract page was actually fetched and read in this session; numbers below come from it.
- **SNIPPET** = only a search-result summary was seen; treat the number as second hand.
- **MEMORY** = from the author's general knowledge, unverified.
- Several PDFs came back as unreadable binary (Levitt-Torrent PDF, LOBPCG-chem PDF, ChASE PDF, ABINIT-GPU PDF). Their HTML or abstract versions were read instead; where only the abstract was read, no timing claim is made.
- Several 2026 arXiv identifiers (2608.24529, 2604.11139, 2609.21759) are newer than the author's training data; they were opened and are cited as found.
- Honest gap: no source found reports an eigensolver comparison for a 8-atom, 20-band, single-core case. Nearly all published gains are for hundreds to thousands of bands, or for parallel scaling. That is the main finding for BAANDI's current regime.

## 1. Ranked table (for BAANDI as it is today: small to medium, 1 core to 1 socket, equal accuracy)

Gain = expected share of si8-class wall time or step count, an estimate by the author (E), not a measurement. Conf = confidence that the gain is real on BAANDI. Ledger = conflicts with a row (premise must be restated).

| Rank | Technique | What it changes | Reported gain (source status) | Est. gain for BAANDI | Conf | Difficulty | Accuracy risk | Ledger conflict |
|---|---|---|---|---|---|---|---|---|
| 1 | fp32 H application and projections only while the eigen tolerance is loose (early SCF steps), fp64 for the tail | Precision of FFT, GEMM, nonlocal in steps where tol >= ~1e-5 | fp32 diagonalization part >30% faster at fp64 final accuracy (Tsuchida, Choe, CPC 2012, SNIPPET); sketch in fp32, H and residual in fp64, accuracy kept (Gubler et al., arXiv 2608.24529, OPENED); R-ChFSI tolerates FP32/TF32 matvecs, 1e-8 residual met (Kodali et al., arXiv 2503.22652, abstract OPENED) | 5-15% of wall if 30-40% of H applications are in loose steps and f32 is ~1.5-1.8x on the Mac/Xeon | low-med | High (generic f32 path through rustfft, faer, projector code, pair trick) | Low if the last steps and the certificate stay fp64; a rounding change needing the P14B approval rule | none (A24 MKL, A35 not about precision) |
| 2 | Per-band residual test only (no all-root Ritz and residual work) via PPCG-style periodic Rayleigh-Ritz | Rayleigh-Ritz every 5-10 iterations on small block problems instead of all-root RR every step | PPCG 1.5-1.7x over QE Davidson (subspace cap 2k) at 2000-2550 bands, Gamma, LDA; RR was 66 of 83 s in Davidson on Li318 (Vecharynski, Yang, Pask, J. Comput. Phys. 290, 73, 2015, arXiv 1407.7506, OPENED) | ~0 at 20 bands (RR is O(j^3), j <= 120, already cheap); real only above ~200 bands | high that it does not help si8 | Medium-high (new solver, sequential inner loop) | Low (same convergence test on residual) | related to A56 gated, B01/P24A (hard locking failed); PPCG used delayed locking and got more A.X (OPENED) |
| 3 | Preconditioner: TPA function and a Davidson-style shifted denominator, per-band scaling | Step count of Davidson | TPA is the standard plane-wave preconditioner (Zhou, Chelikowsky, Gao, Zhou, CiCP 18, 167, 2015, abstract OPENED; formula from search page SNIPPET); VASP uses a similar damped form (MEMORY); LOBPCG vs Davidson with diagonal vs iLU in chemistry: 56/53/47 vs 52/54/50 iterations (Herbst et al., TCA 2023, OPENED) | 0-15% fewer steps, unmeasured; no plane-wave iteration table was found | low-med | Low (a few lines, fully switchable) | None | A14 (fuller diagonal) rejected: the exact mechanism (adding V and nonlocal diagonal) differs from TPA smoothing and from shifting by theta. State this in the ticket |
| 4 | Wavefunction (and density) extrapolation between ionic steps | Initial guess across relaxation or MD steps | "greatly reduces SCF iterations", optimal order exists and grows with time step and tighter tolerance (Fang, Gao, Song, Wang, arXiv 1604.04712, abstract OPENED; no number); projected commutator-DIIS extrapolating density matrix, SNIPPET | 20-40% fewer SCF iterations per relaxation step (E, MEMORY-level) | med | Medium (subspace alignment/Procrustes between steps, needs P20A-P23A hooks) | None (initial guess only) | none; does not help a single-point SCF |
| 5 | Capacity and block experiment at 20 bands (2x to 3x bands vs 6x) | Basis size and dense algebra per step | QE caps the Davidson subspace at 2k; PPCG paper compares against exactly that (OPENED) | 5-15% at 20 bands (E) | low-med | Trivial (option) | None | A08/A13 rejected at 8 bands only; backlog B17 re-opens it for 20+ bands. Needs the evidence in the ticket |
| 6 | Density by real pairs (2 bands per complex FFT), already in H | Density phase | Gamma-only codes "at least twice as fast", half the plane waves (QE slide deck, VASP manual page, both SNIPPET) | density 0.047 s -> ~0.025 s, about 3-4% of wall | high | Low-med | none (rounding at 1e-16) | none; this is `profile-real-standard-si8.md` section 2 follow-up |
| 7 | Randomized (sketched) Gram-Schmidt in Davidson | Orthogonalization cost | 5-11% full-SCF gain over LOBPCG, and sketched GD only ~equal to deterministic GD, 6-30 atoms Gamma, DFTK (Gubler et al., arXiv 2608.24529, OPENED) | negative at this size: crossover near n ~ 2e4 for k = 500 | high that it does not pay | High | Low | related A52 (CholeskyQR2 rejected on conditioning) |
| 8 | Chebyshev filtered subspace iteration (CheFSI, ChASE) | Replaces Davidson by polynomial filter, one RR per SCF step | Wins at scale: 256-atom Ti, ~15 min on 4096 cores vs >1 h LOBPCG on 1024 (Levitt, Torrent, arXiv 1406.4350, OPENED); on 100 bands CG and full-block LOBPCG still competitive, at 200 bands CG clearly better, performance degrades about 1/sqrt(Ecut), no preconditioner (same paper) | negative for si8, 1 core | high | Medium-high | Low | none; A22 (RMM) showed the same pattern "needs 2-10x more H columns" |
| 9 | LOBPCG instead of Davidson | Solver choice | Davidson "remains the best choice for most applications", LOBPCG 20-30% slower when memory is ample, wins on memory and as a fallback (Herbst et al., TCA 2023, arXiv 2305.06668, OPENED, quantum-chemistry matrices, not plane waves) | negative | high | Medium | none | A22 pattern |
| 10 | Neural initial guess for density | SCF count | up to 25% end-to-end wall reduction on unseen PAW structures, perfect initialization would give 40-52% (arXiv 2609.21759, abstract OPENED) | not applicable to the eigensolver; SCF count ~14, atomic-density start already good | low | Very high | none | none |

Not ranked because no usable source or no fit: Jacobi-Davidson, Davidson-Liu variants (MEMORY only: correction-equation solves rarely pay in plane-wave codes because the diagonal preconditioner is already strong), double-grid for local potential (see 5.4), communication-avoiding algebra (no multi-node target).

## 2. Block eigensolver variants (task item 1)

**Davidson (what BAANDI has).** Opened evidence is consistent: in chemistry matrices Davidson is optimal when memory allows and LOBPCG is 20-30% slower in time but uses 356 GB vs 55 GB on the largest CI case (Herbst et al. 2023). In plane-wave SCF with 6-30 atoms Gamma, generalized Davidson is 8-11% faster than DFTK's LOBPCG (Gubler et al. 2026). So the eigensolver that BAANDI uses is already the one that wins in small and medium plane-wave runs. No published argument was found for replacing it at 20 bands.

**LOBPCG.** Three-block basis (X, W, P), RR of size 3N each iteration, O(N^3) with N = bands; instabilities from ill-conditioned bases are fixed by double orthogonalization, shifted Cholesky and locking (Herbst et al. 2023). Multiblock splits the RR by N_blocks^3 (Levitt-Torrent). Relevance: BAANDI already does the same repairs inside Davidson (A52 failure is the same conditioning issue).

**PPCG.** Small 3 x sbsize problems, RR only every 5-10 iterations, Cholesky QR in between. Measured vs QE Davidson (subspace cap 2k), Gamma, norm-conserving LDA, Edison: Li318 480 cores 49 s vs 84 s (1.7x), graphene512 97 vs 144 (1.5x), bulk Si 2550 bands 189 vs 329 (1.7x); tolerances 1e-2 to 1e-3 (differ per system); the authors are the algorithm's developers; orthogonalization every 15 iterations failed to converge once; sbsize 1 loses BLAS3; PPCG spent more in A.X (delayed locking) while saving RR (OPENED). The gain scales with the RR share, so it belongs to the thousands-of-bands regime. A Feng et al. 2024 study compares Davidson, LOBPCG, PPCG and CheFSI in PWDFT up to 42,592 cores (CPC 299, 109135; abstract only via SNIPPET, comparison numbers NOT seen).

**CheFSI and ChASE.** First SCF step needs a real diagonalization, later steps a degree-m filter on the warm subspace (Zhou, Saad, Tiago, Chelikowsky 2006, SNIPPET: "several times faster" than the best eigensolvers in real-space DFT). Plane-wave evidence (Levitt-Torrent, OPENED): Chebyshev beats CG near the bottom of the spectrum with 100 bands but loses at 200; degree n costs n H applications per band per outer step without orthogonalization inside the filter; RR is one per SCF step; scalability is far better than LOBPCG (still scaling at 16384 processors where LOBPCG saturates at 2048). ABINIT GPU paper (arXiv 2604.11139, abstract OPENED, body unreadable): on CPUs higher Chebyshev degree does not pay back, on GPUs it does (statement from a search summary, SNIPPET). ChASE targets sequences of dense Hermitian problems, filter degree optimized to minimize flops, wins when about 20% or less of the spectrum is needed (Winkelmann, Springer, Di Napoli, arXiv 1805.10121, abstract OPENED; helmholtz.software SNIPPET). R-ChFSI makes Chebyshev filtering tolerant to inexact matvecs: GPU speedups up to 2.7x (filter) and 2.1x (whole solve) for finite-element DFT with FP32/TF32 (arXiv 2503.22652, abstract OPENED). Verdict for BAANDI: a parallel-scaling tool for later (many cores, thousands of bands, GPU), not a single-core win.

**RMM-DIIS (VASP, QE).** Kresse and Furthmuller (PRB 54, 11169, 1996, abstract OPENED): DIIS diagonalization with minimal O(N^3) work, iteration count nearly independent of system size, roughly N^2 scaling up to ~1000 electrons. Mechanism details (one orthonormalization per sweep, subspace rotation each SCF step, Davidson for first iterations) are MEMORY. Repository evidence: A22 hybrid RMM needs 2-10x more H columns than Davidson at 8-64 bands; A23 shows cold starts fail. Matches the literature statement that RMM-DIIS pays only for many bands.

**Jacobi-Davidson.** Not opened; not ranked.

**When each wins (synthesis).** Few hundred bands or fewer, one node: block Davidson. Memory-bound: LOBPCG. Thousands of bands: PPCG, RMM-DIIS, CheFSI (RR cost dominates). Tens of thousands of cores or GPUs: CheFSI/ChASE (few global reductions). Metals: more bands and partial occupations; Chebyshev and Davidson both need the buffer-band logic that BAANDI already has (A17). Many k-points: embarrassingly parallel over k, no change in the per-k solver choice.

## 3. Preconditioners (task item 2)

- TPA (Teter, Payne, Allan, PRB 40, 12255, 1989; MEMORY plus citation confirmed by SNIPPET): K(x) = (27 + 18x + 12x^2 + 8x^3) / (27 + 18x + 12x^2 + 8x^3 + 16x^4), x = (kinetic energy of the G component)/(kinetic energy scale of the residual); K -> 1 at low x, -> 1/x^... (inverse kinetic energy) at high x. The scale is chosen per band, for example the mean kinetic energy of the current band. A general family with higher-order approximation was derived by Zhou, Chelikowsky, Gao, Zhou (CiCP 18(1), 167-179, 2015, abstract OPENED; no iteration counts in the abstract).
- Kresse-Furthmuller use an analogous smooth kinetic damping, with a metric tailored to plane waves for mixing (abstract OPENED; formula MEMORY).
- Lu and Yang, preconditioned orbital minimization, a Fermi-operator preconditioner, iteration number reduced to O(1) in the reported tests (SIAM Multiscale Model. Simul. 15, 254, 2017, arXiv 1603.09724, SNIPPET): cost and invasiveness too high for BAANDI now.
- Chemistry comparison: diagonal vs tridiagonal vs iLU preconditioner changes iteration counts only 47-56 on a hard Hessian (OPENED); diagonal is the practical compromise.
- Effect for BAANDI: unknown without measurement; none of the opened sources gives a plane-wave table of iteration counts per preconditioner. A14 shows adding the local potential mean and nonlocal diagonal gave no gain. Two untried variants: (a) TPA smooth form with per-band scale; (b) denominator shifted by the current Ritz value (Davidson-style `1/(h_diag - theta)` with a positive floor, MEMORY of standard Davidson practice). Experiment cost is small; judge by Davidson steps and vectors on the benchmark suite, not Si alone.

## 4. Fewer H applications per SCF step (task item 3)

- **Warm start between SCF steps.** Already adopted (previous states as initial, A06 and A07 on the first step). Chebyshev-type methods and ChASE depend on this alignment (ChASE abstract, OPENED).
- **Adaptive eigen tolerance.** BAANDI has `min(0.1 * density residual, 1e-2)` (A06); factors of 3 to 100 on the density-residual term were rejected (A21); a Ritz-change stopping rule is in the standard profile (P17A/P15A). Literature: DFTK adapts the number of bands during SCF (release note, SNIPPET, not the tolerance); Herbst and Levitt adapt the SCF damping by line search (arXiv 2109.14018, abstract OPENED; no H-application claim). No published inexact-Newton style schedule for plane-wave Davidson was found that beats the schedules already tried. Treat this family as exhausted without new evidence.
- **Early stopping.** A15 rejected the strict-refinement variant; the Ritz-change standard profile adopts the cheaper variant and cut si8 from 345 to 69 steps (Mac observation, uncertified, from backlog). The remaining uncertainty is the equal-accuracy definition, not speed.
- **Mixing that exploits eigensolver state.** Pulay and Kerker variants are in `scf-mixing-strategies.md` and A31; projected commutator-DIIS on density matrices reduces SCF iterations in hybrid-functional plane-wave runs (SNIPPET, hybrid only). Not a Davidson cost lever for LDA.
- **Wavefunction extrapolation.** Helps ionic steps (rank 4), not the SCF loop itself, because inside SCF the states are already close.
- **ML guesses.** Density-based neural initializers reduce wall by up to ~25% in PAW (OPENED abstract) and ~33% fewer SCF iterations in molecular codes (SNIPPET); none addresses Davidson. Credible but outside the current scope.

## 5. Dense algebra (task item 4)

1. **Avoiding all-root Rayleigh-Ritz each step.** Both PPCG (OPENED) and BAANDI's own accounting (A55, A56: Ritz/residual is 56.8% of explicit GEMM work on si8) identify it. Locking failed (P24A), so the PPCG variant (periodic RR plus delayed locking) is the literature's answer, but it costs more H applications (OPENED), exactly the failure side of P24A (+24% H applications). Gain only likely at larger band counts.
2. **Cholesky-based orthonormalization.** Needs double application and shifted Cholesky above condition ~1e8 (Herbst et al. 2023, OPENED). A52 rejected CholeskyQR2 on skewed cells; the literature repairs (shifted CholeskyQR3) are MEMORY-level extra and would change rounding again.
3. **Randomized orthogonalization.** 14njk -> 10njk flops claimed (OPENED), with a generalized eigenproblem H y = theta G y of size j; fp32 sketch; full-SCF gain 5% (GPU) to 10% (CPU) over LOBPCG, similar to plain Davidson; loses at small n. Rank 7.
4. **Mixed precision.** See rank 1. Bound on usefulness: the reduced matrices are tiny; the gain is in the n x j products and FFT. Tsuchida-Choe (CPC 2012) report >30% (SNIPPET, direct-minimization eigensolver). R-ChFSI shows the residual-based formulation is what keeps convergence under inexact products (abstract OPENED).
5. **Why nvecx = 2 x bands in QE.** Not published as a rule in the sources opened; the PPCG paper only states QE's subspace is capped at 2k (OPENED). Repository notes (`qe-study/qe-solvers-and-subspace.md`) hold the analysis. For BAANDI the cap is 6 x bands (backlog B17), a candidate experiment (rank 5).
6. **Crossover GEMM vs FFT.** No published threshold found. In the PPCG tables A.X is about 94-96 s of 190-330 s at 2550 bands while RR plus QR is 59-191 s (OPENED), so algebra overtakes H at thousands of bands; BAANDI's own si8 numbers (profile above) show H plus algebra both significant at 20 bands.

## 6. Fast H application (task item 5)

- **Gamma trick.** Two real bands in one complex FFT, half the plane waves stored; QE materials say "at least twice as fast", VASP exposes a real-to-complex variant, and DFTK has an open issue about it (all SNIPPET). BAANDI has it in H (R2/R4) but the density phase does not use it (`profile-real-standard-si8.md`): rank 6.
- **Pruned FFT.** Adopted (A02, A26); QE's DFTI path does not prune (backlog note 4). No new literature lever.
- **Batched FFT / band-grouped.** A28 and A11 rejected; GPU papers (ABINIT, arXiv 2604.11139, abstract OPENED) batch because kernel-launch overhead is large; that premise does not hold on a CPU core.
- **Nonlocal projector tricks.** Real-space projectors (real-space implementation of nonlocal pseudopotentials, titled in search results, NOT opened) trade accuracy (ghost forces, egg-box) for cost; they violate equal delivered accuracy unless the radius is converged. Reciprocal-space blocked traversal A43 rejected. Not recommended.
- **Double grid / smaller grid for V psi.** MEMORY analysis: the density sphere is 2 Gmax of the wavefunction sphere, and the product V psi needs a grid of about twice that to avoid aliasing into the retained sphere; a coarser grid changes H and thus the answer. Do not pursue without an exactness proof.

## 7. Recent results by group (task item 6)

- **DFTK:** Randomized and deterministic Davidson (GD) 8-11% faster than LOBPCG in full SCF; Hamiltonian application the largest cost and equal across solvers (arXiv 2608.24529, OPENED). LOBPCG robust implementation paper (Herbst, TCA 2023, OPENED). Adaptive damping (arXiv 2109.14018, OPENED).
- **ABINIT:** LOBPCG and Chebyshev on GPUs, comparison by efficiency per H application; ChFSI degree trade-off differs between CPU and GPU (arXiv 2604.11139 abstract OPENED, details from SNIPPET). Chebyshev origins in Levitt-Torrent (OPENED).
- **PWDFT/DGDFT, SIRIUS, Qbox, JDFTx, CASTEP, VASP:** Only PWDFT (Feng et al. 2024) was found with a four-solver comparison, abstract level. No eigensolver-efficiency papers from the others were found in this search; "not found", not "does not exist".
- **Finite-element and real-space DFT (not plane waves):** R-ChFSI (OPENED abstract), ChFSI robustness (SNIPPET), cited for the precision idea only.
- **Quantum ESPRESSO developers:** PPCG study used QE as baseline (OPENED); a QE 7.x ParO/PPCG implementation exists per a SISSA slide and the repository's own notes (SNIPPET).

## 8. Recommendations for the architect

1. Do not change the solver family. Davidson is the published winner for 6-30 atom Gamma plane-wave SCF (rank 9, rank 7, rank 8 negative).
2. Cheap switchable experiments first (hours each): TPA preconditioner and Ritz-shifted denominator (rank 3); capacity 2x/3x at 20 bands (rank 5); density pair trick (rank 6). Each needs the P10G generality set; A14 and A13 require the premise statements written above.
3. Medium bet with measurable upside: fp32 in loose SCF steps (rank 1). First measure fp32 FFT and GEMM speedup on dev08 and the share of H applications in steps with tol >= 1e-5 (from `P14A_WORK_ACCOUNTING`), then decide. Final fp64 steps and certificate keep delivered accuracy.
4. Defer to larger systems or relaxation: PPCG (RR share), CheFSI/ChASE (scaling, GPU), wavefunction extrapolation (ionic steps) as isolated switchable modules per the Generality rule.
5. Largest remaining non-eigensolver items in the si8 profile (engine setup 20%, energy audit 9%) are already documented in `profile-real-standard-si8.md` and should be handled before sophisticated eigensolver work.

## 9. Sources (OPENED unless marked)

- Levitt, Torrent, Parallel eigensolvers in plane-wave DFT, arXiv 1406.4350 (ar5iv HTML).
- Vecharynski, Yang, Pask, PPCG, J. Comput. Phys. 290, 73 (2015), arXiv 1407.7506 (ar5iv HTML; last 26k characters not read).
- Herbst et al., open-source LOBPCG, Theor. Chem. Acc. 2023, arXiv 2305.06668 (ar5iv HTML).
- Gubler, Park, Bussy, Laute, Marzari, Timrov, Grigori, Randomized block Davidson for PW DFT, arXiv 2608.24529 (HTML).
- Kodali, Ramakrishnan, Motamarri, R-ChFSI, arXiv 2503.22652 (abstract page).
- Winkelmann, Springer, Di Napoli, ChASE, arXiv 1805.10121 (abstract page).
- Kresse, Furthmuller, PRB 54, 11169 (1996) (publisher abstract via materialsdesign.com).
- Zhou, Chelikowsky, Gao, Zhou, CiCP 18, 167 (2015) (Cambridge abstract page).
- Fang, Gao, Song, Wang, arXiv 1604.04712 (abstract page).
- Herbst, Levitt, line search for SCF, arXiv 2109.14018 (abstract page).
- Lygatsika et al., ABINIT GPU, arXiv 2604.11139 (abstract page).
- Ærtebjerg, Elsborg, Bhowmik, neural initialization, arXiv 2609.21759 (abstract page).
- SNIPPET only: Zhou, Saad, Tiago, Chelikowsky, J. Comput. Phys. 219, 172 (2006) and cond-mat/0703239; Tsuchida, Choe, Comput. Phys. Commun. (2012), arXiv 1108.4509; Feng et al., Comput. Phys. Commun. 299, 109135 (2024); Lu, Yang, arXiv 1603.09724; QE slide deck on Gamma tricks; VASP manual page; DFTK issue 36; helmholtz.software/software/chase.
- Internal: docs/performance/attempts-ledger.md rows A08, A13, A14, A15, A21, A22, A23, A28, A43, A52, A55, A56 and P24A; backlog B17; profile-real-standard-si8.md.
