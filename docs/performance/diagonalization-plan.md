# Plan: making the eigensolver step faster (architect, 2026-11-07)

Inputs: docs/research/diagonalization-literature.md (sources marked OPENED, SNIPPET or MEMORY; no published comparison exists for 8 atoms, 20 bands, one core), docs/performance/profile-real-standard-si8.md, tests/reference/perf/qe_iteration_data_relayed.md, tests/reference/perf/dev08_materials_qe_relayed.md, attempts-ledger.md. Labels: M = measured on dev08 (relayed), O = Mac observation, E = estimate.

## What the evidence says

1. The solver family stays Davidson. Literature (OPENED): at 6 to 30 atoms Gamma, Davidson beats LOBPCG by 8 to 11 percent in DFTK; PPCG, Chebyshev filtering and randomized orthogonalization pay off only at hundreds to thousands of bands or on many cores. Locking-like ideas (PPCG delayed locking) are on the side of what failed in P24A.
2. Step counts are no longer the gap. With the standard profile BAANDI needs 77 Davidson steps and 772 H vectors on si8 (M, pure, T=1), QE complex Gamma at conv_thr 1e-8 needs 83 steps and about 2138/2 = 1069 band FFT pairs (M, relayed QE data; QE counts a block per step, BAANDI counts vectors, so compare steps and FFT calls, not raw "ops"). So the remaining gap on dev08 (QE real Gamma 0.367 s versus BAANDI standard 1.32 s at equal accuracy class, M) is cost per H application and per algebra step, plus SCF overhead already attacked in P28A/B, not the number of iterations.
3. After the speedups already merged (real Gamma, standard profile, setup, audit) the diagonalization is about 70 percent of si8 wall on the Mac (O). Its composition (FFT, nonlocal, subspace GEMM and eigen, residual and correction work) has NOT been profiled on dev08 with at least 30 samples for the real+standard configuration. That profile is the prerequisite for everything below.

## Phases

### D0. Measure (no code; dev08, through the orchestrator)
- si8 Gamma, 1 thread, pure native, gamma_real on, standard profile, 30 runs: share of wall in (a) local FFT pair applications, (b) nonlocal projection GEMMs, (c) subspace algebra (reduced matrix updates, small eigen solve, Ritz rotation), (d) residual and correction construction, (e) everything else; plus the same for k444s with symmetry on. Include the QE side for the same case from its timers (h_psi, vloc_psi, fftw, cegterg sub-timers, already relayed for the complex path; ask for the real-gamma path timers).
- Extract the SCF iteration history of si8 MKL standard (27 iterations versus 15 on pure): open anomaly, requested from the dev08 session.
- Gate to continue: the D0 profile decides the order of D1 to D3; if FFT plus nonlocal is below 60 percent, algebra work moves up.

### D1. Cheap switchable experiments (hours each, default off, generality set required)
- D1a capacity: Davidson subspace capacity 2x and 3x the band count for 20 and more bands (backlog B17; A08/A13 were tested at 8 bands only). Expected 5 to 15 percent (E). Acceptance: equal delivered accuracy, standard and strict.
- D1b preconditioner: Teter-Payne-Allan function with per-band scale and a Ritz-shifted denominator, as an alternative to the kinetic diagonal with floor 0.1 (A14 added potential diagonals; this is different and the ticket must say so). Expected 0 to 15 percent fewer steps (E).
- D1c pair density accumulation for the real path (R2 routine exists): about 3 percent of wall (E, profile-based).
- D1d real-path algebra: real GEMM tile and threshold tuning, MKL dgemm/dsyevd adapters for the real kernels (design R8) only if D0 shows algebra above 25 percent.

### D2. Per-application cost (after D0)
- If the local FFT dominates: compare our pruned fused FFT against MKL DFTI per band pair at the si8 grid on dev08 (the plan R2 microbenchmark requested), then consider batching bands per plane (A28 rejected band-grouped FFT in the complex path; new evidence is the pair structure, say so) and 2D-then-1D pruning variants.
- If nonlocal dominates: band-blocked projector GEMMs (P15B measured on dev08 once master runs (6); A43 row-blocked traversal lost, P15B is the bit-identical complex alternative and the real operator already uses GEMMs).
- Threads: 8 to 16 core scaling with the larger benchmark cases (si32 and si16-k222, P25A) before any claim for large systems.

### D3. Rounding-changing options (need a user decision, not started)
- fp32 H application and projections only in loose SCF steps (tolerance at or above about 1e-5), fp64 for the final steps and the certificate. Literature gain over 30 percent in the fp32 part (SNIPPET, OPENED abstracts), our estimate 5 to 15 percent of wall (E, low-medium confidence), high implementation cost. It changes rounding, so it is a decision for the user; start only after D0 shows the share of H applications in loose steps.
- Skipping the energy audit in non-converged standard iterations (profile C7): 8.7 to 18 percent, changes the output history. User decision.

### Deferred (isolated modules for larger systems, per the generality rule)
Chebyshev-filtered iteration (ChASE-like) for many cores or GPU, PPCG for 1000+ bands, randomized orthogonalization above about 2e4 plane waves, wavefunction extrapolation between ionic steps (helps relaxation, 20 to 40 percent fewer SCF iterations per step (E), needs forces/relaxation hooks).

## Expected effect (E, honest)
D1a to D1c together 10 to 25 percent of the diagonalization if they stack, which is 7 to 18 percent of wall: useful but not the missing factor 2 to 2.5 against QE real Gamma. The larger lever is per-application cost (D2), whose size is unknown until D0 is measured. No technique here claims to close the gap alone.

## Rules kept
No tolerance relaxed; every option default off until dev08 equal-accuracy evidence; rejected ideas (A28, A33, A35, A43, A52, P24A) are not repeated without the new evidence stated in the ticket; Mac numbers are observations; at least 30 samples for named-site claims.
