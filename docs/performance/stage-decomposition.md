# Stage-by-stage decomposition of QE 7.6 and BAANDI (plan, 2026-11-08)

Goal (user, 2026-11-08): decompose every step of a QE run and of a BAANDI run, profile each stage on the same cases, and read the optimization points from the per-stage gap. Existing material this builds on: [qe-stage-by-stage.md](../research/qe-study/qe-stage-by-stage.md) (QE timers for si8 and si16, standard accuracy), the dev08 call-site profile of 2026-11-08 (`tests/reference/perf/dev08_materials_qe_relayed.md`), the attempts ledger. What is new: one common stage list for both codes, the same cases and accuracy classes on the same node, per-stage time AND work counters, and the gap table per stage.

## 1. Common stage list (a time in a stage is wall seconds, 1 core, pinned)
| ID | Stage | QE timer(s) | BAANDI phase / call-site bucket |
|---|---|---|---|
| S0 | Input, setup, tables, FFT grid, structure factors, pseudopotential radial tables | init_run, init_us_1, init_vloc, setup | setup_s, UPF parse, preflight, basis |
| S1 | Starting density and starting wavefunctions | atomic rho, wfcinit | first-iteration initial guess |
| S2 | H application: local potential (FFT pair work) | vloc_psi, fftw/fft | FFT + local multiply (pair plan, apply_local) |
| S3 | H application: nonlocal projections and back-application | calbec, add_vuspsi | nonlocal add / projections (hamiltonian::add, add_block) |
| S4 | Davidson: reduced matrix update | cdiaghg/rdiaghg inputs, overlap products | reduced_matrix_into |
| S5 | Davidson: small dense eigenproblem | cdiaghg / rdiaghg | evd (tridiag, divide and conquer) |
| S6 | Davidson: Ritz vectors and products (rotation) | update step products | real::rotate / rotate_pair |
| S7 | Davidson: correction, preconditioner, orthonormalization of new vectors | g_psi, MGS-like | residual, preconditioner, projection (gram, MGS) |
| S8 | Convergence bookkeeping, restarts, certificate (BAANDI strict only) | (none; QE tests Ritz changes) | driver, certificate applications |
| S9 | Density from wavefunctions (sum_band) | sum_band | density |
| S10 | Potential: Hartree and XC | v_of_rho | potential (hartree, lda_pz/pbe) |
| S11 | Mixing | mix_rho | mixing |
| S12 | Energy terms, Ewald, certificate of energies | delta_e, descf, ewald | energy audit |
| S13 | End of run: output files, forces/stress if requested | punch, forces, stress | output writers, forces, stress pass |
Work counters per stage where they exist: SCF iterations, Davidson steps, H vectors (h_psi columns), FFT calls, dgemm shapes and counts.

## 2. Cases and classes
Gamma si8, si16, si32 (QE `K_POINTS gamma` real path and BAANDI real path); k222s and k444s (symmetry on, both codes); standard class (QE conv_thr 1e-8 / BAANDI standard threshold) and strict class (QE 1e-10 / BAANDI strict certified); n as in the existing dev08 protocols. One extra cell shape (skewed si-skew) for generality. Same node, load below 0.3, same thread budget.

## 3. Measurements to collect (dev08 master, measurement only)
A. QE: full timer table of the run output (all `PWSCF` and sub-timers with calls) for each case and class; no `-nk`/`-nd` flags beyond what the earlier protocol used; record the command, version, input and UPF SHA-256.
B. BAANDI: stage shares by sampling (eu-stack, 20 ms, 3 to 10 repeats, call-site attribution as in the 2026-11-08 profile) PLUS the engine phase timers of `results.json`; MKL native and pure native builds both for the cases above; real path on.
C. Work counters: QE (iterations, `Davidson diagonalization with overlap` steps, h_psi calls) from its output; BAANDI (SCF/Davidson/H) from results.json.
D. A small microbenchmark set at the dominant shapes (single-core MKL dgemm, FFT pair transform at the grid sizes, nonlocal projection GEMM shapes) as ceilings.

## 4. Analysis the orchestrator does with the data
1. Fill the stage table: QE seconds, BAANDI seconds, ratio, share of the gap (gap_i = BAANDI_i - QE_i divided by the total gap).
2. Split each stage ratio into work ratio (counter) and rate ratio (seconds per unit of work), as in exceed-qe-plan section a.1, so that "does more work" and "does it slower" are separated.
3. Rank optimization points by gap share times plausible reduction; every candidate gets a ledger-style entry (what, expected gain, risk, rejected-idea check against the attempts ledger).
4. Re-run the table after each adopted change; the table is a standing artifact (`stage-decomposition-results.md`).

## 5. Status
Plan written; dev08 request sent 2026-11-08; no data yet.
