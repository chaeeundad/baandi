# Performance backlog: ideas not yet tried, or only partly tried

Compiled 2026-11-05. Ranked by (justified gain x chance it works), highest first. Ledger row IDs (`A..`) refer to [attempts-ledger.md](attempts-ledger.md); read the rows named under "builds on / avoids" before starting. Every gain below is an ESTIMATE or a CEILING unless stated as measured; the basis is given in each item. Acceptance is always dev08 via the orchestrator ([AGENT-BRIEF.md](AGENT-BRIEF.md)). Labels: **ready** (can be implemented now, the first experiment is already known to be cheap), **needs-measurement** (a cheap measurement must justify the work), **blocked** (needs a decision, a case, or a prerequisite that does not exist).

Common gates for every item (details in the brief): energy within 1e-12 Ha of the stored baseline on gamma, k222s, k222, k444s, si8 plus the generality set (cubic, skewed, elongated, slab; 4/8/20 atoms; synthetic projectors); SCF count and certificate unchanged for realistic fixtures; density bit-identical across thread counts within a backend; no tolerance relaxed; no new convergence failure.

## Ranking

| Rank | ID | Idea | Label | Gain basis |
|---|---|---|---|---|
| 1 | B01 | Cheaper tail steps: restrict Ritz, residual and reduced-problem work to unconverged and guard roots | needs-measurement | ceiling 57% of remaining si8 Davidson GEMM work |
| 2 | B02 | Real-valued wavefunctions at Gamma (two bands per FFT, half plane waves, real GEMM) | needs-measurement | estimate 1.7-2.3x on Gamma-only runs |
| 3 | B03 | Band-level parallel decomposition with fewer, larger regions | needs-measurement | utilisation 23% vs 33% needed for 2x |
| 4 | B04 | Fewer Davidson steps through an energy-error stopping test and the energy-weighted SCF threshold (H6) | needs-measurement | step ratio 3.5x vs QE (ceiling) |
| 5 | B05 | Faster convergence of the slow tail root (previous-direction or doubled correction vectors) | needs-measurement | 55% of si8 steps in 5 of 15 SCF iterations |
| 6 | B06 | Serial remainder at Gamma: parallel potential, energy audit, mixing; faster EVD | ready (audit/potential), needs-measurement (EVD) | about 5-10% of si8 T=8 wall |
| 7 | B07 | `NonlocalOperator::add` as GEMM (projector application) | needs-measurement | 6-7.5% of profile |
| 8 | B08 | FFT planning and layout (batched stick transforms with MKL, plane buffers) | needs-measurement | FFT about 18% (stale) |
| 9 | B09 | k-point parallel efficiency (dynamic ordering, per-k cost model) | needs-measurement | 7.1x of 16 on k444s |
| 10 | B10 | Cheaper first SCF iteration on many-k (better starting vectors) | needs-measurement | 12.6% of k444s H vectors |
| 11 | B11 | RMM-DIIS at 128 or more bands in a real SCF | blocked | no real case above 64 bands |
| 12 | B12 | MPI (and later GPU) | blocked | out of scope per AGENTS.md |

## Items

### B01 Cheaper tail steps (locking or restricted Ritz/residual work)
- **Expected gain.** CEILING: after P14B (A55) the all-root Ritz/residual products are 5386241040 of 9489226838 2mnk (56.8%) on si8 and 946782528 of 2258007134 (41.9%) on k444s (N4, p14b/adopted/acceptance.md). Tail steps carry few vectors (si8 SCF 12: 52 steps for 192 vectors; SCF 4, 9, 10, 12, 13 hold 188 of 340 steps with 3.7-6.4 vectors per step, p14a work-tables) yet pay all-root work. ESTIMATE of attainable gain: if GEMM is about half of si8 T=1 wall (derived from 1.78 s summed GEMM time vs about 3.7 s wall, different runs) and 60% of the ceiling is removable in tail steps, about 8-15% of wall. Not measured.
- **Risk.** Under the unlocked stopping rule the guaranteed removable work is zero because a converged root can fail after the next rotation (A56, N4). A lock rule changes the trajectory (rounding) and risks the fixed-point multiplicity problem (ledger C.11); the final certificate must still use true residuals for all roots.
- **Prerequisites.** P14A collector (exists, source-only) to count re-crossings; architect decision on a rounding policy like P14B.
- **Validate.** Gates above plus: certificate unchanged, per-root residual history compared to baseline, 20-atom skewed and elongated fixtures (the cases where P14A/P10M broke).
- **Builds on / avoids.** Builds on A55, A56, A53. Avoids the exact idea in A15 (Ritz-value-change test with a strict refinement pass) and the CholeskyQR2 route (A52).
- **First experiment.** From the existing P14A collector on si8, k444s, skewed 20-atom and elongated 4-atom: count how often a root falls below tolerance and later exceeds it again (the table only records permanence). If zero re-crossings occur when roots are locked at 0.1 x tolerance, prototype "recompute residual of locked roots every n steps and at certificate".

### B02 Real-valued Gamma path (two bands per FFT, half plane waves, real GEMM)
- **Expected gain.** ESTIMATE, labelled: for Gamma-only runs the wavefunction coefficients satisfy c(-G) = c(G)*, so half the plane waves are independent (the engine uses the full complex sphere by design, docs/physics-conventions.md, N2 section 3). Two real bands share one complex FFT (QE note, N2 item 3: the QE comparison inputs used complex Gamma "not the real-wavefunction gamma trick", tests/reference/perf/qe_iteration_data_relayed.md). Using the stale Mac pre-P10I profile (FFT about 18%, nonlocal and local about 10%, GEMM and projection about 55%): FFT and nonlocal work halves, real GEMM costs about a quarter to a half of zgemm per independent coefficient; total time factor about 0.43-0.58, i.e. 1.7-2.3x for Gamma runs. The profile is stale (before P10I to P14B), so re-profile first. No gain for k-point runs.
- **Risk.** High. New wavefunction representation and a second Hamiltonian, density, Davidson and orthonormalisation path (real arithmetic, G=0 handling, normalisation factor 2 for paired coefficients); needs an ADR (physics conventions) and a QE complex-Gamma validation; not bit-identical to the complex path; hits the fixed-point multiplicity problem when comparing trajectories. Only Gamma-only calculations benefit; QE's complex-Gamma ratios in README stay valid only if the README states which path was timed.
- **Prerequisites.** ADR; ticket that keeps the complex path as the reference oracle; real-valued GEMM adapter (faer, MKL dgemm).
- **Validate.** Energy within 1e-12 Ha of the complex path on gamma and si8 and the 4/8/20-atom synthetic Gamma fixtures; eigenvalues and density against the stored QE complex-Gamma comparisons; determinism across T; certificate retained.
- **Builds on / avoids.** Builds on A26 (fused FFT, zero-weight skip), A36, A55. Does not conflict with rejected rows; it is not in the ledger.
- **First experiment.** Time the two-bands-per-complex-FFT trick in the existing fused operator at 36^3 on dev08 (microbenchmark through `FusedOperator`), and re-profile si8 T=1 on dev08 for the FFT, nonlocal and GEMM shares. If FFT plus nonlocal is below 25% of T=1 wall, rank drops.

### B03 Band-level parallel decomposition with fewer, larger regions
- **Expected gain.** Basis measured: worker utilisation 23.0% (mean 1.84 of 8 busy), serial 9.2%, ideal 4.24x vs actual 1.40x at si8 T=8 (D08 part 4); the P10L ticket states 33% utilisation is needed for 2x and 55% for 3x. All scheduling-only attempts since P10I gave about 5% each on the Mac (A40, A45, A48). A structurally different decomposition (parallelise over the 20 bands with per-band serial kernels for residual, preconditioning, correction-projection, and fuse small phases into one region) has not been tried; estimate of gain unknown, 1.3-2x at T=8 would meet the P10L targets. Estimate only.
- **Risk.** Changes GEMM shapes and so rounding and the Davidson trajectory (A49 shows finer tiles changed Davidson work 408 to 362). Per-band kernels are many small GEMV-like operations; may lose single-thread efficiency (A35 was 6% slower at T=1). Needs an accepted rounding policy.
- **Prerequisites.** A dev08 profile of si8 T=8 attributing idle time per phase (A48 counters exist: `worker-analysis`), architect decision on rounding.
- **Validate.** si8 T=8 dev08 scaling at least 2x (T1/T8) with T=1 not slower by more than 1%; 30-pair gate for T=1 on si8 and k444s; the 20-atom skewed synthetic for iteration counts.
- **Builds on / avoids.** Builds on A40 (pool), A45, A48 counters. Avoids A33/A34 (scoped threads, packed-lane FFT threads), A46 (capped leaves), A49 (finer GEMM tiles), A50.
- **First experiment.** Use the existing `worker-analysis` build and the P10L region CSV (Mac) plus a dev08 run: list, per Davidson iteration, the serial chain length (MGS, EVD, residual norms, preconditioner) versus the parallel work. If the serial chain is above 40% of Davidson time at T=8 (P10L: 42.2% outside scheduled regions on the Mac), prototype the per-band region for the residual/precondition/projection chain only.

### B04 Fewer Davidson steps: energy-error stopping test and H6
- **Expected gain.** CEILING: QE needs 83 steps vs 290-340 for BAANDI on si8 (3.5x), accounting for 3.49 of the 8.48 observed time ratio (ledger D). Attainable gain unknown. The previous attempts failed for specific reasons: the L2-factor schedule raises SCF count (A21), Ritz-change test with strict refinement erases its saving (A15), preconditioner and capacity changes did not help (A13, A14).
- **Risk.** QE's threshold is an energy-error estimate (Hartree-metric density residual divided by electron count, with a first-iteration redo), our test is a residual norm in Ha (N1 section 3). A residual r gives eigenvalue error about r^2/gap but wavefunction error linear in r; density accuracy needs the latter. A looser eigenvalue test may therefore raise SCF count again (A21). The first-step redo adds one extra diagonalisation in the worst case.
- **Prerequisites.** P14A collector for per-iteration steps; none else.
- **Validate.** Gates plus: SCF count unchanged on si8, k444s, skewed and elongated fixtures; total H vectors and Davidson steps reduced; fixed-point certificate unchanged.
- **Builds on / avoids.** Builds on A21 (does not repeat the multiplied-residual family), A53. Avoids A15 (do not add a strict refinement pass; certificate only at the end), A14.
- **First experiment (offline, no engine change).** With the stored P14A per-iteration histories (dr2 proxies, thresholds, steps), compute for si8 and k444s the Hartree-metric estimate `dr2` instead of the plain L2 density residual and the schedule `min(0.1*dr2/nelec, 1e-2)` with the first-iteration redo rule; compare the implied tolerance sequence with the actual one and with QE's (qe_iteration_data_relayed.md). Only if the implied sequence is looser in steps 2-8 but not in the last 3, implement behind a switch.

### B05 Faster convergence of the slow tail root
- **Expected gain.** Basis: 55% of si8 steps occur in five SCF iterations with one late root each (p14a early/late 19/1) and few vectors per step; QE shows the same spikes (12, 11, 20 steps) but only after thresholds tighten and stays at 1-2 otherwise. A method that shortens the tail (adding the previous correction direction to the subspace, or two correction vectors for roots that stall) could reduce step count in those iterations. No estimate beyond the ceiling in B04.
- **Risk.** More vectors per step grows basis and capacity use (A13: larger basis reduces applications but costs memory and Gram work); near-degenerate clusters (A17) are the known hard case.
- **Prerequisites.** B01 data on which roots are slow (buffer roots? top strict root?).
- **Validate.** Gates; skewed and elongated fixtures; degenerate-cluster synthetic from P7PC.
- **Builds on / avoids.** Builds on A17, A53. Avoids A13 (smaller subspace), A15.
- **First experiment.** From the P14A collector report, which band indices are the late roots on si8 and k444s (strict vs buffer) and how many steps each late root needs after all others are below threshold.

### B06 Serial remainder at Gamma: potential, energy audit, mixing, EVD
- **Expected gain.** Measured Mac si8 phases (P10I final): potential 0.019, mixing 0.014 and total-energy audit 0.059 s stay constant for T=1..8 (serial at Gamma), about 0.092 s of 1.77 s at T=8 (5%, derived); EVD 0.139 s of 1.41 s (10%) at T=8 and is sequential by design (OL P10I, P10L). dev08 serial part is 9.2% (0.382 s) of 4.151 s at T=8.
- **Risk.** Reductions must keep fixed association (bit-identical across T): use fixed chunking independent of T (as in A35 and A40). EVD: a different EVD library changes rounding; ADR 0019 reported EVD about equal between MKL and faer (ticket P10H wanted zheev through MKL, never done).
- **Prerequisites.** None for audit/potential; ticket for EVD routing.
- **Validate.** Bit-identical energy and density for the full-grid parallelisation (fixed chunk tree); EVD change needs the energy gate and the 30-pair gate.
- **Builds on / avoids.** Builds on A10 (per-k phases already parallel), A40. Avoids parallel reductions with thread-dependent order (A09 rejected list, A10).
- **First experiment (ready).** On dev08 si8 T=8, profile the 0.382 s serial part (which of audit, potential, mixing, setup, EVD, copies), then parallelise the largest with a fixed-chunk tree and run the 30-pair gate.

### B07 `NonlocalOperator::add` as GEMM
- **Expected gain.** ESTIMATE: nonlocal add is 4.8% of the post-O1 profile, 93 of about 1550 samples (6%) in the Mac pre-P10I profile, 0.474 s (7.5%; 17 of 227 samples) on dev08 MKL T=1 at fe0492b. Best case 5-7% of wall, only if the GEMM form is faster than the current loop at 8-20 atoms.
- **Risk.** Changes atom/channel association order, so rounding and iteration count; the blocked traversal already failed on dev08 (A43).
- **Prerequisites.** Rounding policy; dev08 caller-resolved profile of the function with at least 30 samples.
- **Validate.** Gates; isolated dev08 microbenchmark with 4/8/20 atoms and 1/3/9 channels (shapes used in A43).
- **Builds on / avoids.** Avoids A43 (blocked traversal, measured 2.0x worse on dev08), A41. Builds on A41 and the existing independent traversal oracle test.
- **First experiment.** dev08 profile with at least 30 samples of `add`, split by atom and channel loop; compare against a GEMM of the same shape in the P10J microbenchmark.

### B08 FFT planning and layout
- **Expected gain.** Not estimated: FFT was about 18% of the Mac pre-P10I profile (stale); dev08 FFT share after P14B is not recorded. QE-side measurements (N2 items 2 and 5: stick layout, per-band 2.3 ms down to 0.31 ms on their hardware at 72^3) are not transferable; our P10C gained only 6% on the apply path and 1.75x on density (A26). MKL batched stick transforms are not used (A30).
- **Risk.** Plane-fused z/y/x ordering changes bits (A27); MKL descriptors are plan-owned.
- **Prerequisites.** dev08 profile share of FFT for si8 and k444s.
- **Validate.** Gates; per-call timing 25^3, 36^3, 15x20x24; slab fixture (A26 got slower).
- **Builds on / avoids.** Avoids A11 (lane batching), A28 (band grouping), A34 (FFT-lane threads). Builds on A26, A29, A30, N2 S-A to S-E.
- **First experiment.** dev08 `perf` profile of si8 and k444s (pure and MKL) to get the FFT share after P14B; only if above 20%, benchmark MKL batched stick transforms in the existing fused plan.

### B09 k-point parallel efficiency
- **Expected gain.** Basis: k444s T=1 to T=16 gives 7.1x pure and 6.4x MKL on dev08 (fe0492b); k222s and k222 stop scaling at T=4 to T=8. First-step Davidson iterations differ per k point (31-49 on k444s), so per-k cost varies. Estimate: up to about 10-25% at T=16 from dynamic ordering of k solves by expected cost; not measured.
- **Risk.** Low for ordering (reductions already fold in fixed k order, A09/A10).
- **Prerequisites.** Per-k solve time distribution on dev08.
- **Validate.** Gates; density identical across T.
- **Builds on / avoids.** Builds on A09, A10, A40. Avoids completion-order reductions.
- **First experiment.** dev08 k444s T=16: per-k Davidson time and idle time at the end of the k-solve region (use the `worker-analysis` build); compute the achievable gain from largest-first ordering.

### B10 Cheaper first SCF iteration on many-k (starting vectors)
- **Expected gain.** Basis: first SCF iteration is 4.3% of si8 H vectors (104 of 2416) and 12.6% of k444s (1410 of 11181), derived from p14a work-tables. QE starts from atomic plus random wavefunctions; we start from kinetic coordinate vectors (A07). Estimate: if atomic starts halve that iteration, up to 6% on k444s, 2% on si8. Depends on pseudo-atomic wavefunctions being present in UPF files (real files may lack them).
- **Risk.** Low; changes only the first iteration's trajectory.
- **Prerequisites.** UPF PSWFC parsing or a synthetic atomic guess.
- **Validate.** Gates; k444s and Gamma.
- **Builds on / avoids.** Builds on A07. Avoids nothing rejected.
- **First experiment.** Count first-iteration steps with the kinetic start versus a start from a converged neighbouring k point's eigenvectors (offline in the accounting harness).

### B11 RMM-DIIS at 128 or more bands in a real SCF
- **Expected gain.** None measured. A22: no real SCF crossover to 64 bands (hybrid 3.7x slower at 64 bands); fixed-H pure RMM warm start at 128 bands converged in two sweeps.
- **Label: blocked.** Needs a real or realistic SCF case with at least 128 bands (benchmark suite has `syn-100-bands`) and a measured Davidson share that justifies it.
- **First experiment.** Run `syn-100-bands` with the existing opt-in and compare H columns and wall.

### B12 MPI and GPU
- **Label: blocked.** Out of scope (AGENTS.md: CUDA, Metal, MPI are out unless a ticket says so); GPU design note exists (A32). First step would be an ADR on data distribution; nothing to measure now.

## Not on the backlog (closed, see ledger)
Capacity changes (A08, A13), fuller preconditioner (A14), Ritz-change test with refinement (A15), tolerance multipliers (A21), fewer or more extra bands as default (A19), FFT lane batching/threads (A11, A34), scoped-thread schemes (A33, A35), CholeskyQR2 (A52), finer GEMM tiles (A49), blocked Nonlocal traversal (A43), hybrid RMM as a default (A22). Re-open only with new evidence stated as in ledger section A.
