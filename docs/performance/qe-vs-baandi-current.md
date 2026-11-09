# QE 7.6 vs BAANDI, consolidated speed table at 8114447 (2026-11-09)

**Parallel models differ: the QE multi-core columns are MPI (`mpirun -np N`, OMP and MKL threads 1) and the BAANDI multi-core columns are shared-memory threads (`threads = N`, no MPI layer exists yet). Only the 1-core cells compare like with like; the 4- and 8-core ratios mix two parallel models and include the MPI-versus-threads difference (QE 4.26x on 8 ranks against BAANDI 2.55x on 8 threads at si32 strict). The like-for-like multi-core comparison needs BAANDI's MPI layer (ADR 0041, M1).**

Source: dev08 master session (Xeon Gold 6230), relayed 2026-11-09; raw data /scr/mino/rustdft/{qenp_out.txt, final8_out.txt, probe_out.txt, stageA/}. BAANDI at 8114447: all current defaults (gamma_real auto, MKL dgemm R1/R2, blocked MGS on MKL, FIN1, AUD1, dsyevd on MKL), no gamma_real key, real Gamma path used in every run, n = 5/5/3, all converged, density SHA stable. Builds: cli-mkl d7a5a5279c1bead8, cli-pure d0a80be02ec08e85. QE 7.6 (orig-mpi-dfti), stage-A inputs and the same FFT grids, `K_POINTS gamma`, nosym/noinv, `mpirun -np N`, no -nk/-nd, OMP/MKL threads 1. Mapping: QE conv_thr 1e-8 <-> BAANDI standard profile, 1e-10 <-> strict (certified); this is a convention, not an equivalence (see the offset below). Walls only: QE work (H vectors) was not measured, no work ratio is quoted. BAANDI's `operator_applications` counts columns (band vectors).

## QE vs BAANDI per (case, class, cores), seconds
"QE wall" is the measured mpirun wall; "QE PWSCF" is QE's own timer (the launcher adds up to 0.15 s, which matters only at si8).
| case | class | cores | QE wall | QE PWSCF | BAANDI MKL | BAANDI pure | better | QE faster (wall) | (PWSCF) |
|---|---|---|---|---|---|---|---|---|---|
| si8 | standard | 1 | 0.469 | 0.36 | 1.005 | 1.054 | MKL | 2.1x | 2.8x |
| si8 | standard | 4 | 0.262 | 0.15 | 0.792 | 0.818 | MKL | 3.0x | 5.3x |
| si8 | standard | 8 | 0.295 | 0.14 | 0.803 | 0.838 | MKL | 2.7x | 5.7x |
| si8 | strict | 1 | 0.552 | 0.45 | 1.827 | 2.089 | MKL | 3.3x | 4.1x |
| si8 | strict | 4 | 0.282 | 0.17 | 1.343 | 1.475 | MKL | 4.8x | 7.9x |
| si8 | strict | 8 | 0.284 | 0.17 | 1.409 | 1.502 | MKL | 5.0x | 8.3x |
| si16 | standard | 1 | 1.413 | 1.31 | 3.139 | 2.939 | pure | 2.1x | 2.2x |
| si16 | standard | 4 | 0.609 | 0.50 | 1.955 | 1.794 | pure | 2.9x | 3.6x |
| si16 | standard | 8 | 0.565 | 0.44 | 1.870 | 1.684 | pure | 3.0x | 3.8x |
| si16 | strict | 1 | 1.831 | 1.72 | 7.523 | 8.871 | MKL | 4.1x | 4.4x |
| si16 | strict | 4 | 0.703 | 0.60 | 4.362 | 4.951 | MKL | 6.2x | 7.3x |
| si16 | strict | 8 | 0.688 | 0.53 | 4.103 | 4.602 | MKL | 6.0x | 7.7x |
| si32 | standard | 1 | 5.028 | 4.94 | 14.037 | 16.002 | MKL | 2.8x | 2.8x |
| si32 | standard | 4 | 1.961 | 1.89 | 6.776 | 8.142 | MKL | 3.5x | 3.6x |
| si32 | standard | 8 | 1.700 | 1.59 | 5.904 | 7.093 | MKL | 3.5x | 3.7x |
| si32 | strict | 1 | 6.700 | 6.59 | 40.753 | 46.477 | MKL | 6.1x | 6.2x |
| si32 | strict | 4 | 2.440 | 2.35 | 18.677 | 21.184 | MKL | 7.7x | 7.9x |
| si32 | strict | 8 | 1.984 | 1.91 | 16.001 | 17.501 | MKL | 8.1x | 8.4x |
QE is faster in all 18 cells. MKL is the better BAANDI build in 15 of 18 cells; the three si16 standard cells favour the pure build by 6 to 10 percent (a trajectory effect: the pure runs take fewer Davidson steps there).

## Thread scaling (best BAANDI build; QE wall)
| case | BAANDI T=1 / 4 / 8 | speedup | QE 1 / 4 / 8 | speedup |
|---|---|---|---|---|
| si8 standard | 1.005 / 0.792 / 0.803 | 1.25x | 0.469 / 0.262 / 0.295 | 1.22x |
| si8 strict | 1.827 / 1.343 / 1.409 | 1.30x | 0.552 / 0.282 / 0.284 | 1.58x |
| si16 standard | 2.939 / 1.794 / 1.684 | 1.75x | 1.413 / 0.609 / 0.565 | 2.32x |
| si16 strict | 7.523 / 4.362 / 4.103 | 1.83x | 1.831 / 0.703 / 0.688 | 2.50x |
| si32 standard | 14.037 / 6.776 / 5.904 | 2.38x | 5.028 / 1.961 / 1.700 | 2.91x |
| si32 strict | 40.753 / 18.677 / 16.001 | 2.55x | 6.700 / 2.440 / 1.984 | 3.32x |
Neither code scales linearly at 8 cores; QE's scaling advantage is modest. The gap is already present at one core (si32 strict 6.1x), so single-core cost, not scaling, is the main factor. si32 strict T=1 fell from 50.92 s to 40.753 s with the dsyevd default (ratio 7.7x to 6.1x).

## Energy offset between the codes (unequal targets)
| case | QE 1e-10 (Ha) | BAANDI default (Ha) | difference | per atom | with pseudo.local_mesh=qe_like and radial_even_rule=simpson_quadratic_end |
|---|---|---|---|---|---|
| si8 | -31.321813785 | -31.321803601 | +1.018e-5 | +1.273e-6 | unchanged (1e-14 Ha) |
| si16 | -63.129264785 | -63.129244422 | +2.036e-5 | +1.273e-6 | unchanged |
| si32 | -126.704038825 | -126.703998096 | +4.073e-5 | +1.273e-6 | unchanged |
The offset is +1.273e-6 Ha per atom to four figures in every case (perfectly size-extensive) and is not touched by the radial-convention keys, which closed the PseudoDojo PBE cases. It points to a per-atom constant in the atomic or pseudopotential treatment of Si.pz-vbc.UPF (not a grid, convergence or SCF-target effect). The QE-vs-BAANDI ratios above therefore compare two runs converged to numbers 4.1e-5 Ha apart (si32); the tolerance mapping is a convention. An isolated-atom probe would isolate the term (open).

## Reading
1. QE is 2.1x to 3.5x faster in the standard class and 3.3x to 8.1x in the strict class (wall); the strict class includes BAANDI's residual certificate, which QE does not do.
2. The ratio grows with the system size and with the number of cores in the strict class (si32 strict 6.1x at 1 core, 8.1x at 8 cores).
3. Unperturbed single starts; a perturbed-start table for si32 (seeds 501 to 505, both codes) is pending.

## Update 2026-11-09 (later): QE best start (startingwfc = atomic) is the headline
dev08 found that QE's default start (atomic+random wavefunctions, free-atom potential) is not its best: `startingwfc = 'atomic'` is faster in all 18 cells and freezes the SCF count (si32 1e-10: 8 SCF and 35 h_psi calls at 1, 4 and 8 cores; default start 13 SCF, 186 h_psi). The project target is the comparison against QE's best configuration, so the ratios against QE best are the headline; the default-start ratios above are the out-of-the-box line. BAANDI columns unchanged (8114447).
| case | class | cores | QE default | QE atomic (best) | BAANDI best | ratio vs QE best | ratio vs QE default |
|---|---|---|---|---|---|---|---|
| si8 | standard | 1 / 4 / 8 | 0.469 / 0.262 / 0.295 | 0.456 / 0.245 / 0.233 | 1.005 / 0.792 / 0.803 | 2.2x / 3.2x / 3.4x | 2.1x / 3.0x / 2.7x |
| si8 | strict | 1 / 4 / 8 | 0.552 / 0.282 / 0.284 | 0.533 / 0.281 / 0.265 | 1.827 / 1.343 / 1.409 | 3.4x / 4.8x / 5.3x | 3.3x / 4.8x / 5.0x |
| si16 | standard | 1 / 4 / 8 | 1.413 / 0.609 / 0.565 | 1.080 / 0.451 / 0.402 | 2.939 / 1.794 / 1.684 | 2.7x / 4.0x / 4.2x | 2.1x / 2.9x / 3.0x |
| si16 | strict | 1 / 4 / 8 | 1.831 / 0.703 / 0.688 | 1.330 / 0.535 / 0.463 | 7.523 / 4.362 / 4.103 | 5.7x / 8.2x / 8.9x | 4.1x / 6.2x / 6.0x |
| si32 | standard | 1 / 4 / 8 | 5.028 / 1.961 / 1.700 | 3.000 / 0.953 / 0.722 | 14.037 / 6.776 / 5.904 | 4.7x / 7.1x / 8.2x | 2.8x / 3.5x / 3.5x |
| si32 | strict | 1 / 4 / 8 | 6.700 / 2.440 / 1.984 | 3.518 / 1.106 / 0.827 | 40.753 / 18.677 / 16.001 | 11.6x / 16.9x / 19.4x | 6.1x / 7.7x / 8.1x |
Against QE best: minimum 2.2x, median 5.0x, maximum 19.4x. QE's atomic start also improves its thread scaling (si32 1e-10, 1 to 8 cores 4.26x, identical work at every core count) against BAANDI's best of 2.55x. QE's SCF count does not survive hard geometries (si-skew 51 / 180 iterations, si-elongated 137 and a non-converged 1e-10 run, si-slab non-converged at 200 iterations, fixed occupations, ecut 8 Ry robustness cells); BAANDI's numbers on those cells are being measured.
