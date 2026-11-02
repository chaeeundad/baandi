# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in qe_si30_ours.in > qe_si30_ours.out 2> qe_si30_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:31:59
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:31:59Z
- input_sha256: 017e91ea6dc1c1d0e382f3a8158a52248f155379ae360ceb0ef6f3d3f2554065
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.07
- returned_by: dft5090 session (Claude Remote Control), relayed by the 'master' bridge; executed on dev08
- notes: wall_time_s is QE's own 'PWSCF : 0.05s CPU 0.07s WALL' line; the runner's shell estimate was about 0.9 s with 1 s resolution. stdout returned verbatim as gzip+base64 (1 chunk) and verified on receipt by chunk, gzip and stdout SHA-256. stderr (759 bytes: PMIx compression-library warning and an IEEE_DIVIDE_BY_ZERO signalling notice) was not returned. Pseudopotential Si.pz-vbc.UPF: author, date and license unknown, provenance unverified; used with the user's approval for internal validation. Pseudo directory /home/mino/fq/pseudo; input executed unmodified.
- is_format_fixture: False
- uptime_before:  23:31:15 up 268 days,  4:29,  1 user,  load average: 0.00, 0.00, 0.00
- mpi_launcher_options: --map-by core --bind-to core
- pw_x_path: ~/fq/builds/orig-mpi-dfti/bin/pw.x (symlink to ../PW/src/pw.x, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1)

QE functional: SLA  PZ   NOGX NOGC
FFT grids: ours [25, 25, 25]; QE [25, 25, 25]
Tolerances used: {'energy': 0.0001, 'term': 0.0001, 'eigenvalue': 0.001, 'gap': 0.001, 'electrons': 1e-08}
Nondefault tolerances: {}
Informational decomposition reason: None

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.6477652341054636 | -3.6477665149999998 | 1.28089454e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7917834812316733 | 2.7917809249999999 | 2.55623167e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83604605526868414 | 0.836046065 | 9.73131586e-09 | Ha | 0.0001 | PASS |
| xc | -2.5238881382016523 | -2.523888135 | 3.20165228e-09 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718850000005 | 1.84903683e-08 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19011399644747087 | -0.19011526841131599 | 1.27196385e-06 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25863418733601878 | 0.2586343796078247 | 1.92271806e-07 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25863418733659288 | 0.2586343796078247 | 1.92271232e-07 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.2586341873370463 | 0.2586343796078247 | 1.92270778e-07 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33714642191919358 | 0.33714563150389404 | 7.904153e-07 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.3371464219196072 | 0.33714563150389404 | 7.90415713e-07 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33714642192003197 | 0.33714563150389404 | 7.90416138e-07 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37844791825995638 | 0.37844819469711266 | 2.76437156e-07 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078512234582147289 | 0.078511251896069334 | 9.82686078e-07 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 6.6613381477509392e-16 | 6.66133815e-16 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 4.3152503285188984e-09 | 4.31525033e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: 0e832112fd9e03bbf78c11ebbf5c95c7e6a599e466ae18076eb73e92b8d9cb42
QE output SHA-256: 9ee1a3a2feb9f5ddc26e5271ca1084a3b732fde363cbeaeb2ebefe308add2d41
Maximum differences:
- eigenvalues: eigenvalue band 1, 1.27196385e-06 Ha
- energy terms: one_electron, 2.55623167e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
