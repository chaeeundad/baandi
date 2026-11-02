# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e40_ours.in > e40_ours.out 2> e40_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:29
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:29Z
- input_sha256: ec901eb4ebc1e3d459a04a00fbf11d9e9ee6341ba072c03d15c2d5abdce3365b
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.08
- returned_by: dft5090 session (Claude Remote Control), relayed by the 'master' bridge; executed on dev08
- notes: Cutoff-scan run. Input produced by editing prefix, outdir, ecutwfc, ecutrho and nr1-3 of the verified 30 Ry input; size and sha256 verified by the runner before running. wall_time_s is QE's own WALL value. stdout returned inside one tar+gzip+base64 stream of 6 outputs (8460 characters), verified by piece, whole stream, gzip, tar and per-file SHA-256 on receipt (several pieces had to be resent because of transcription errors on the receiving side). stderr not returned. Pseudopotential Si.pz-vbc.UPF: provenance unverified; used with the user's approval for internal validation.
- is_format_fixture: False
- uptime_before_load_average: 0.10, 0.03, 0.01
- mpi_launcher_options: --map-by core --bind-to core
- pw_x_path: ~/fq/builds/orig-mpi-dfti/bin/pw.x (symlink to ../PW/src/pw.x, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1)

QE functional: SLA  PZ   NOGX NOGC
FFT grids: ours [30, 30, 30]; QE [30, 30, 30]
Tolerances used: {'energy': 0.0001, 'term': 0.0001, 'eigenvalue': 0.001, 'gap': 0.001, 'electrons': 1e-08}
Nondefault tolerances: {}
Informational decomposition reason: None

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.647964754770058 | -3.6479660350000001 | 1.28022994e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7912795360615954 | 2.791276995 | 2.5410616e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83622083813818893 | 0.83622081000000004 | 2.81381889e-08 | Ha | 0.0001 | PASS |
| xc | -2.5239580172302687 | -2.523957995 | 2.22302687e-08 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718850000005 | 1.84903683e-08 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19014295490211169 | -0.19014466786905648 | 1.71296694e-06 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25860279569069516 | 0.25860130521786662 | 1.49047283e-06 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.2586027959218764 | 0.25860130521786662 | 1.49070401e-06 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25860279610514175 | 0.25860130521786662 | 1.49088728e-06 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33710602391219191 | 0.33710520724950077 | 8.16662691e-07 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.33710602408684837 | 0.33710520724950077 | 8.16837348e-07 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33710602417709834 | 0.33710520724950077 | 8.16927598e-07 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37802592024446907 | 0.37802557749209265 | 3.42752376e-07 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078503227807050158 | 0.078503902031634154 | 6.74224584e-07 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 7.7715611723760958e-16 | 7.77156117e-16 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 8.0353204315056767e-09 | 8.03532043e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: 517218c3c3a95b9c2f8b3c36b5043ea08a002da6ce09ec9ec05beed13a8984d7
QE output SHA-256: ba2de060faedf8ddd2fab404e1020ea5d41e37ddce878dc1df2796b8572624c8
Maximum differences:
- eigenvalues: eigenvalue band 1, 1.71296694e-06 Ha
- energy terms: one_electron, 2.5410616e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
