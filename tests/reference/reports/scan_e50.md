# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e50_ours.in > e50_ours.out 2> e50_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:29
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:29Z
- input_sha256: 1a9e94e7b4750593c609ab2fcd6e05e42f31dfd4c1ce9ea4cbbca91f08bf056a
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.09
- returned_by: dft5090 session (Claude Remote Control), relayed by the 'master' bridge; executed on dev08
- notes: Cutoff-scan run. Input produced by editing prefix, outdir, ecutwfc, ecutrho and nr1-3 of the verified 30 Ry input; size and sha256 verified by the runner before running. wall_time_s is QE's own WALL value. stdout returned inside one tar+gzip+base64 stream of 6 outputs (8460 characters), verified by piece, whole stream, gzip, tar and per-file SHA-256 on receipt (several pieces had to be resent because of transcription errors on the receiving side). stderr not returned. Pseudopotential Si.pz-vbc.UPF: provenance unverified; used with the user's approval for internal validation.
- is_format_fixture: False
- uptime_before_load_average: 0.10, 0.03, 0.01
- mpi_launcher_options: --map-by core --bind-to core
- pw_x_path: ~/fq/builds/orig-mpi-dfti/bin/pw.x (symlink to ../PW/src/pw.x, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1)

QE functional: SLA  PZ   NOGX NOGC
FFT grids: ours [36, 36, 36]; QE [36, 36, 36]
Tolerances used: {'energy': 0.0001, 'term': 0.0001, 'eigenvalue': 0.001, 'gap': 0.001, 'electrons': 1e-08}
Nondefault tolerances: {}
Informational decomposition reason: None

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.6480095273224071 | -3.6480108100000002 | 1.28267759e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7911572275461158 | 2.7911546650000001 | 2.56254612e-06 | Ha | 0.0001 | PASS |
| hartree | 0.8362676954314241 | 0.83626769999999995 | 4.56857585e-09 | Ha | 0.0001 | PASS |
| xc | -2.5239721111127227 | -2.5239720999999999 | 1.11127227e-08 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718850000005 | 1.84903683e-08 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.1901586165423127 | -0.19015936759792676 | 7.51055614e-07 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25860016692219695 | 0.25860130521786662 | 1.13829567e-06 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25860016692247068 | 0.25860130521786662 | 1.1382954e-06 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25860016692289284 | 0.25860130521786662 | 1.13829497e-06 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33709670388941537 | 0.3370978573850657 | 1.15349565e-06 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.33709670388988544 | 0.3370978573850657 | 1.15349518e-06 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33709670389030411 | 0.3370978573850657 | 1.15349476e-06 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37778486507046871 | 0.37778303196573326 | 1.83310474e-06 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078496536966522534 | 0.078496552167199085 | 1.52006766e-08 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 1.1546319456101628e-14 | 1.15463195e-14 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 5.7042593082345806e-09 | 5.70425931e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: aca9f6a8c9b89abc52b844c9fee906a10da036c779367b2190dab775807676c0
QE output SHA-256: 1c54224bd7f444d9b57538d3f5142ed0217c4fe39c31e8221f209ca0be395d3c
Maximum differences:
- eigenvalues: eigenvalue band 8, 1.83310474e-06 Ha
- energy terms: one_electron, 2.56254612e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
