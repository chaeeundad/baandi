# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e60_ours.in > e60_ours.out 2> e60_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:30
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:30Z
- input_sha256: dfc60b0d13a2a9ba48cd77610d58f65c83d23a01ac93e4324c107f5260354074
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.1
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
| total energy per atom | -3.6480193039437112 | -3.6480205775000001 | 1.27355629e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7911273479267633 | 2.79112479 | 2.55792676e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83627869641678709 | 0.83627870000000004 | 3.58321295e-09 | Ha | 0.0001 | PASS |
| xc | -2.5239727857213405 | -2.5239727799999998 | 5.72134073e-09 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718699999993 | 3.49036711e-09 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19016299489924579 | -0.19016304253014432 | 4.76308985e-08 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25860016552089449 | 0.25860130521786662 | 1.13969697e-06 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25860016552223103 | 0.25860130521786662 | 1.13969564e-06 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25860016552845805 | 0.25860130521786662 | 1.13968941e-06 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33709653907369064 | 0.3370978573850657 | 1.31831138e-06 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.337096539074285 | 0.3370978573850657 | 1.31831078e-06 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33709653907529036 | 0.3370978573850657 | 1.31830978e-06 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.3777584347527404 | 0.3777573074402103 | 1.12731253e-06 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078496373545232589 | 0.078496552167199085 | 1.78621966e-07 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 2.3980817331903381e-14 | 2.39808173e-14 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 8.0565597955287423e-09 | 8.0565598e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: b3074a1752bbbb33a2ad716fa0f17681e2fa6038841462498eeb2c7335c86813
QE output SHA-256: 28a260d3835d7d1e91a248cffc92c5c177de8236495900059ca1479d24a64fe0
Maximum differences:
- eigenvalues: eigenvalue band 5, 1.31831138e-06 Ha
- energy terms: one_electron, 2.55792676e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
