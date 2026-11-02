# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e20_ours.in > e20_ours.out 2> e20_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:28
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:28Z
- input_sha256: 0ed41ae68e43f483b1c5290adf71a3fa39dc1713ebe38d8606664e0ab39c4a69
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.05
- returned_by: dft5090 session (Claude Remote Control), relayed by the 'master' bridge; executed on dev08
- notes: Cutoff-scan run. Input produced by editing prefix, outdir, ecutwfc, ecutrho and nr1-3 of the verified 30 Ry input; size and sha256 verified by the runner before running. wall_time_s is QE's own WALL value. stdout returned inside one tar+gzip+base64 stream of 6 outputs (8460 characters), verified by piece, whole stream, gzip, tar and per-file SHA-256 on receipt (several pieces had to be resent because of transcription errors on the receiving side). stderr not returned. Pseudopotential Si.pz-vbc.UPF: provenance unverified; used with the user's approval for internal validation.
- is_format_fixture: False
- uptime_before_load_average: 0.02, 0.01, 0.00
- mpi_launcher_options: --map-by core --bind-to core
- pw_x_path: ~/fq/builds/orig-mpi-dfti/bin/pw.x (symlink to ../PW/src/pw.x, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1)

QE functional: SLA  PZ   NOGX NOGC
FFT grids: ours [24, 24, 24]; QE [24, 24, 24]
Tolerances used: {'energy': 0.0001, 'term': 0.0001, 'eigenvalue': 0.001, 'gap': 0.001, 'electrons': 1e-08}
Nondefault tolerances: {}
Informational decomposition reason: None

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.6455446820759461 | -3.6455459600000002 | 1.27792405e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7967108717580555 | 2.7967082900000002 | 2.58175806e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83456395995937305 | 0.83456401000000002 | 5.0040627e-08 | Ha | 0.0001 | PASS |
| xc | -2.522892329359689 | -2.5228923449999998 | 1.56403108e-08 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718749999997 | 8.49036752e-09 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19000082285081901 | -0.19000134551257145 | 5.22661752e-07 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25914178438830149 | 0.25914152025384873 | 2.64134453e-07 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.259141784429953 | 0.25914152025384873 | 2.64176104e-07 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25914178443390651 | 0.25914152025384873 | 2.64180058e-07 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33752354115915856 | 0.33752414952230325 | 6.08363145e-07 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.33752354116209227 | 0.33752414952230325 | 6.08360211e-07 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33752354118934103 | 0.33752414952230325 | 6.08332962e-07 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37840457453866416 | 0.37840409551050191 | 4.79028162e-07 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.07838175672525205 | 0.078382629268454518 | 8.72543202e-07 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 5.5511151231257827e-16 | 5.55111512e-16 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 8.987514296722391e-09 | 8.9875143e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: ceac65dc56ec13e612426da1d58f6a387eed9938f9a2db2764db1a3a8bb4658f
QE output SHA-256: 4bef1213148031c0b0b3224bfd2bbbe2ef063fa6346a619699a59c5921ae875b
Maximum differences:
- eigenvalues: eigenvalue band 5, 6.08363145e-07 Ha
- energy terms: one_electron, 2.58175806e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
