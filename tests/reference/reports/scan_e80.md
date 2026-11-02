# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e80_ours.in > e80_ours.out 2> e80_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:30
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:30Z
- input_sha256: f8aa10f3e2ad81799bdf4b211923fb765704570c582a5d8ebd266e8228e92d1b
- upf_sha256: e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217
- wall_time_s: 0.18
- returned_by: dft5090 session (Claude Remote Control), relayed by the 'master' bridge; executed on dev08
- notes: Cutoff-scan run. Input produced by editing prefix, outdir, ecutwfc, ecutrho and nr1-3 of the verified 30 Ry input; size and sha256 verified by the runner before running. wall_time_s is QE's own WALL value. stdout returned inside one tar+gzip+base64 stream of 6 outputs (8460 characters), verified by piece, whole stream, gzip, tar and per-file SHA-256 on receipt (several pieces had to be resent because of transcription errors on the receiving side). stderr not returned. Pseudopotential Si.pz-vbc.UPF: provenance unverified; used with the user's approval for internal validation.
- is_format_fixture: False
- uptime_before_load_average: 0.10, 0.03, 0.01
- mpi_launcher_options: --map-by core --bind-to core
- pw_x_path: ~/fq/builds/orig-mpi-dfti/bin/pw.x (symlink to ../PW/src/pw.x, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1)

QE functional: SLA  PZ   NOGX NOGC
FFT grids: ours [45, 45, 45]; QE [45, 45, 45]
Tolerances used: {'energy': 0.0001, 'term': 0.0001, 'eigenvalue': 0.001, 'gap': 0.001, 'electrons': 1e-08}
Nondefault tolerances: {}
Informational decomposition reason: None

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.648020455029771 | -3.6480217275000002 | 1.27247023e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7911239167958297 | 2.791121365 | 2.55179583e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83627998053457875 | 0.83627996999999998 | 1.05345788e-08 | Ha | 0.0001 | PASS |
| xc | -2.5239729408803186 | -2.5239729249999998 | 1.58803188e-08 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718650000006 | 1.50963153e-09 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19016319263843559 | -0.19016304253014432 | 1.50108291e-07 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25860005342295639 | 0.25860130521786662 | 1.25179491e-06 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25860005342342041 | 0.25860130521786662 | 1.25179445e-06 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25860005342580666 | 0.25860130521786662 | 1.25179206e-06 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33709629264978402 | 0.33709418245284811 | 2.11019694e-06 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.33709629265181085 | 0.33709418245284811 | 2.11019896e-06 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33709629265216601 | 0.33709418245284811 | 2.11019932e-06 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.3777521427933655 | 0.37775363250799276 | 1.48971463e-06 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.07849623922397736 | 0.078492877234981495 | 3.361989e-06 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 1.2967404927621828e-13 | 1.29674049e-13 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 8.5837679083237925e-09 | 8.58376791e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: 02b093a30c4684e8dc8735a9d7a6218ce7125c2c4f70797f5708fb0ba17d2923
QE output SHA-256: dc45e4dee619146e18a27f606528118d226f92ddea23b7a3679861c0cee88064
Maximum differences:
- eigenvalues: eigenvalue band 7, 2.11019932e-06 Ha
- energy terms: one_electron, 2.55179583e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
