# QE comparison
Pseudopotential provenance is unverified

PASS

## Provenance

- command: cd ~/rustdft_p6
source ~/fq/env.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 14400 mpirun -np 4 --map-by core --bind-to core ~/fq/builds/orig-mpi-dfti/bin/pw.x -in e25_ours.in > e25_ours.out 2> e25_ours.err
- qe_version: Program PWSCF v.7.6 starts on  2Nov2026 at 23:45:29
- host: dev08 (192.168.0.228), Rocky Linux 9.7
- mpi_processes: 4
- omp_threads: 1
- date_utc: 2026-11-02T14:45:29Z
- input_sha256: 30efd70d170c3f49a5ed6d6248a21c5ac61e053714ed1c343313bd4456c88c48
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
| total energy per atom | -3.6471340591845864 | -3.6471353375 | 1.27831541e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7930807125344295 | 2.793078145 | 2.56753443e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83562511479689772 | 0.83562513999999999 | 2.52031023e-08 | Ha | 0.0001 | PASS |
| xc | -2.5235020791908673 | -2.5235020800000001 | 8.09132761e-10 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718749999997 | 8.49036752e-09 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19013779944881418 | -0.19013731800462136 | 4.81444193e-07 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25880255956765891 | 0.25880342648983268 | 8.66922174e-07 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25880255956774079 | 0.25880342648983268 | 8.66922092e-07 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.25880255956780746 | 0.25880342648983268 | 8.66922025e-07 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33733102654181535 | 0.33732937811477232 | 1.64842704e-06 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.33733102654191149 | 0.33732937811477232 | 1.64842714e-06 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33733102654205677 | 0.33732937811477232 | 1.64842728e-06 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37856953978227759 | 0.3785694674602923 | 7.23219853e-08 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078528466974007893 | 0.078525951624939638 | 2.51534907e-06 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 2.3314683517128287e-15 | 2.33146835e-15 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 4.1442355434818105e-09 | 4.14423554e-09 | Ha | 1e-07 | PASS |

Our results SHA-256: a6ffee4178ee8641b559f6a85cb5e953a480e70746583b1e872c2fffb9d1237b
QE output SHA-256: 1559e78f73b8fd6305053fd90626326c428cf50ae278d24ec9afdb2a60d57762
Maximum differences:
- eigenvalues: eigenvalue band 7, 1.64842728e-06 Ha
- energy terms: one_electron, 2.56753443e-06 Ha

Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.
QE electron count is configured count, not an independent density integral.
1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.
