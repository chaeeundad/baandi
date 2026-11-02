# QE reference environment (dev08)

Source: readiness report from the `dft5090` session, 2026-11-02. **Reported by that session, not verified by the rust-dft session.** Re-verify with a real run before relying on it in docs/validation.md.

| Item | Reported value |
| --- | --- |
| Access | `dft5090` session reaches dev08 over key-based ssh (address held by that session; not stored here) |
| pw.x to use | `~/fq/builds/orig-mpi-dfti/bin/pw.x`, PWSCF v.7.6, built 2026-10-19, OpenMPI 5.0.8, MKL 2026.1 (gnu_thread + libgomp linked). Use only this build as the reference; other builds (`r3-*`, `opt-*`, `*-omp-*`) are not references |
| Environment | `source ~/fq/env.sh` first (gcc 13.3, MKL 2026.1, OpenMPI 5.0.8); runs use `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` |
| Runner | `~/fq/run1.sh <orig|r3> mpi <nk> <input> <outdir>`, 4 h timeout, logs wall time and maxrss to `~/fq/runs.log` |
| Hardware | 2x Xeon Gold 6230 (40 physical cores), 187 GB RAM, no scheduler; check `uptime` and top processes before each run and do not run if the node is busy |
| Disk | about 16 GB free in home; prefer `disk_io='low'` for cutoff scans |
| File exchange | scp through the `dft5090` session |
| Pseudopotential | `~/fq/pseudo/Si.pz-vbc.UPF`, sha256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`; NC, `relativistic=no`, LDA PZ, `z_valence=4`. **Author, date and license are blank in the file: provenance unverified** |

## Protocol for a reference run (P6)
1. rust-dft sends: input file(s), their SHA-256, exact command, expected UPF SHA-256.
2. `dft5090` copies to dev08, runs, and returns unmodified: pw.x output, the exact command, QE version, input and UPF SHA-256, wall time, maxrss.
3. No numbers are estimated or rewritten. Results are stored under `tests/reference/` with this provenance.
4. Do not run anything for rust-dft before step 1 is satisfied.
