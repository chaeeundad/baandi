# dev08 request: MPI stage M2d (batched transposes and distributed SCF stages), 2026-11-10

Branch: the M2d worktree branch named in the ledger row MPI-M2D (commit in the hand-off). Build exactly like the M1 gate binaries: `--features mpi`, MKL (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`), Open MPI 5.0.8, `MKL_NUM_THREADS=1`, one binary for every configuration. Adoption rule (AGENTS.md, user 2026-11-10): any measurable improvement of the median wall over M1 on dev08 is adopted, even a few percent; the change is rejected only for physical hard-coding or a correctness failure (accuracy classes strict 1e-10 Ha, standard 5e-9 Ha, certificate, convergence, determinism). The percentage gates of the M2e plan are superseded by that rule for adoption; they stay as reported numbers.

## What changed (all only with `BAANDI_MPI_M2=on` and 2+ ranks; default and np=1 paths unchanged)
- Batched FFT transposes: one alltoall per direction for all column pairs of an H application (was one per pair). Alltoall calls per apply of m columns: 2 + 2 mirror exchanges (was 2*ceil(m/2) + 2); results bit-identical. `BAANDI_M2_BATCH=off` restores the per-pair transposes for an A/B.
- Density: built from the distributed eigenvectors on each rank's FFT planes (stick FFT, one batched transpose), planes assembled with one allgatherv. Per-band plane contributions feed the certificate screen (one tree_sum).
- Potential: LDA (PZ, PW) XC split over ranks by contiguous grid blocks (bit-identical potential; energy and potential integral by tree_sum, one allgatherv). Hartree stays replicated.
- Energies: audit nonlocal energy from the distributed projections (one tree_sum); output-density XC energy as above.
- Not distributed: mixing (0.12 s of 41 s in MPI-M1-GATE4), Hartree FFT, kinetic energy sum, symmetrization.

## Runs
Cell: si32 strict (`benchmarks/cases/si32-gamma`), perturbed seeds 601-610 (+-1e-7 Bohr as in MPI-M1-GATE3). Also si32 standard for seeds 601-603 (accuracy class check). Load-gated (load below 0.30), 8 physical cores of socket 0, interleaved per seed.

| Config | env | ranks x threads |
| --- | --- | --- |
| baseline | `BAANDI_MPI_M2` unset | np=1 x T=1 (same binary) |
| M1 | `BAANDI_MPI_M2` unset | np=4 x T=1, np=8 x T=1 |
| M2d | `BAANDI_MPI_M2=on` | np=4 x T=1, np=8 x T=1, np=2 x T=4, np=4 x T=2 |
| M2d unbatched (A/B) | `BAANDI_MPI_M2=on BAANDI_M2_BATCH=off` | np=8 x T=1, seeds 601-603 only |
| threads-only | `BAANDI_MPI_M2` unset | np=1 x T=8 |

## Proof that the M2 branch ran (every M2d run; reject the row without it)
`results.json` `parallel.model == "mpi_stick_rows"`, `parallel.m2.enabled == true`, and in `parallel.m2.counters`: `rows_solves > 0`, `batched_applies > 0` (`unbatched_applies > 0` instead for the A/B row), `density_builds > 0`, `band_density_builds > 0`, `xc_evaluations > 0`, `nonlocal_energy_audits > 0`; rank directories present; energy differs from np=1 at 1e-13 or below (nonzero is expected) and lies within 1e-10 Ha (strict) / 5e-9 Ha (standard).

## What to send back
Per run: wall, `phases` timers, `parallel.stages` (M1 Davidson stages), `parallel.m2.stages` (engine `density`, `xc`, `nonlocal_energy`, `screen_integrals`; operator `op.apply_mirror_exchange`, `op.apply_pack_unpack_kinetic`, `op.apply_fft_transposes`, `op.apply_nonlocal`, `op.band_densities`, `op.nonlocal_energy`, `op.gather_planes`), `parallel.collectives` (calls and seconds per op), `parallel.m2.counters`, `convergence.work`, energy (17 digits), density SHA-256 of rank 0 and rank 1, exit code. Determinism: one seed twice at np=4 M2d (density SHA and energy bit-identical). Binary sha256, commit, build commands, uptime before and after.

## Decision
Median wall over seeds per configuration. M2d is adopted as the MPI default (`BAANDI_MPI_M2` default on) when its best configuration beats M1 at the same rank count with no correctness failure; the batched transposes stay on unless the A/B row shows them slower. The M2e numbers (np=4 >= 2.2x, np=8 >= 3.5x over np=1, communication below 25 percent of wall at np=8, best MPI configuration against threads-only np=1 x T=8) are reported for the record.
