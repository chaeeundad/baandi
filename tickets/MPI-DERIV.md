TASK_ID: MPI-DERIV
GOAL: Under `--features mpi`, the post-SCF forces and stress of the reference path are distributed over ranks by k block (each rank evaluates its own k blocks, one allreduce of the per-k sums), instead of every rank evaluating all k blocks redundantly. Also write the energy/eigenvalue part of results (an early results file or a console energy line) BEFORE derivatives start, so a run stopped during derivatives still yields energy and gap. Approved by the user 2026-11-26 (option A), started after D3-DISPERSION merges.

EVIDENCE: dev08 MPQE entry 17 (B12, 12 atoms, 172 irreducible k, 54^3 grid, 22 bands, MPI8): SCF + verification done in ~9 min, then >36 min in forces/stress until the 45 min cap; results.json never written. QE: 84 s total. dft_engine::forces / stress take no communicator and work on replicated coefficients of all k blocks (forces.rs DerivativeBlock, "complete replicated coefficients").

OWNED_FILES: crates/dft-engine/src/forces.rs, crates/dft-engine/src/stress*.rs, crates/dft-engine/src/refscf/** derivative/runtime pieces that evaluate forces and stress (reference_derivatives), crates/dft-cli/src/run.rs (derivative call sites, early energy record), crates/dft-mpi/** (only if a reduction helper is missing), tests under crates/dft-cli/tests/ and crates/dft-engine/tests/ for the new behaviour, docs/roadmap/mpi-* (one evidence paragraph), docs/performance/attempts-ledger.md (one row).

READ_ONLY_CONTEXT: AGENTS.md, ownership.toml, crate AGENTS.md files, docs/decisions/0041-mpi-parallel-layer.md (determinism contract), ADR 0043/0047, how the SCF already distributes k blocks and reduces under MPI.

CONTRACT:
- Threads-only default build: bit-identical results (RG-ID baselines, `bench_scf.py --require-identical`).
- MPI: forces and stress equal to the replicated evaluation within 1e-12 Ha/Bohr and 1e-12 Ha/Bohr^3 (summation order may change; state the measured max difference). Result independent of rank count within that bound (np 1, 2, 3 on a k-mesh case with uneven k split).
- Symmetrization stays after the reduction (same as now).
- Early record: results contain energies and eigenvalues even if derivatives are interrupted; final results.json unchanged in schema except an optional status field.

TEST_COMMANDS: scripts/check.sh; local mpirun np 1/2/3 on a small k-mesh fixture with forces+stress if an MPI toolchain exists on the Mac, else the orchestrator runs it on dev08.

ACCEPTANCE: dev08 B12 rerun with forces+stress under MPI8 finishes, derivative wall reported; agreement numbers above.

NON_GOALS: G-vector distribution of derivatives inside one k block; speeding up the serial derivative kernels.

GIT: leave changes uncommitted; orchestrator commits.
RETURN_FORMAT: docs/schemas/report.schema.json
