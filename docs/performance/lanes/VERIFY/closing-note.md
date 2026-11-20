# LANE-VERIFY Phase 2 closing note

Status: **CLOSED — no qualifying new lever established**, 2026-11-20.
This uses the task's explicit closing-note alternative. No optimization is
implemented and no switch, acceptance condition, attempt trigger, tolerance,
continuation or default is changed.

The six relayed dev08 Si64 MPI8 runs show a median verification cost of
1.0475 s out of 8.5466 s external wall. The first attempt fails and the second
passes in every run. Source inspection establishes that default HALF already
feeds failed verification orbitals into the next ordinary map; the MPI kept
start can update H X for the changed potential rather than reapply the full
operator. A failed attempt is therefore not simply discarded. Cache eligibility
is conditional; the relay does not establish how often the kept start was taken.

The next potential defines a different Hamiltonian. Retained orbitals/H X do
not retain the full Davidson search space, and the eigenvalue-change convergence
test restarts. Five iterations and 388 H columns in the second attempt do not
prove redundant work: widths, corrections and operator calls vary. Exact
attribution requires per-iteration evidence absent from this relay.

The reported density, energy and potential components sum to 0.1676 s, below
the requested 0.2 s gain even before identifying necessary fresh work. These
are independently summarized phase maxima, not a coherent wall-time bound.
Energy alone is 0.0499 s. A qualifying solve-only gain would require about
22.23% of the reported solve cost. Full search-space recycling might address
that budget, but no supplied evidence supports that expected gain, and it
requires solver ownership and validity checks beyond the owned verification
paths. Relaxing the solve or delaying acceptance is excluded by the task.

See [diagnosis.md](diagnosis.md) for the contract, product lifetime analysis,
candidate comparison, source provenance and aggregation caveats. The historical
VERIFY-REUSE request predates adoption of HALF; it is not a new opportunity.

## Validation and handoff

This phase changes documentation only. Runtime source and scientific results
are unchanged by construction; no new numerical equivalence result or speedup
is claimed. Mac checks are documentation/integration checks only, not profiling.
Executed on the Mac:

- `CHECK_BASE=HEAD CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast`: exit **1**. The ledger dependency selected dft-solvers and its direct dependents (dft-engine, dft-physics, dft-cli). Build, fmt and clippy passed; 66 test binaries reported **617 passed, 3 failed, 157 ignored**. The three CLI cache failures are the known in-repository TMPDIR rejection (`pseudopotential cache must be outside repositories`). The sandbox also denied the informational `ps` load query; the gates continued.
- `git diff --check`: exit **0**.

No code was changed to bypass these environment failures. These are fresh
Phase 2 check results; Phase 0 numerical identity evidence is historical.

No A/B request is issued because there is no candidate implementation. Phase 0
instrumentation remains available. Reopening needs evidence of avoidable solver
work (including effective switches, kept-start taken/declined counts and
correction widths), or a separately scoped solver proposal. The queued one-core
profile must be assessed in its own parallel model; no one-core equivalent is
extrapolated from MPI. This run stops with changes uncommitted.
