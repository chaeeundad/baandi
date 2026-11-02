# MVP review (Claude Code, 2026-11-03)

Independent review of the MVP-CPU scope defined in the project specification: Gamma-only, spin-unpolarized, norm-conserving, LDA plane-wave Kohn-Sham SCF for diamond Si. This is a review record, not an approval: **final MVP approval belongs to the user**.

## Verdict
All spec gates P0 to P6 are met with evidence recorded in [validation.md](validation.md) and the dashboard data in `artifacts/`. Nothing in the evidence contradicts the claim that the engine reproduces QE 7.6 for this system and this pseudopotential within the agreed tolerances. The claim does not extend beyond that scope (see "What is not established").

## Checklist (spec section 13)
| Gate | Status | Evidence |
| --- | --- | --- |
| P0 environment, toolchain, workspace | met | Rust 1.99.1 pinned; scripts/check.sh; the MCP bridge to Codex no longer exists in Codex CLI, so `codex exec` plus scripts/delegate.sh was used (user approved) |
| P1 lattice, G-vectors, FFT mathematics | met | physics-conventions.md; tests with stated tolerances; Si 30 Ry: 749 plane waves, 5985 density G, FFT 25^3 |
| P2 real Si UPF and local/nonlocal terms | met | Real file parsed after three real-world deviations were fixed (P2C, P5C); projector norm identity in an isolated cell to 3e-13; translation phases, Hermitian D |
| P3 Hamiltonian, Hartree, XC | met | Matrix-free H against an independently built matrix (9e-16); XC derivative test differentiates the energy density; Gaussian electrostatic cross-check of the G=0/Ewald bookkeeping (48 cases, 3.3e-14) |
| P4 eigensolver | met | Dense (faer) versus block Davidson, degenerate and near-degenerate cases, explicit non-convergence errors |
| P5 SCF and total energy | met | Real Si converges (11 iterations); convergence needs density, energy, eigen residual and charge criteria plus a fixed-point certification step |
| P6 QE comparison and cutoff convergence | met | R3 to R5: total energy 1.28e-6 Ha/atom (tolerance 1e-4), eigenvalues <= 2.1e-6 Ha (1e-3), same grid, same iteration count; six cutoffs 20 to 80 Ry; density in Fourier space |
| Documentation | met except license | README, architecture, dependencies, conventions, UPF support, CLI, validation, ADR 0001 to 0010 |
| Claude final review | this document | |
| User approval | **pending** | |

## Measured facts used in the verdict
- Rust: 9392 lines (src and tests) in 9 crates; 70 fast tests, 9 slow release tests (all pass), 90 Python tests (all pass); no `unsafe` and no `extern "C"` in the workspace sources.
- QE 7.6 (dev08, MPI x4): see validation.md R3 to R5 for every number and its provenance. QE's automatic FFT grid equals ours at all six scanned cutoffs.
- One honest correction recorded along the way: the constant 2.56e-6 Ha offset is consistent with the mesh-origin segment of the G=0 quadrature (predicted 2.548e-6 Ha); QE's internal procedure was not read.

## What is not established (read before approving)
1. **Pseudopotential provenance is unverified.** `Si.pz-vbc.UPF` has blank author, date, generator input and license. It was used with the user's approval for internal validation only; it is not redistributed and not in the repository. A published, verifiable Si NC LDA file would remove this caveat; the engine needs no change for it.
2. **One cell, one pseudopotential, Gamma only.** No other element, cell, k-point or spin state was run. The three real-world UPF deviations found on the first real file show that a second pseudopotential will probably expose more.
3. **Density comparison is partial**: 12 sampled Fourier coefficients and aggregates (QE's full file was not transferred); the full-file mode is implemented and tested on synthetic data only.
4. **Orthogonality** is verified inside the engine's tests but was not part of the QE comparison (neither side exports overlaps).
5. **Performance** was not optimized or benchmarked beyond wall times in the validation record (1 to 6 s per run on one thread); GPU and parallel work are out of scope.
6. **The reference runs came through relayed sessions.** Their provenance records are attestations by the runner (commands, versions, hashes), checked by SHA-256 on receipt; they are not cryptographic proof of execution.
7. **Real-format compatibility of the verifier** was only established on one real QE output (a regression test protects it).

## Decisions that belong to the user
- Engine license (none selected; no GPL code was copied; QE source was never read).
- Whether to accept the MVP with the pseudopotential caveat above, or to first obtain a verifiable pseudopotential.
- When to start P7 (k-points, PBE through a library, GPU backend design). Nothing of P7 has been started.

## Process notes
25 delegated Codex runs are recorded in `artifacts/sessions.json` (19 done, 3 failed, 2 blocked, 1 killed); every failure had a recorded cause (service capacity error, an orchestrator command error, wrapper bugs found by real use). Code review found and removed four defects that the agents' own tests did not catch: value-dependent shortcuts in the FFT (P1A), the real-world UPF deviations (P2C, P5C), and the verifier label mismatch (P6B).
