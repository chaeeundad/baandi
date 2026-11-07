# QE source study notes

Policy: AGENTS.md, "QE source study policy". These notes are written by Claude Code from reading Quantum ESPRESSO 7.6 (read-only clone outside the repository) with the purpose of finding performance tuning points for rust-dft. Prose and mathematics only; no code, no excerpts. Each note ends with the measurable tuning hypotheses it supports, which are tested on our own code before anything is adopted.
Index (added as notes are written):
- [davidson-and-scf-thresholds.md](davidson-and-scf-thresholds.md): block Davidson structure, preconditioner, SCF threshold rule, and the tuning hypotheses H1 to H8.
- [fft-call-structure.md](fft-call-structure.md): FFT call structure of the wavefunction path (sticks, plane-fused pipeline, descriptor settings, what did not help), from the user's fastqe Phase 1 report.
- [davidson-gemm-accounting.md](davidson-gemm-accounting.md): note 4, Codex analysis of rust-dft only; GEMM shape accounting and the adopted incremental-Gram optimization (no QE source read).
- [qe-pseudo-treatment-differences.md](qe-pseudo-treatment-differences.md): ranked candidates for the BAANDI-versus-QE bulk energy offset (PW92 prefactor convention first), with BAANDI-side diagnostics.
- [qe-radial-endpoint-convention.md](qe-radial-endpoint-convention.md): how QE chooses radial integration limits, odd counts and the end rule, and the distinct variants to score against the MgO residual.
- [qe-stage-by-stage.md](qe-stage-by-stage.md): note 12, one SCF run walked stage by stage (setup, loop, H application, regterg step structure, density and mixing, parallel design, end of run) with measured QE timers (si8, si16) and the BAANDI comparison; facts that extend or correct earlier notes are listed in its section 9. Plan derived from it: [../../performance/qe-parity-plan.md](../../performance/qe-parity-plan.md).
