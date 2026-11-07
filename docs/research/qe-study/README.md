# QE source study notes

Policy: AGENTS.md, "QE source study policy". These notes are written by Claude Code from reading Quantum ESPRESSO 7.6 (read-only clone outside the repository) with the purpose of finding performance tuning points for rust-dft. Prose and mathematics only; no code, no excerpts. Each note ends with the measurable tuning hypotheses it supports, which are tested on our own code before anything is adopted.
Index (added as notes are written):
- [davidson-and-scf-thresholds.md](davidson-and-scf-thresholds.md): block Davidson structure, preconditioner, SCF threshold rule, and the tuning hypotheses H1 to H8.
- [fft-call-structure.md](fft-call-structure.md): FFT call structure of the wavefunction path (sticks, plane-fused pipeline, descriptor settings, what did not help), from the user's fastqe Phase 1 report.
- [davidson-gemm-accounting.md](davidson-gemm-accounting.md): note 4, Codex analysis of rust-dft only; GEMM shape accounting and the adopted incremental-Gram optimization (no QE source read).
- [qe-pseudo-treatment-differences.md](qe-pseudo-treatment-differences.md): ranked candidates for the BAANDI-versus-QE bulk energy offset (PW92 prefactor convention first), with BAANDI-side diagnostics.
