# Validation record

Only measured results are recorded here, each with its command, inputs and checksums. Nothing is copied from memory or from a reference that was not actually run. QE reference numbers will be added only from real pw.x runs returned with their provenance (see docs/reference-environment.md).

## Pseudopotential provenance (applies to every run below)
`Si.pz-vbc.UPF`, sha256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`, 74552 bytes. Author, date, generator input and license are **blank or unknown** in the file; provenance is **unverified**. Use approved by the user on 2026-11-02 for internal validation only, no redistribution. The file is not in the repository. Facts: docs/candidate-si-upf.md.

## R1: Si diamond primitive cell, Gamma only, LDA-PZ, 30 Ry (first real-file run)
- Engine commit: `4fc299d3215cbfc7248105d0be5a886591f69009` (release build, Rust 1.99.1, single thread, faer sequential).
- Command: `./target/release/dft-cli run ../rust-dft-wt/runs/si_30ry/input.toml`
- Input sha256 `30d9dc99f2439f3997d9583bf67762d32b1397e6137300c0edb3a3abc2dcf2bd`; a = 5.43 Angstrom, atoms at (0,0,0) and (1/4,1/4,1/4); ecutwfc 30 Ry; Pulay mixing beta 0.5; tolerances: density 1e-7, energy 1e-8 Ha, eigen 1e-8 Ha.
- Basis: 749 plane waves, 5985 density G, FFT 25x25x25 (our rule, not necessarily QE's); 8 electrons, 8 bands (4 occupied, 4 empty).
- Result: converged in 11 iterations, wall time 1.2 s. Total energy -7.295530468210927 Ha = -14.591060936421854 Ry.
- Terms (Ha): kinetic 4.143257239948187; local -2.950156276146177 (including the G=0 alpha-Z part -0.145624326070663); nonlocal 1.598682517429664; Hartree 0.836046055268684; XC -2.523888138201652; Ewald -8.399471866509632. Harris-Foulkes functional -7.295530468210838 Ha (difference 9e-14 Ha); double-counting integral -1.620206350107395 Ha.
- Eigenvalues at Gamma (eV): -5.1733, 7.0378 (x3, degenerate to 1e-12 Ha), 9.1742 (x3), 10.2981; HOMO/LUMO gap 2.1364 eV.
- Observations (not claims of agreement): the triple degeneracies are preserved to 1e-12, valence bandwidth at Gamma (Gamma_1 to Gamma_25') is 12.21 eV, Gamma_2' to Gamma_25' is 3.26 eV. Gamma-only sampling of a two-atom cell is a coarse Brillouin-zone sample; the absolute energy is not comparable to k-converged literature values.
- Status: **not yet compared with QE.** See P6.

## Real-world UPF deviations found by this run (fixed)
1. `relativistic="no"` (legacy value) instead of the specification's "nonrelativistic" (P2C).
2. Numbers separated by newlines only; the `upf` crate splits on spaces (P2C).
3. Functional label `SLA  PZ   NOGX NOGC` instead of exactly `SLA PZ` (P5C).
These were invisible to synthetic fixtures; the lesson is recorded in docs/progress.md.
