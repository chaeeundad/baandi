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

## QE comparison protocol

P6 tooling produces no QE reference values. Generate inputs with
`python3 -I scripts/make_qe_input.py`, and validate returned bundles with
`python3 -I scripts/verify_qe.py`. The authoritative
[bundle format and comparison rules](../tests/reference/README.md) specify exact
CLI options, parser patterns, required provenance, and report semantics. The
remote execution protocol remains [reference-environment.md](reference-environment.md).

A bundle contains the executed `qe.in`, verbatim `qe.out`, and
`provenance.json` with command, printed version, host, MPI/OMP counts, UTC date,
input/UPF SHA-256, measured wall time, runner, notes, and fixture marker. Input
hash, UPF hash, and printed version must match before any comparison.
FORMAT FIXTURES always refuse comparison. No real reference is fabricated or
inferred from an internally consistent synthetic test.

Conditions must match: Bohr cell and positions (1e-12 Bohr), cutoffs in Ry,
electrons, bands, plane waves, fixed occupations, full complex Gamma sampling,
and reported LDA-PZ functional. Explicit FFT dimensions must match both our
results and QE output. Automatic QE grids are recorded beside ours with a
warning and require convergence-study interpretation. QE's density cutoff is
four times its wavefunction cutoff. Geometry uses 17 significant digits with
our matrix columns transposed to QE vector rows. The UPF's own functional is
used; `input_dft` is not set.

The total-energy tolerance is 1e-4 Ha/atom; eigenvalues and gap use 1e-3 Ha;
electron count uses 1e-8 relative. Decomposition terms use 1e-4 Ha and map
one-electron to kinetic + local + nonlocal, plus Hartree, XC, and Ewald directly.
An explicit documented convention difference may make the decomposition
informational, but cannot exempt total energy. The project also requires
orthogonality 1e-8 and eigen residual 1e-7 Ha. The verifier checks the exported
occupied eigen residual; it explicitly marks orthogonality unassessed because
these files lack an overlap matrix. It does not claim density comparison.
All used tolerances and nondefault overrides appear in reports. Tolerances
are never loosened to pass; trace units, FFT normalization, G=0, projectors,
double counting, and XC parameterization before changing physics.

Exit codes: 0 comparison passes; 1 usage/our-results/I/O error; 2 provenance
incomplete or mismatched; 3 unsupported/missing parser record or conditions
mismatch; 4 format fixture; 5 completed comparison exceeds a tolerance.
JSON and Markdown are written after all prerequisite checks, with provenance,
units, absolute differences, maxima and attaining items. For the candidate
Si digest, every report states "Pseudopotential provenance is unverified".
Only real returned, provenance-complete runs may later populate reference data.

## Convergence studies

```sh
python3 -I scripts/convergence_scan.py --input examples/si_lda.toml \
  --dft-cli target/release/dft-cli --out-dir runs/scan --scan all \
  --emit-qe-inputs --pseudo-dir /remote/pseudo --upf-name Si.pz-vbc.UPF
```

Supply a valid input with an accessible pseudopotential path; the example path
is a placeholder. The driver does not alter the source. Each point receives a
derived TOML, isolated output directory, stdout/stderr files, and a timed
`dft-cli run INPUT --output-dir DIR` command. Existing scan directories are
refused to prevent reuse of stale results. A nonzero exit, timeout (default
14400 seconds, `--timeout` override), or uncertified result aborts the scan and
records a failed point in `manifest.json`. Partial runs do not produce a success
table. Python 3.9 is supported without a TOML dependency: changed tables require
documented bare keys and single-line values; pseudopotential paths require
single-line quoted strings. Unrelated text is preserved, relative pseudo paths
are made absolute, and the real CLI validates derived TOML.

The independent series are:

- `ecut`: 20,25,30,40,50,60,80 Ry by default (`--ecuts` overrides). Remove any
  grid override so each cutoff gets its own automatic alias-free grid. Other
  settings remain those of the source input.
- `grid`: at the source cutoff, first run the automatic grid. Enumerate the
  density sphere using the first result's reciprocal matrix and cutoff (density
  radius squared = 4*ecutwfc_ry, relative boundary allowance 1e-12, matching the
  engine). Its maximum absolute Miller indices give the just-sufficient odd
  dimensions `2*max_abs_m+1`. Also run two larger grids, rounding each component
  of 1.25 and 1.5 times the automatic dimensions upward to a 2/3/5-smooth integer.
  Explicit minimum may equal automatic; retain both to exercise the override.
  `--grids 24x24x24,30x30x30` substitutes explicit candidates after automatic.
- `tolerance`: density residual 1e-5,1e-6,1e-7,1e-8,1e-9 Bohr^(-3/2).
  Set both energy tolerance in Ha and eigen tolerance in Ha to 0.1 times each
  density value. Cutoff and optional grid remain the source settings. These
  solver thresholds are recorded scan parameters, not relaxed acceptance gates.

`scan.json`, `scan.csv`, and `scan.md` record point ID, parameter, total Ha,
Ha/atom, signed total-energy difference to the finest point of that series,
HOMO-LUMO gap eV (null without an empty band), iterations, and CLI wall seconds.
Finest means largest cutoff, largest grid volume, or smallest density tolerance;
JSON/CSV retain its point ID and ordering rank. For arbitrary nonnested grids,
volume is only an ordering convention, not a proof of quadrature improvement.

`manifest.json` records source input SHA-256, invocation, each exact subprocess
argument list, derived input SHA-256, result SHA-256, exit status, elapsed wall
seconds, and the binary's recorded `tool.git_commit` (including `unknown` if
that is what the binary reports). A commit identifies the build's HEAD and does
not certify a clean worktree. With `--emit-qe-inputs`, every point emits both
explicit-our-grid and automatic-QE-grid inputs plus their sidecars. Manifest
entries bind each input hash to its matching results path/hash and expected UPF
digest. No QE execution or comparison occurs during a scan.
