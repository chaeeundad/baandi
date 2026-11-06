# dev08 request batch (to be relayed through the dft5090 session once its approval prompt is answered)

Prepared 2026-11-06 by the orchestrator. Nothing here has been run on dev08. Repository HEAD at preparation: see `git log` of this file. All commands run only inside /scr/mino/rustdft (approved scope: pseudopotential transfer from the official source with hash check, QE FFTW3 build). No server-wide changes, no bulk installs.

## A. Build and identity
1. Build BAANDI release with `scripts/build_native.sh` (pure Rust, target-cpu=native) and the default MKL build; record `git rev-parse HEAD`, rustc version, flags.
2. Pseudopotentials: PseudoDojo NC-SR v0.4 PBE standard (Si, O, Mg, Ni, Al) from the official source, record SHA-256 per file (docs/data/nc-datasets.json lists expected names); Si.pz-vbc.UPF SHA-256 as already recorded.
3. QE 7.6 with FFTW3 (not MKL FFT) under /scr/mino/rustdft; record configure line and version.

## B. Real Gamma path (ticket R0 of docs/performance/real-gamma-design.md)
For si8 Gamma, 1 thread, pure native build, at least 30 samples per timed site (or 30 repeated runs, medians with min/max):
- gamma_real = off versus on, strict profile and standard profile (four combinations), wall time, Davidson steps, H vectors, SCF iterations, energy; energy difference off/on (expected below 1e-10 Ha strict).
- bench cases si8-gamma, si16-gamma, and the new si32-gamma (SCF phase) in the same four combinations, T = 1 and T = 8 (if the machine is free of other jobs).
- QE reference with `K_POINTS gamma` (real wavefunction path) for si8, si16, si32: command, input, version, UPF SHA-256, conv_thr ladder 1e-8, 1e-9, 1e-10 Ry; timer resolution rule (cases below 0.3 s invalid).
- QE with `K_POINTS automatic 4 4 4 1 1 1` and symmetry ON for k444s versus BAANDI k444s with symmetry on and off.

## C. Equal-accuracy comparison
Each code against its own tight reference (QE conv_thr 1e-12 Ry; BAANDI strict tighter): delivered energy error, density difference and wall time per row, as in tests/reference/perf/p17a/final-table.md. Label standard-profile rows `standard`, never `certified`.

## D. Output format
Plain-text summaries (tables), not base64. Include the raw numbers, commands and hashes. Do not paste anything that could contain credentials.

## E. Not requested
No change to other users' jobs, no installs outside /scr/mino/rustdft, no GPU work.
