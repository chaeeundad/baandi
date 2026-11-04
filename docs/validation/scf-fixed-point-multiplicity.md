# SCF fixed-point multiplicity in tiny synthetic cases (syn-s2, syn-s3, syn-z14)

Date: 2026-11-04. Trigger: dev08 measurement session reported that `syn-s2` converges to different total energies in the pure Rust and the MKL build (E = 0.6304095283 Ha vs 0.6348260419 Ha, both `converged = true`, difference 4.42e-3 Ha), while the five protocol cases and the seven silicon suite cases agree to 1e-12 Ha.

## Finding
The difference is not a backend defect. `syn-s2` (one atom, 4 bohr cubic cell, 3 Ry cutoff, two coupled s projectors) has only 7 plane waves, a 5x5x5 FFT grid and 6 bands (2 occupied, 3 buffer). Its SCF map has several distinct self-consistent solutions, and which one is reached depends on rounding-level differences in the trajectory. The pure Rust build alone reproduces both dev08 values on the Mac (aarch64, faer backend):

| change to the same input (pure Rust, Mac) | converged total energy (Ha) | iterations |
| --- | --- | --- |
| unchanged (seeds 1 to 12 identical) | 0.6304095283 | 17 |
| `nbands = 2` or `extra_bands = 0` | 0.6348260419 (equals the dev08 MKL value to all printed digits) | 14 |
| `ecutwfc_ry = 3.1` | 0.6354522014 | - |
| atom x coordinate shifted by 1e-14 | 0.6354522014 | 24 |
| shifted by 1e-12 | 0.6304095283 | 13 |
| shifted by 1e-10 | 0.6348260419 | 14 |
| shifted by 1e-8 | 0.6304095283 | 15 |
| shifted by 1e-6 | 0.6354522014 | 16 |

A displacement of 1e-14 (far below any physical meaning) changes the reached solution by up to 5 mHa: three different converged energies appear. Within one build and one input the result is bit-identical (determinism contract holds); across builds that differ at rounding level (MKL versus faer FFT and GEMM) the trajectories diverge and may settle on different fixed points.

## Consequences
- The 1e-12 Ha cross-backend gate is valid only for inputs whose SCF solution is unique. It passes for the five protocol cases and all `si-*`, `si8-gamma` and the 11 other synthetic cases (largest difference 4.6e-13 Ha, si8). It must not be applied to `syn-s2`, `syn-s3`, `syn-z14`, which are recorded as "non-unique fixed point".
- `syn-s3` and `syn-z14` do not converge within 100 iterations; their backend differences (2.9e-4 and 4.6e-3 Ha) are consistent with the same behaviour but are not evidence either way.
- No MKL code change is indicated. The earlier statement "accuracy gate passed" refers to the protocol cases only.
- Product concern (not addressed here): a converged SCF can be a non-minimal stationary point. A stability or lowest-state check for the converged density would detect it. This belongs to a later ticket.

## Reproduce (Mac, pure Rust)
`python3 scripts/bench_suite.py --binary target/release/dft-cli --upf <Si UPF> --cases syn-s2 --repeats 1 --output out.json --work-dir work`, then edit `work/syn-s2/input.toml` (`extra_bands = 0` under `[scf]`, or shift the atom coordinate) and run `dft-cli run`. Evidence files: none committed; the table above is from runs of that command on 2026-11-04.
