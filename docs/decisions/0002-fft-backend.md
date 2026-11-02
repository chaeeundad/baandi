# ADR 0002: FFT backend selection

Status: Accepted (2026-11-02, Claude Code). Revisit if P1 tests or benchmarks contradict the unverified items below.

## Context
The engine needs complex-to-complex 3D FFTs on small grids (tens of points per axis) with sizes that are products of small primes (2, 3, 5). FFT must sit behind a narrow trait so it can be swapped (see ADR 0001). Accelerate is only a BLAS/LAPACK backend here.

## Facts checked (crates.io API, 2026-11-02)
| Crate | Latest | Last release | License | Note |
| --- | --- | --- | --- | --- |
| rustfft | 6.4.1 | 2025-09-18 | MIT OR Apache-2.0 | 1D, arbitrary lengths, ~34.7M downloads |
| realfft | 3.5.1 | 2026-10-26 | MIT | real-to-complex wrapper around rustfft; not needed for complex Gamma-only |
| phastft | 0.4.1 | 2026-07-31 | MIT OR Apache-2.0 | power-of-two lengths only (per its README); too restrictive for 3D DFT grids |
| fftw (binding) | 0.8.0 | 2023-04-29 | non-standard (GPL-family via FFTW) | stale; GPL conflicts with possible closed distribution |

## Decision
1. Default backend: `rustfft` (permissive license, arbitrary lengths, widely used).
2. Implement the 3D transform ourselves as three batched 1D passes over axes, in our own crate, with plan and scratch reuse. rustfft provides no 3D API.
3. Keep FFTW and Apple vDSP out of the default build. A vDSP or FFTW backend may be added later behind the same trait; FFTW only after a license review.
4. External crates are fetched by Claude (Codex sandbox has no network) and reviewed for license before being added.

## Unverified, to be checked in P1
- NEON acceleration on aarch64 for the pinned rustfft version (check by build and micro-benchmark, do not assume).
- Last release is 14 months old; acceptable for a stable numeric kernel, but record the exact pinned version and checksum in Cargo.lock.
- Accuracy of 3D round trip and Parseval with non-power-of-two sizes such as 15, 20, 24, 27.

## Consequences
Single permissive dependency, no GPL exposure. Cost: we own the 3D axis loop and its normalization tests.
