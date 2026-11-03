# Optimization log (newest first)

Every entry: profile finding, change, numbers before and after (same machine, `scripts/bench_scf.py`), equivalence evidence, rejected ideas. Rules: [README.md](README.md).

## 2026-11-03: profile of the starting point (macOS Apple Silicon, release build with symbols, 8-k-point case, `sample` for 8 s, 6445 samples)
Top of stack by category: FFT in `ndrustfft` lane transforms 48.0 percent; allocation (`malloc`, `free`, `memset`, `bzero`) 15.3 percent plus 15.0 percent attributed to `mach_absolute_time` (called inside the macOS `free` path, so it is allocation overhead too); `ndarray` iterator copies 5.9 percent; Hamiltonian apply and nonlocal projector accumulation 3.4 percent; `hypot` 2.1 percent (complex norm in a hot loop); orthonormalization 1.8 percent; GEMM 1.3 percent. About 30 percent of the time is allocation overhead, and the stacks show it comes from `ndrustfft::ifft_lane` and the per-call temporaries (`Fft::process` allocates its scratch on each call). Diagonalization is 91 to 97 percent of the engine time (R7).
Baseline (commit `6049217`, release, 1 thread): see baseline.json: Gamma 1.30 s, 4 k points 6.55 s, 8 k points 12.1 s (median of 3).
Planned steps (each gated by `bench_scf.py --require-identical` unless it changes the algorithm): O1 allocation-free FFT plan (preallocated scratch, `process_with_scratch`); O2 Hamiltonian apply without per-call temporaries and without `hypot`; O3 pruned FFT that skips lanes that are exactly zero for the wavefunction sphere (about 40 percent fewer lane transforms, bit-identical); O4 Davidson work: fewer Hamiltonian applications (tolerance schedule, subspace, restarts) only with the energy and certificate gates and measured equivalence; then re-profile.

## 2026-11-03: O1 / P7P1-FFT-ALLOC — persistent FFT scratch

Starting evidence: the pre-change profile above (6445 samples, 8 s) attributes
48.0% to FFT lanes and about 30% to allocation, including `Fft::process` scratch
and ndrustfft lane temporaries. This task uses that supplied profile and the
stored baseline at `6049217`; it does not claim a fresh pre-change profile.
Implementation worktree starts at `55c2565c68d60ec2225cd7058720095513e7b6b3`,
macOS arm64 (Apple Silicon), release build, no threading changes.

Change: plan both directions once with rustfft; reserve the maximum in-place
scratch requirement across all six plans and one maximum-axis-length lane.
Copy input to output once, gather/transform/scatter strided axes, and transform
contiguous lanes directly. Axis order remains 0, 1, 2. Inverse multiplication
by `1 / axis_length` after each axis and final multiplication by grid size
reproduce the former ndrustfft operation sequence. Forward division by grid
size is unchanged. The two full-grid work arrays are no longer needed.
Public API, manifests, lockfile and numerical tolerances are unchanged.
ndrustfft remains only as the test oracle; its manifest disposition is deferred
to Claude (see dependency inventory).

Allocation verification uses the permitted documented argument plus debug
assertions, not a counting global allocator (which would require unsafe):
all Vec allocation and planning occurs in `new`; `transform` uses only borrowed
slices, fixed-length copies, indexing, arithmetic and rustfft's
`process_with_scratch`. Scratch length is asserted against each plan's
`get_inplace_scratch_len`; lane length is asserted against the current axis.
No resizing, planning, logging or I/O occurs in a transform. The regression test
also checks scratch/lane pointers, lengths and capacities across repeated calls;
these checks alone are not claimed to count allocations inside dependencies.

Validation: the historical ndrustfft implementation is retained in the test
module. Seeded random complex inputs on 1x1x1, 1x3x5, 3x1x5, 3x5x1, 7x9x11,
13x5x17, 15x20x24, 24^3, 25^3, 27^3 and 36^3 produce bit-identical real and
imaginary components in both directions, twice per direction (zero ulp).
All existing tests are unchanged. `cargo test -p dft-numerics -- --nocapture`
exited 0: 11 tests passed in 0.55 s; no heavy cases need ignoring. Maximum
relative round-trip error was 1.7644611924521969e-15; maximum relative Parseval
error was 3.5528748979513846e-15 (both dimensionless).

Rejected/deferred: collapsing inverse normalization into one scale would change
rounding; retained the original sequence. A counting allocator was not added
because this ticket requires no unsafe. Pruning, threading and Hamiltonian or
solver changes remain outside O1. README figures and artifact JSON are outside
this ticket's explicit owned files and were not modified.

Verification commands (all exit 0):

```sh
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
cargo test -p dft-numerics -- --nocapture
cargo test --release -p dft-numerics legacy_bitwise_equivalence_and_workspace_reuse -- --nocapture
scripts/check.sh dft-numerics
scripts/check.sh
cargo build --release -p dft-cli
python3 -I scripts/bench_scf.py --binary /Users/mino/Develop/rust-dft-wt/P7P1-FFT-ALLOC/target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s --repeats 3 --baseline docs/performance/baseline.json --require-identical --record /tmp/P7P1-FFT-ALLOC-bench.json
```

Both check scripts passed build, test, fmt and Clippy with warnings denied.
Release oracle test passed in 0.01 s. Benchmark ran after local build/check jobs
finished, on macOS 27.0.1 arm64; the provided baseline identifies this machine
as an Apple Silicon Mac mini. Hardware-model sysctl was sandbox-denied.
Other jobs may run on this machine: Claude re-measures on a quiet machine at
acceptance. Ratios below are observations against the stored same-machine
baseline, not controlled speedup claims. Baseline was not modified.

| Case (30 Ry) | Before median wall (s) | After median wall (s) | Before engine / diag (s) | After engine / diag (s) | Iterations before / after | Energy difference (Ha) | Density | Observed ratio |
| --- | ---: | ---: | --- | --- | --- | ---: | --- | ---: |
| Gamma | 1.304109792 | 0.960475208 | 1.278177542 / 1.207626833 | 0.925481750 / 0.864972624 | 11 / 11 | 0 | identical | 1.36x |
| Shifted 2x2x2 (4 k) | 6.548742542 | 4.622876292 | 6.759375458 / 6.588469207 | 4.722267875 / 4.586680250 | 11 / 11 | 0 | identical | 1.42x |

Engine and diagonalization times are from the final repeat (script behavior),
not medians. After wall samples (s): Gamma `[1.268224375, 0.960475208,
0.929638417]`; shifted `[4.609614333, 4.622876292, 4.726622250]`.
Final energies (Ha): Gamma `-7.295530468210927`, shifted
`-7.9194234124690706`. Density SHA-256 values, identical to baseline:

- Gamma: `c4fde645365f4e7ed609486d76233335b3894540239b0443880351eea82257cd`
- Shifted: `7b7e9c5aaf457f4d6ca5de52c8f2bc9aca531e27c4f0cf7f334cad08e2515c72`

The UPF was read only from the supplied external path, never copied or sent.
These are rust-dft baseline comparisons, not new QE references. Changes are
left uncommitted for Claude's review and acceptance.
