# P10C fused wavefunction FFT study

The latest integration section at the end records completion against main
`9d6fef1`, including single-thread fused density. Earlier measurements and
permission checkpoints are preserved as historical evidence, not current blockers.

Implementation based on commit `0142ebf`, derived from the prose note
`docs/research/qe-study/fft-call-structure.md`. No QE source was read. Results
below compare rust-dft implementations; they are not new QE references.

## Outcome and scope

The default `OrderedFusedPlan` preserves the old inverse **and** forward x/y/z
axis order, every per-axis inverse scale, the final inverse multiplication by
N, forward division by N, and `weight * norm_sqr / volume` density arithmetic.
It is wired into Hamiltonian workers and engine density workers. The additive
physics adapter `from_orbitals_fused` is tested against the existing generic
`Fft3d` API. Zero-weight bands skip transforms after input validation; the
ascending-band density reduction and existing thread budget are unchanged.

The requested z-stick/xy-plane pipeline is available as `PlaneFusedPlan` for
explicit backend experiments. Its inverse z/y/x order changes bits. Measured
errors meet the 1e-13 normwise and 1e-12 Ha SCF bounds, but the density hashes
change. It is not the default. A temporary measurement build selected it by
aliasing the imports in Hamiltonian and density; those aliases were restored.

## Fused operator API

See [ADR 0026](../decisions/0026-fused-wavefunction-operators.md) and rustdoc in
`crates/dft-numerics/src/fused.rs`. A plan is constructed from a grid and the
ordered sphere indices. `apply_local` writes P V psi in that same order;
`accumulate_density` adds weight times squared modulus divided by volume.
Real arrays use C order `(x*ny+y)*nz+z`. Potentials are in Ha, volumes in Bohr^3,
and densities in electrons/Bohr^3. Coefficients describe the periodic part of a
Bloch state. The constructor supports unsorted, duplicate (last-write-wins),
empty and dense supports. Invalid shapes/lengths and invalid density weight or
volume fail before writing. Finite coefficient/potential checks remain in the
physics caller. Zero weight preserves even signed-zero accumulator bits.

Descriptors, lane storage and scratch are per-worker; there is no allocation,
planning, logging, I/O, internal threading or global cache in either kernel.
The vendor stage-2 implementation can batch sticks with its own stride and
spacing; the trait does not expose a Rust-specific buffer layout. No MKL files
are touched. The only crate-root edit is the orchestrator-approved export.
`docs/architecture/fft-backends.md` was absent during the first run. Main now
provides it; its appended "Fused operator API" section documents integration.

## Arithmetic-preserving plane fusion

A symmetric pipeline cannot simply reverse the existing inverse axis order
and claim bit identity. The final default uses an asymmetric decomposition:

- Inverse x on compact (y,z) sticks, then inverse y only in active z planes.
- For each real y, inverse z into a single xz plane. Multiply V there, then
  forward x immediately, retaining only requested x columns in a compact
  intermediate.
- Forward y in those columns, then forward z only on requested (x,y) sticks.

No full real-space wavefunction grid is formed, copied, scattered or zeroed.
The compact initial sticks are reset, and small lane buffers are filled as
needed. Inverse and forward intermediates become dense only at dense support.
Density reuses the inverse prefix and streams the final z lanes directly into
its output. The initial ordered prototype used a full real grid; it was
rejected in favor of this decomposition. The final forward z pass transforms
contiguous storage directly, avoiding a needless gather/scatter.

The symmetric experiment uses contiguous batched z transforms, active-column
y transforms, and an x pass into a distinct xy plane. Only active x columns
are gathered on the return pass; inactive columns stay zero. At >60% active
columns it evaluates all y lanes, still as separable 1D Rust FFTs. Tests compare
that dense control with pruning; no claim of a vendor-style native 2D kernel
is made.

## Verification protocol

Oracle tests retain `FftPlan3d` full inverse for density and pruned transforms
for local application. Two seeded complex vectors per support exercise scratch
reuse; coefficient order is reversed. Gamma and shifted supports cover 24^3,
25^3, 27^3, 36^3, 15x20x24 and 7x9x11, with empty/zero-radius, small, medium,
nearly full and full support. Active-column fractions range from a single
column (about 3-7%) to 100%. Additional tests cover singleton axes, duplicates,
invalid shapes/lengths, zero weights, signed zeros and errors without writes.
Complete Hamiltonian tests use Gamma and shifted bases, skew cells, and 0/3/5
synthetic projector channels over six grid shapes. Existing Hamiltonian
Hermiticity, nonlocal, SCF and thread-determinism tests remain unchanged.

Pointer/length/capacity tests cover every numerical scratch buffer through
repeated local/density calls and errors; physics/engine worker reuse tests
remain in place. Together with inspection of fixed-length slice operations and
rustfft's `process_with_scratch`, this establishes allocation-free call
structure. Pointer stability alone is not presented as a counting allocator.

The large sweep and timing experiments are ignored in debug and document their
release commands. No unsafe counting allocator, dependency or tolerance change
was introduced.

## Environment and limitations

Local host: macOS arm64, `Minhoui-Macmini.local`, Darwin 27.0.0; pinned Rust
1.99.1, rustfft 6.4.1, release with `CARGO_PROFILE_RELEASE_DEBUG=true`.
Numerical kernels use one thread; SCF benchmarks explicitly request one thread.
Other jobs share this host. Samples and min/median/max are reported rather than
asserting a quiet-machine speedup. No remote server was changed.

A fresh `sample <owned-child-pid> 1 -file /tmp/p10c-evidence/sample.txt` failed
with exit 255: sandbox process inspection was denied. The owned k444s SCF
completed with exit 0. The supplied post-O4c profile (ADR 0019: FFT 54%) and
O1/O3/P7P9 history motivate the experiment; no fresh profile is claimed.

The unchanged binary failed the requested default si8 SCF with CLI exit 4:
`k index 0: not converged after 400 iterations; max residual
1.90890249000622e-4`. The benchmark wrapper exited 3. No tolerance, band count
or solver setting was changed to make it converge. Its complete Hamiltonian
is measured separately; it is not labeled a successful si8 SCF. Two distinct
synthetic cells supplement the converged Si cases. Quiet dev08 remeasurement
with pinned cores remains an orchestrator acceptance step.

Artifacts, tickets, baselines, manifests, Cargo.lock and README are unchanged
as explicitly required. Local raw logs, binaries and inputs are under
`/tmp/p10c-evidence`; persistent numerical summaries and commands follow here.

## Numerical results

Ordered kernels and complete Hamiltonians match the oracle component bits (including signed zero). For the experimental symmetric backend, maximum relative infinity-norm error was `6.12098052203517064e-16` for local application and `1.15185970029673752e-15` for accumulated density. The denominator is the maximum output magnitude, not a pointwise relative error near zero.

The independent normalized-orbital test measured integrated-charge error `5.10702591327572009e-15` electrons, coefficient orthogonality error `1.21314423183292764e-17` (dimensionless), and local-operator Hermiticity error `4.75568200259655837e-18` Ha. The existing multi-projector synthetic SCF schedule audit passed: scheduled/legacy energies `-4.6845611842778503` / `-4.6845611842778307` Ha, difference `1.9539925233402755e-14` Ha. This latter difference compares schedules, not before/after FFT implementations.

## Per-call timings

Five release samples, 64 calls/sample after warmup, microseconds per call; each entry is min / median / max. “Local” includes scatter, both transforms, multiply and gather. “Density” includes scatter, inverse and weighted accumulation. Planning and allocation are outside timing.

| Shape | Local oracle | Local ordered | Local symmetric | Density oracle | Density ordered | Density symmetric |
| --- | --- | --- | --- | --- | --- | --- |
| [24, 24, 24] | 64.360 / 64.522 / 68.343 | 53.612 / 58.830 / 60.612 | 54.225 / 54.383 / 54.979 | 55.565 / 56.213 / 61.466 | 33.824 / 34.753 / 37.369 | 32.150 / 33.932 / 37.173 |
| [25, 25, 25] | 144.832 / 145.975 / 146.316 | 136.457 / 138.511 / 139.374 | 138.025 / 142.268 / 142.819 | 126.268 / 129.367 / 130.613 | 73.433 / 78.674 / 91.156 | 73.166 / 81.850 / 104.347 |
| [27, 27, 27] | 178.392 / 182.131 / 185.055 | 165.561 / 170.196 / 171.560 | 168.157 / 168.941 / 174.153 | 157.166 / 157.747 / 162.579 | 87.879 / 90.135 / 90.272 | 94.172 / 94.280 / 94.386 |
| [36, 36, 36] | 397.361 / 399.018 / 402.101 | 371.878 / 373.790 / 380.460 | 381.930 / 388.858 / 393.689 | 341.971 / 343.839 / 355.774 | 196.115 / 196.438 / 202.230 | 205.569 / 214.104 / 230.283 |
| [15, 20, 24] | 36.081 / 38.737 / 40.272 | 29.094 / 30.097 / 34.799 | 32.328 / 34.663 / 41.220 | 34.695 / 37.934 / 39.196 | 17.713 / 20.521 / 21.234 | 24.359 / 25.872 / 27.725 |

## Complete Hamiltonian timings

Five paired samples, 128 applications/sample. The old path includes input/output finite checks, kinetic and nonlocal terms; the new path uses the public application method. Si fixtures read the real UPF but use a deterministic synthetic local grid potential (not a converged SCF potential). The skew fixture has three synthetic Hermitian projector channels. All outputs matched bitwise. Units: microseconds/application, min / median / max.

| System | Grid / plane waves | Before | Ordered |
| --- | --- | --- | --- |
| Si-Gamma | [25, 25, 25] / 749 | 152.949 / 165.853 / 170.806 | 144.262 / 146.687 / 151.815 |
| Si-shifted | [25, 25, 25] / 754 | 153.074 / 156.880 / 170.169 | 144.718 / 146.576 / 166.845 |
| Si8 | [36, 36, 36] / 2969 | 490.139 / 505.183 / 558.246 | 454.986 / 456.658 / 522.735 |
| synthetic-skew-3-projectors | [15, 20, 24] / 29 | 27.302 / 28.786 / 30.131 | 22.297 / 23.722 / 24.536 |

## Full SCF: alternating unchanged/final binaries

Five runs of each binary per system, reversing A/B order on alternate pairs. One thread. The comparison checks density hashes and all bits of energies, eigenvalues, basis metadata and convergence histories excluding timings. All five systems match. At this measurement point the single-thread engine still selects the generic density fallback (which now skips zero weights); fused density workers are exercised by the separate thread-budget test.

| System | Before wall s (min/median/max) | After wall s (min/median/max) | Grid / occupied-or-buffer bands | E (Ha) | SCF steps |
| --- | --- | --- | --- | --- | --- |
| gamma | 0.156 / 0.160 / 0.182 | 0.138 / 0.148 / 0.420 | [25, 25, 25] / 8 | -7.2955304682107718 | 10 |
| k222s | 0.522 / 0.527 / 0.555 | 0.479 / 0.486 / 0.500 | [25, 25, 25] / 8 | -7.919423412469067 | 10 |
| k222 | 0.974 / 0.992 / 1.022 | 0.893 / 0.904 / 0.999 | [25, 25, 25] / 8 | -7.8312300571211448 | 10 |
| synthetic-skew | 0.005 / 0.005 / 0.006 | 0.005 / 0.005 / 0.006 | [8, 5, 8] / 5 | -1.4978409056181081 | 15 |
| synthetic-slab | 0.005 / 0.005 / 0.005 | 0.005 / 0.005 / 0.005 | [5, 5, 12] / 5 | -1.671232045348878 | 11 |

| System | Final max eigen residual (Ha) | Max relative charge error, all steps | Certificate density residual | Certificate occupied eigenvalue change (Ha) |
| --- | --- | --- | --- | --- |
| gamma | 3.904607423e-09 | 1.953992523e-14 | 6.932080507e-08 | 3.221221473e-08 |
| k222s | 1.107173661e-09 | 1.953992523e-14 | 7.951406994e-10 | 1.419280624e-09 |
| k222 | 1.714299600e-09 | 1.953992523e-14 | 3.663697094e-10 | 3.108434843e-09 |
| synthetic-skew | 4.182257197e-10 | 3.552713679e-15 | 4.860253962e-10 | 1.922229764e-09 |
| synthetic-slab | 1.506730781e-10 | 1.776356839e-15 | 1.142370925e-08 | 4.625753830e-08 |

Certificate density residuals follow the engine's volume-weighted L2 convention (Bohr^-3/2). No new claims of physical accuracy follow from these numerical equivalence checks.

## Actual stage-grouping experiment

Test-only grouping advances all bands through each stage, reuses one potential load per plane point across the group, and retains separate plan scratch per band. This is distinct from the integration benchmark's simpler scheduling-group control. Five samples, four groups of 32 total bands per sample; every grouped output is bitwise identical to the single-band symmetric operator. Units: microseconds/band, min / median / max.

| Grid | Group size | Time |
| --- | --- | --- |
| [25, 25, 25] | 1 | 166.276 / 168.793 / 177.718 |
| [25, 25, 25] | 2 | 166.429 / 168.667 / 171.550 |
| [25, 25, 25] | 4 | 165.382 / 167.179 / 170.846 |
| [25, 25, 25] | 8 | 163.291 / 164.629 / 166.206 |
| [25, 25, 25] | 16 | 164.424 / 166.079 / 169.485 |
| [25, 25, 25] | 32 | 173.422 / 176.678 / 185.627 |
| [36, 36, 36] | 1 | 438.259 / 441.583 / 448.591 |
| [36, 36, 36] | 2 | 436.351 / 444.050 / 463.976 |
| [36, 36, 36] | 4 | 425.991 / 434.147 / 448.816 |
| [36, 36, 36] | 8 | 423.257 / 433.595 / 444.507 |
| [36, 36, 36] | 16 | 475.161 / 499.169 / 584.235 |
| [36, 36, 36] | 32 | 516.851 / 532.020 / 601.700 |
| [15, 20, 24] | 1 | 36.204 / 38.590 / 47.454 |
| [15, 20, 24] | 2 | 37.007 / 38.988 / 43.241 |
| [15, 20, 24] | 4 | 34.878 / 38.979 / 49.339 |
| [15, 20, 24] | 8 | 47.785 / 51.165 / 53.884 |
| [15, 20, 24] | 16 | 48.345 / 48.843 / 50.713 |
| [15, 20, 24] | 32 | 35.618 / 37.616 / 48.489 |

No consistent grouping win spans the three grids. Small apparent improvements overlap variability, while large groups can increase cost substantially. Grouping remains a test-only rejected variant; no many-band 3D or plane-batching production mode was added.

## Plane pruning experiment

Only <=60% active-column shares were compared with the dense separable control. Both controls use identical axis order and compare bitwise. Five samples, 64 calls/sample. Units: microseconds/local application, min / median / max.

| Grid | Active x columns | Dense | Pruned |
| --- | --- | --- | --- |
| [24, 24, 24] | 1/24 | 56.685 / 59.306 / 60.929 | 33.152 / 35.481 / 36.892 |
| [24, 24, 24] | 7/24 | 59.069 / 60.832 / 92.951 | 41.390 / 44.309 / 45.966 |
| [24, 24, 24] | 14/24 | 63.058 / 65.732 / 71.621 | 55.312 / 57.490 / 57.740 |
| [25, 25, 25] | 1/25 | 156.979 / 157.425 / 159.854 | 86.145 / 88.271 / 88.342 |
| [25, 25, 25] | 7/25 | 167.747 / 167.965 / 169.479 | 113.230 / 116.000 / 126.545 |
| [25, 25, 25] | 15/25 | 176.010 / 181.156 / 227.572 | 148.419 / 152.245 / 152.618 |
| [27, 27, 27] | 1/27 | 191.146 / 191.635 / 193.590 | 108.104 / 109.919 / 111.473 |
| [27, 27, 27] | 8/27 | 212.234 / 213.934 / 214.817 | 145.053 / 147.717 / 149.492 |
| [27, 27, 27] | 16/27 | 224.896 / 227.104 / 227.189 | 187.643 / 187.688 / 188.352 |
| [36, 36, 36] | 1/36 | 425.081 / 426.786 / 432.721 | 238.421 / 242.038 / 246.301 |
| [36, 36, 36] | 10/36 | 462.296 / 472.848 / 482.729 | 325.861 / 336.664 / 355.893 |
| [36, 36, 36] | 21/36 | 523.426 / 543.426 / 571.705 | 415.117 / 429.458 / 457.372 |
| [15, 20, 24] | 1/15 | 41.053 / 45.311 / 48.952 | 18.889 / 19.656 / 29.372 |
| [15, 20, 24] | 4/15 | 44.732 / 48.115 / 51.570 | 25.277 / 25.986 / 29.149 |
| [15, 20, 24] | 9/15 | 60.926 / 66.255 / 70.925 | 51.436 / 52.141 / 60.039 |

Pruning helps these sparse synthetic supports. This does not establish 60% as the optimal Rust threshold or support a native 2D claim. The experiment respects the requested <=60% evaluation range; no new threshold is applied to the arithmetic-preserving default.

## Symmetric alternative SCF

The temporary symmetric-Hamiltonian binary passes the energy gate on Gamma/k222s/k222 with absolute differences `2.6645352591003757e-15`, `8.881784197001252e-16`, `8.881784197001252e-16` Ha against baseline-after-o2. All three density hashes differ. Synthetic skew/slab energy differences are `0` and `2.220446049250313e-16` Ha; their density hashes also differ. Single-thread runs still used the generic density fallback. Normwise density-kernel evidence is independent in the oracle tests. The symmetric variant is not selected because the ordered implementation avoids this numerical tradeoff.

## Persistent raw timing evidence

All arrays below are actual samples; no cross-machine ratios or historical QE timings are mixed in.

<details><summary>Kernel, Hamiltonian, grouping, pruning and SCF samples</summary>

```json
{
  "kernels_us": {
    "[24, 24, 24]": {
      "old_local": [
        66.97981250000001,
        64.360015625,
        64.46745312499999,
        64.522125,
        68.343109375
      ],
      "ordered_local": [
        53.611984375000006,
        60.131515625,
        60.611984375000006,
        56.70571875,
        58.830078125
      ],
      "plane_local": [
        54.38346875,
        54.225265625,
        54.445968750000006,
        54.28125,
        54.979171875
      ],
      "old_density": [
        61.465500000000006,
        56.212890625,
        55.565109375,
        55.76171875,
        59.083984375
      ],
      "ordered_density": [
        34.604171875,
        33.82421875,
        34.75325,
        37.369140625,
        36.806640625
      ],
      "plane_density": [
        32.149734375,
        32.160796874999996,
        36.04621875,
        37.1731875,
        33.932296875
      ]
    },
    "[25, 25, 25]": {
      "old_local": [
        145.974609375,
        145.67839062500002,
        146.281890625,
        146.31640625,
        144.83203125
      ],
      "ordered_local": [
        136.45703125,
        137.57160937499998,
        138.511078125,
        139.203765625,
        139.373703125
      ],
      "plane_local": [
        142.27865625,
        138.53060937499998,
        138.025390625,
        142.267578125,
        142.819015625
      ],
      "old_density": [
        129.3671875,
        126.267578125,
        130.61264062499998,
        126.66862499999999,
        129.951171875
      ],
      "ordered_density": [
        73.43293750000001,
        76.895828125,
        91.15559375000001,
        81.89779687500001,
        78.673828125
      ],
      "plane_density": [
        104.347015625,
        81.850265625,
        76.9296875,
        73.166015625,
        83.91276562499999
      ]
    },
    "[27, 27, 27]": {
      "old_local": [
        185.05532812500002,
        178.391921875,
        182.31901562500002,
        180.04231249999998,
        182.130859375
      ],
      "ordered_local": [
        165.560546875,
        170.19596875,
        167.82682812500002,
        170.938796875,
        171.55989062499998
      ],
      "plane_local": [
        174.15298437500002,
        168.94076562499998,
        168.156890625,
        168.486984375,
        170.89907812500002
      ],
      "old_density": [
        159.1888125,
        157.74675,
        157.166015625,
        157.643875,
        162.57878125
      ],
      "ordered_density": [
        89.8125,
        87.87890625,
        90.21549999999999,
        90.134765625,
        90.272140625
      ],
      "plane_density": [
        94.2858125,
        94.25910937500001,
        94.3860625,
        94.171875,
        94.27995312499999
      ]
    },
    "[36, 36, 36]": {
      "old_local": [
        397.3606875,
        399.249359375,
        402.10092187500004,
        399.018234375,
        398.55728124999996
      ],
      "ordered_local": [
        373.79035937500004,
        373.68620312499996,
        377.107421875,
        371.878265625,
        380.45962499999996
      ],
      "plane_local": [
        393.68879687500004,
        381.93034374999996,
        389.78515625,
        388.85806249999996,
        383.34635937499996
      ],
      "old_density": [
        343.8391875,
        342.19075000000004,
        348.407546875,
        341.971359375,
        355.77409375
      ],
      "ordered_density": [
        196.431640625,
        197.09635937500002,
        196.115234375,
        202.23046875,
        196.4375
      ],
      "plane_density": [
        205.56900000000002,
        212.675125,
        219.0546875,
        230.283203125,
        214.103515625
      ]
    },
    "[15, 20, 24]": {
      "old_local": [
        38.056640625,
        38.73696875,
        40.272140625,
        36.08071875,
        40.097015625
      ],
      "ordered_local": [
        34.799484375,
        31.046234375,
        30.097015624999997,
        29.09440625,
        30.085281249999998
      ],
      "plane_local": [
        37.592453125,
        41.220046874999994,
        32.861328125,
        32.328125,
        34.66340625
      ],
      "old_density": [
        39.19596875,
        38.630218750000004,
        37.934234374999996,
        34.6953125,
        37.60221875
      ],
      "ordered_density": [
        17.712890625,
        17.958328125,
        20.52084375,
        21.208984375,
        21.234375
      ],
      "plane_density": [
        25.282546875,
        25.872406249999997,
        27.724609375,
        26.542328125,
        24.358734375
      ]
    }
  },
  "hamiltonian_us": {
    "Si-Gamma": {
      "before": [
        170.805984375,
        166.605140625,
        165.8525390625,
        152.9488984375,
        154.7701796875
      ],
      "after": [
        146.687171875,
        151.81478125,
        147.4182890625,
        145.8505859375,
        144.262375
      ]
    },
    "Si-shifted": {
      "before": [
        153.07356249999998,
        170.1689453125,
        156.8802109375,
        153.5537109375,
        159.1002578125
      ],
      "after": [
        150.33528124999998,
        146.57584375,
        144.7177734375,
        146.1064453125,
        166.84538281250002
      ]
    },
    "Si8": {
      "before": [
        505.1829453125,
        490.139,
        491.96940624999996,
        558.246421875,
        557.8841171875
      ],
      "after": [
        456.658203125,
        454.9860078125,
        458.422203125,
        522.7353515625,
        456.145828125
      ]
    },
    "synthetic-skew-3-projectors": {
      "before": [
        27.519859375,
        30.130859375,
        27.3020859375,
        28.7864609375,
        29.587234374999998
      ],
      "after": [
        22.74609375,
        23.819984375,
        24.536460937500003,
        22.297203125000003,
        23.7223359375
      ]
    }
  },
  "grouped_us": [
    {
      "shape": "[25, 25, 25]",
      "batch": 1,
      "samples_us": [
        177.718421875,
        170.921875,
        166.27571874999998,
        166.62890625,
        168.79328906249998
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "batch": 2,
      "samples_us": [
        166.429359375,
        168.7568359375,
        168.6673203125,
        171.550453125,
        168.5205078125
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "batch": 4,
      "samples_us": [
        167.1793671875,
        168.50032812499998,
        170.846359375,
        167.177734375,
        165.3818359375
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "batch": 8,
      "samples_us": [
        163.29134374999998,
        165.8798828125,
        164.62890625,
        163.7900390625,
        166.2060546875
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "batch": 16,
      "samples_us": [
        164.424484375,
        164.50553125,
        166.0791015625,
        169.484703125,
        166.974609375
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "batch": 32,
      "samples_us": [
        185.626953125,
        181.554359375,
        176.677734375,
        174.8860703125,
        173.421875
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 1,
      "samples_us": [
        445.1376953125,
        438.2587890625,
        441.5826875,
        439.01106250000004,
        448.5911484375
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 2,
      "samples_us": [
        436.3509140625,
        439.1253203125,
        444.05012500000004,
        453.3411484375,
        463.97590625
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 4,
      "samples_us": [
        442.83040625,
        425.9912109375,
        426.4368515625,
        434.147140625,
        448.816078125
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 8,
      "samples_us": [
        444.50715625,
        435.29296875,
        423.2568359375,
        426.722984375,
        433.5947265625
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 16,
      "samples_us": [
        562.6409453125001,
        584.2346953125,
        499.16861718750005,
        479.62728125,
        475.1608046875
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "batch": 32,
      "samples_us": [
        532.01953125,
        517.0872421875,
        539.0364609375,
        516.8505859375,
        601.699875
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 1,
      "samples_us": [
        47.453773437500004,
        37.640625,
        36.204429687499996,
        38.590171874999996,
        39.4733125
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 2,
      "samples_us": [
        39.8177109375,
        43.2412109375,
        37.0595703125,
        37.007484375000004,
        38.987625
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 4,
      "samples_us": [
        38.978515625,
        34.877609375000006,
        38.1796875,
        49.119468749999996,
        49.3391953125
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 8,
      "samples_us": [
        53.8841171875,
        51.73209375,
        48.121421874999996,
        47.784828125,
        51.1647109375
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 16,
      "samples_us": [
        48.843421875000004,
        49.583656250000004,
        50.71321875,
        48.774734374999994,
        48.345046875
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "batch": 32,
      "samples_us": [
        35.7181015625,
        38.16340625,
        37.615882812500004,
        35.618484375,
        48.4892578125
      ]
    }
  ],
  "pruning": [
    {
      "shape": "[24, 24, 24]",
      "columns": "1/24",
      "dense_us": [
        56.68490625,
        59.51171875,
        57.411453125,
        60.92903125,
        59.305984375
      ],
      "pruned_us": [
        35.481109374999996,
        35.3671875,
        36.8919375,
        36.557296875,
        33.15234375
      ]
    },
    {
      "shape": "[24, 24, 24]",
      "columns": "7/24",
      "dense_us": [
        59.069015625,
        60.660796875,
        60.83203125,
        73.548171875,
        92.951171875
      ],
      "pruned_us": [
        41.38996875,
        41.78125,
        45.966140625,
        45.04621875,
        44.30859375
      ]
    },
    {
      "shape": "[24, 24, 24]",
      "columns": "14/24",
      "dense_us": [
        67.880859375,
        71.62109375,
        65.326828125,
        63.057937499999994,
        65.732421875
      ],
      "pruned_us": [
        55.31184375,
        57.067703124999994,
        57.490234375,
        57.521484375,
        57.740234375
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "columns": "1/25",
      "dense_us": [
        157.425125,
        157.66340625,
        159.853515625,
        157.015625,
        156.97917187500002
      ],
      "pruned_us": [
        88.341796875,
        88.298828125,
        88.27084375,
        86.14453125,
        87.0481875
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "columns": "7/25",
      "dense_us": [
        167.96484375,
        167.747390625,
        167.961578125,
        168.09635937500002,
        169.478515625
      ],
      "pruned_us": [
        114.50065624999999,
        116.01821875,
        116.0,
        113.23046875,
        126.544921875
      ]
    },
    {
      "shape": "[25, 25, 25]",
      "columns": "15/25",
      "dense_us": [
        176.01042187500002,
        215.958328125,
        227.572265625,
        178.662109375,
        181.15559375
      ],
      "pruned_us": [
        152.294921875,
        148.99284375,
        148.419265625,
        152.617828125,
        152.244796875
      ]
    },
    {
      "shape": "[27, 27, 27]",
      "columns": "1/27",
      "dense_us": [
        191.201828125,
        191.634765625,
        191.145828125,
        193.59048437500002,
        191.8125
      ],
      "pruned_us": [
        111.4733125,
        111.376953125,
        108.10415625,
        109.919265625,
        109.3079375
      ]
    },
    {
      "shape": "[27, 27, 27]",
      "columns": "8/27",
      "dense_us": [
        213.94465624999998,
        213.319015625,
        214.817046875,
        213.93425,
        212.234375
      ],
      "pruned_us": [
        147.716796875,
        149.491546875,
        149.401703125,
        146.018234375,
        145.052734375
      ]
    },
    {
      "shape": "[27, 27, 27]",
      "columns": "16/27",
      "dense_us": [
        224.89582812499998,
        226.650390625,
        227.103515625,
        227.1315,
        227.1888125
      ],
      "pruned_us": [
        188.3515625,
        187.79946875,
        187.662109375,
        187.642578125,
        187.688140625
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "columns": "1/36",
      "dense_us": [
        425.1484375,
        432.13996875,
        432.720703125,
        426.7858125,
        425.08137500000004
      ],
      "pruned_us": [
        242.83203125,
        238.42056250000002,
        240.516921875,
        242.03842187499998,
        246.30078125
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "columns": "10/36",
      "dense_us": [
        482.729171875,
        462.295578125,
        472.8483125,
        475.42121875,
        471.727859375
      ],
      "pruned_us": [
        333.929046875,
        343.755203125,
        325.860671875,
        336.6640625,
        355.893234375
      ]
    },
    {
      "shape": "[36, 36, 36]",
      "columns": "21/36",
      "dense_us": [
        571.705078125,
        543.4264375,
        538.660796875,
        523.4264375,
        566.8522187499999
      ],
      "pruned_us": [
        455.572265625,
        429.4576875,
        457.37240625,
        415.1171875,
        428.96875
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "columns": "1/15",
      "dense_us": [
        45.3111875,
        44.268234375000006,
        48.951828125,
        46.535796875,
        41.053390625
      ],
      "pruned_us": [
        18.889328125000002,
        19.65625,
        26.718109375,
        29.372406249999997,
        19.645828125
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "columns": "4/15",
      "dense_us": [
        44.731765625,
        47.487624999999994,
        48.115234375,
        48.9798125,
        51.56965625
      ],
      "pruned_us": [
        27.949875,
        25.27734375,
        29.149078125000003,
        25.425140624999997,
        25.986328125
      ]
    },
    {
      "shape": "[15, 20, 24]",
      "columns": "9/15",
      "dense_us": [
        60.926421875,
        62.927078125,
        66.25456249999999,
        70.92512500000001,
        68.51171875
      ],
      "pruned_us": [
        52.240234375,
        52.12890625,
        52.140625,
        60.0390625,
        51.436203125
      ]
    }
  ],
  "scf": {
    "gamma": {
      "before_wall_s": [
        0.165185833,
        0.15696329200000003,
        0.18182370800000003,
        0.15566204100000003,
        0.15956995800000007
      ],
      "after_wall_s": [
        0.41973291700000004,
        0.13804745899999993,
        0.14803295799999994,
        0.14816504200000002,
        0.16313670800000013
      ],
      "before_engine_s": [
        0.144899209,
        0.151161583,
        0.174733125,
        0.150798624,
        0.152766208
      ],
      "after_engine_s": [
        0.131841292,
        0.133408042,
        0.142890791,
        0.143253376,
        0.157702001
      ],
      "density_sha256": "968fdd04bce7c2bfe30e227a9ea35e95a8971b203e20f18a412e552eff4454eb"
    },
    "k222s": {
      "before_wall_s": [
        0.5415676249999999,
        0.5265940840000001,
        0.5221292080000004,
        0.555190584,
        0.5254420829999997
      ],
      "after_wall_s": [
        0.4829527500000004,
        0.47931570899999976,
        0.4998051659999998,
        0.49030974999999977,
        0.4855184999999995
      ],
      "before_engine_s": [
        0.536015667,
        0.521837875,
        0.517551791,
        0.550206666,
        0.520543333
      ],
      "after_engine_s": [
        0.477595042,
        0.474362416,
        0.495149084,
        0.485217625,
        0.480499957
      ],
      "density_sha256": "40b373200b0be1668f781a809e81fb08fc50ded555948bf040742ea7b873874f"
    },
    "k222": {
      "before_wall_s": [
        0.9915838750000008,
        1.021929708,
        0.9973337919999992,
        0.9783327499999999,
        0.9739246669999986
      ],
      "after_wall_s": [
        0.8930074999999995,
        0.9988099170000009,
        0.8951654590000011,
        0.9040573330000008,
        0.9101023749999992
      ],
      "before_engine_s": [
        0.986575708,
        1.016698125,
        0.992296375,
        0.971728167,
        0.968816
      ],
      "after_engine_s": [
        0.888225916,
        0.993954124,
        0.889558583,
        0.899074292,
        0.904951833
      ],
      "density_sha256": "402c5135f02d806c2420c5ef79260fbab09b90a8e103a581817a431b4f04112d"
    },
    "synthetic-skew": {
      "before_wall_s": [
        0.006248250000002287,
        0.005433542000002234,
        0.005771750000000964,
        0.0051285409999977105,
        0.005320334000000315
      ],
      "after_wall_s": [
        0.005436958000000658,
        0.005366584000000785,
        0.0057037919999984865,
        0.005476374999997091,
        0.005470917000000242
      ],
      "before_engine_s": [
        0.002219458,
        0.002391834,
        0.002555626,
        0.002245709,
        0.002430875
      ],
      "after_engine_s": [
        0.002219625,
        0.002358959,
        0.0024735,
        0.002334417,
        0.002361125
      ],
      "density_sha256": "5ae53acb15baa60a0c281f63973fdbd7b2b780c7c6f0c5b56a2895cc7691ac15"
    },
    "synthetic-slab": {
      "before_wall_s": [
        0.004500583000002223,
        0.005453374999998317,
        0.004871624999999824,
        0.005056540999998305,
        0.005205166999999733
      ],
      "after_wall_s": [
        0.004606792000000581,
        0.005063040999999657,
        0.004551291999998597,
        0.005357499999998794,
        0.004636540999999994
      ],
      "before_engine_s": [
        0.001454458,
        0.001905084,
        0.001761833,
        0.001897042,
        0.001893208
      ],
      "after_engine_s": [
        0.001577874,
        0.001689791,
        0.001396624,
        0.001830501,
        0.001513542
      ],
      "density_sha256": "4d4ba2a4e5f01c95bbe096723c807c307f906f9844edb94887fca1aae2eba4e9"
    }
  }
}
```

</details>

## Reproduction and input provenance

Si UPF: read-only external `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`; SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. No UPF is committed.

Synthetic fixtures use two-valence Gaussian local ions without projectors (different from the Si projector structure), following the existing engine-test fixture formula. The exact generation and measurement helpers used here are preserved below. The generated UPF SHA-256 is `c7f2e6d00f3824cdbd17aabeb900da9fff506a1b69fb3f7fc209819012f53271`. Its cutoff is 3 Ry; lattice vectors are in Bohr. These are numerical fixtures, not physical reference elements.

<details><summary>Recreate /tmp/p10c-evidence/synthetic.py</summary>

```python
import math, pathlib, subprocess, json, hashlib, time
root=pathlib.Path('/tmp/p10c-evidence')
def array(v): return ' '.join(f'{x:.17e}' for x in v)
points=401; rc=1.; z=2.; dx=math.log1p(12/.02)/(points-1)
r=[.02*math.expm1(i*dx) for i in range(points)]; rab=[dx*(v+.02) for v in r]
local=[-4*z/(rc*math.sqrt(math.pi)) if v==0 else -2*z*math.erf(v/rc)/v for v in r]
xml=f'''<UPF version="2.0.1"><PP_INFO>SYNTHETIC Gaussian ion, no real data.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="{z}" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="0" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{array(r)}</PP_R><PP_RAB>{array(rab)}</PP_RAB></PP_MESH><PP_LOCAL>{array(local)}</PP_LOCAL><PP_NONLOCAL></PP_NONLOCAL><PP_RHOATOM>{array([0.]*points)}</PP_RHOATOM></UPF>'''
(root/'synthetic.upf').write_text(xml)
for key,cell in [('synthetic-skew',[[6,0,0],[1,5,0],[.3,.2,7]]),('synthetic-slab',[[4,0,0],[0,4,0],[0,0,10]])]:
 text=f'''[system]
name = "SYNTHETIC Gaussian ion"
lattice_unit = "bohr"
lattice_vectors = {cell}
[[system.atoms]]
species = "Xx"
fractional = [0.13,0.21,0.32]
[pseudopotentials]
Xx = "{root}/synthetic.upf"
[basis]
ecutwfc_ry = 3.0
kpoints = [[0.0,0.0,0.0]]
[xc]
functional = "lda-pz"
[scf]
threads = 1
'''
 (root/f'{key}.toml').write_text(text)
def run(binary,tag):
 records={}
 for key in ['synthetic-skew','synthetic-slab']:
  entries=[]
  for repeat in range(3):
   out=root/f'{tag}-{key}-{repeat}'; t=time.perf_counter()
   p=subprocess.run([binary,'run',str(root/f'{key}.toml'),'--output-dir',str(out)],capture_output=True,text=True)
   dt=time.perf_counter()-t
   if p.returncode: print(key,p.returncode,p.stderr); break
   r=json.loads((out/'results.json').read_text())
   entries.append({'wall_s':dt,'energy_ha':r['total_energy']['ha'],'sha256':hashlib.sha256((out/'density.bin').read_bytes()).hexdigest(),'result':r})
  records[key]=entries
 (root/f'{tag}-synthetic.json').write_text(json.dumps(records,indent=2))
 print(tag,[(k,[(r['wall_s'],r['energy_ha'],r['sha256']) for r in v]) for k,v in records.items()])
if __name__=='__main__':
 import sys
 run(sys.argv[1],sys.argv[2])
```

</details>

<details><summary>Recreate /tmp/p10c-evidence/paired.py</summary>

```python
import importlib.util,pathlib,subprocess,time,json,hashlib,statistics,struct
root=pathlib.Path('/tmp/p10c-evidence')
spec=importlib.util.spec_from_file_location('bench','scripts/bench_scf.py');bench=importlib.util.module_from_spec(spec);spec.loader.exec_module(bench)
for case in ['gamma','k222s','k222']:
 (root/f'{case}.toml').write_text(bench.case_input('/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF',case,1))
def canonical(x):
 if isinstance(x,float):return struct.pack('!d',x).hex()
 if isinstance(x,list):return [canonical(v) for v in x]
 if isinstance(x,dict):return {k:canonical(v) for k,v in x.items() if k!='timings'}
 return x
records={}
for case in ['gamma','k222s','k222','synthetic-skew','synthetic-slab']:
 entries={'before':[],'ordered':[]}
 for rep in range(5):
  for label in (['before','ordered'] if rep%2==0 else ['ordered','before']):
   output=root/f'paired-{case}-{rep}-{label}'
   cmd=[str(root/f'{label}-cli'),'run',str(root/f'{case}.toml'),'--output-dir',str(output)]
   t=time.perf_counter();p=subprocess.run(cmd,capture_output=True,text=True);wall=time.perf_counter()-t
   if p.returncode: raise RuntimeError((cmd,p.returncode,p.stderr))
   report=json.loads((output/'results.json').read_text())
   entries[label].append({'wall_s':wall,'engine_s':report['timings']['engine_total_s'],'phases':report['timings']['phases'],
     'energy_ha':report['total_energy']['ha'],'density_sha256':hashlib.sha256((output/'density.bin').read_bytes()).hexdigest(),
     'eigenvalues':report['eigenvalues'],'convergence':report['convergence'],'basis':report['basis']})
 before=entries['before'][0];after=entries['ordered'][0]
 for row in entries['before']+entries['ordered']:
  for key in ['energy_ha','density_sha256','eigenvalues','convergence','basis']:
   assert canonical(row[key])==canonical(before[key]),(case,key)
 print(case,'identical; medians before/after',*[statistics.median(e['wall_s'] for e in entries[label]) for label in ['before','ordered']],flush=True)
 records[case]=entries
(root/'paired.json').write_text(json.dumps(records,indent=2))
```

</details>

Commands (release builds set `CARGO_PROFILE_RELEASE_DEBUG=true`; put
`/opt/homebrew/opt/rustup/bin` on PATH):

```sh
# Preserve the unchanged binary before editing; do not rebuild it from the
# modified source and call it "before".
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli
cp target/release/dft-cli /tmp/p10c-evidence/before-cli
python3 -I scripts/bench_scf.py --binary /tmp/p10c-evidence/before-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 3 --baseline docs/performance/baseline-after-o2.json --require-identical --record /tmp/p10c-evidence/before.json
# After implementation:
scripts/check.sh
CARGO_PROFILE_RELEASE_DEBUG=true cargo test --release -p dft-numerics --test fused -- --include-ignored --nocapture --test-threads=1
CARGO_PROFILE_RELEASE_DEBUG=true cargo test --release -p dft-numerics grouped_band_study -- --ignored --nocapture
CARGO_PROFILE_RELEASE_DEBUG=true P10C_SI_UPF=/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF cargo test --release -p dft-physics complete_hamiltonian_timing -- --ignored --nocapture
CARGO_PROFILE_RELEASE_DEBUG=true cargo test --release -p dft-engine --test scf full_schedule_equivalence -- --ignored --nocapture
CARGO_PROFILE_RELEASE_DEBUG=true cargo build --release -p dft-cli
cp target/release/dft-cli /tmp/p10c-evidence/ordered-cli
python3 -I /tmp/p10c-evidence/synthetic.py /tmp/p10c-evidence/before-cli before
python3 -I /tmp/p10c-evidence/paired.py
python3 -I scripts/bench_scf.py --binary /tmp/p10c-evidence/ordered-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 3 --baseline docs/performance/baseline-after-o2.json --require-identical --record /tmp/p10c-evidence/final-bench.json
python3 -I scripts/bench_scf.py --binary /tmp/p10c-evidence/ordered-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 16 --repeats 3 --baseline docs/performance/baseline-after-o2.json --require-identical --record /tmp/p10c-evidence/threads16.json
git diff --check
```

The final release commands above exited 0. An initial workspace check caught
an incorrect `PlaneWaveBasis::len()` use in a newly added test (exit 1); it was
corrected to `basis.vectors().len()`, and the complete check subsequently
passed. The initial synthetic-input attempt omitted the required explicit
kpoints, so CLI exit 1 was corrected by adding Gamma kpoints before collecting
any accepted samples. No numerical tolerance was adjusted.

The 16-thread-budget benchmark runs the fused density worker path for all
three cases; energies and hashes match the one-thread baseline in every repeat.
This is a determinism check, not an FFT-internal threading mode or a
cross-thread-count speedup claim.

## Historical integration checkpoint (superseded below)

All permitted implementation files are complete and the full workspace check
passes (exit 0). The existing single-thread constructor in
`crates/dft-engine/src/scf.rs:607` creates `BandDensity` only when `inner > 1`.
Consequently the final one-thread SCF measurements above exercise the fused
Hamiltonian plus zero-occupation skipping in the generic density fallback;
they do **not** establish single-thread fused-density integration.

At that checkpoint, a permission request was pending for the following minimal wiring change, which
is outside the ticket's explicit owned files. The existing `inner` budget and
ordered reduction remain unchanged. No implementation workaround ignores the
caller-supplied `Fft3d` backend, introduces a global cache, or rebuilds FFT plans
inside each density call.

```diff
-            let band_density = if inner > 1 {
-                Some(crate::band_density::BandDensity::new(
-                    &k_grid, bands, inner,
-                )?)
-            } else {
-                None
-            };
+            let band_density = Some(crate::band_density::BandDensity::new(
+                &k_grid, bands, inner,
+            )?);
```

The follow-up authorization was subsequently granted; the connection is now
applied and the latest integration evidence below supersedes this checkpoint. The 16-thread benchmark already
exercises the fused Hamiltonian and fused density together for all three Si
cases and is bit-identical to baseline-after-o2 in all three repetitions.

The symmetric experimental build was also tested with a 16-thread budget,
activating both reordered local and density kernels: absolute total-energy
changes were `2.6645352591003757e-15`, `0`, and
`1.7763568394002505e-15` Ha for Gamma/k222s/k222, respectively; the energy gate
passed (exit 0) and all density hashes differed. It remains experimental.


## Latest integration: main 9d6fef1 and single-thread fused density

The orchestrator authorized the minimal `scf.rs` constructor change. Main was
pinned to `9d6fef14ce86dce821bd359c6f58fef0aecad8ab` (including P7PC, P10A,
P10B, P10D and the CLI/output additions). `git merge --no-commit --no-ff main`
failed with exit 128 because the sandbox cannot create
`/Users/mino/Develop/rust-dft/.git/worktrees/P10C-FUSED-FFT/ORIG_HEAD.lock`.
The original task files were backed up under `/tmp/p10c-resume/before-integration`.
A three-way content merge imported 157 upstream paths, preserved the task edits,
and combined the append-only optimization log without unresolved conflicts.
The original Git branch/index metadata remains unchanged; this is source-tree
integration, not a recorded Git merge. A separate index under `/tmp/p10c-resume`
compares the resulting files against the pinned main for review. The generated
main-relative patch is the handoff for the orchestrator's actual merge/commit.

Only the authorized constructor condition changed in `scf.rs`. Solver choices,
P7PC band tolerances, RMM-DIIS, mixing, thread counts and convergence thresholds
are inherited unchanged from main. Every k state now constructs `BandDensity`,
including an inner budget of one. It owns its plans once and reuses the existing
per-band contribution buffers and ascending reduction; those persistent arrays
remain O(bands * grid size). The single-thread path no longer creates per-band
FFT temporaries through the generic fallback. The generic physics API remains
available for other callers.

P10B introduced compile-time backend selection in main. To preserve it,
`FusedPlan` aliases `OrderedFusedPlan` normally and
`GridFusedPlan<FftPlan3d>` with `mkl`. The latter reuses the selected native 3D
backend and full-grid scratch; it is a compatibility path, not a claimed vendor
fusion. Its generic implementation is covered by local oracle and scratch
reuse tests. No changes relative to main were made to MKL files, manifests,
Cargo.lock, tickets, artifacts or README. The new ADR was renumbered to 0026
to avoid main's existing volumetric-output ADR 0021.

`cargo check --workspace --all-targets --features dft-numerics/mkl` exited 101
at an existing main error in `crates/dft-numerics/src/linalg.rs`: `Accum` and
`matmul` imports are gated off for MKL while the new `gemm_into` implementation
uses them unconditionally. That file is byte-identical to pinned main and is
outside P10C ownership. The default workspace check passes. No successful MKL
feature or linked execution is claimed; the upstream import fix and native
validation remain with the backend owner.

The unmodified pinned-main binary was built by temporarily restoring the four
P10C-modified production files to their exact main contents, building release,
copying the executable to `/tmp/p10c-resume/main-cli`, and restoring task files
in a `finally` block. The final executable is `/tmp/p10c-resume/fused-cli`.
Both use `CARGO_PROFILE_RELEASE_DEBUG=true`, pinned Rust 1.99.1 and one-thread
benchmark inputs. Tests/compilation are excluded from the paired timing run.

### Baseline provenance

Even the unmodified main fails `--require-identical` against the requested
`baseline-after-o4c.json` (exit 5): density hashes already differ, while absolute
energy differences are 0, 1.7763568394002505e-14 and 3.552713678800501e-15 Ha
for Gamma/k222s/k222. The final implementation must therefore pass two separate
gates: energy within 1e-12 Ha against that historical baseline, and bitwise
identity against unmodified current main. No baseline was edited or silently
substituted. Si8 is absent from the historical baseline and is compared to the
fresh main snapshot instead. Unlike the original task base, pinned main's
P7PC implementation converges Si8 in 15 SCF iterations.

### Final current-main comparisons

Three alternating main/fused pairs per system, reversing the order on alternate
pairs. Every stored energy, eigenvalue, basis entry and convergence-history
float (apart from timings) is bitwise equal, and density files have identical
SHA-256 values. The paired run used no concurrent local build/test command;
other users/jobs on the shared host remain uncontrolled. Times below are
seconds, min / median / max. These are observed samples, not quiet-node claims.

| Case | Main wall s | Fused wall s | E (Ha), identical | Steps |
| --- | --- | --- | --- | --- |
| gamma | 0.125023000 / 0.125359875 / 0.135435417 | 0.113186167 / 0.113771583 / 0.123638167 | -7.2955304682107691 | 10 |
| k222s | 0.400644667 / 0.404884750 / 0.406616083 | 0.357913833 / 0.359316791 / 0.361822041 | -7.919423412469083 | 10 |
| k222 | 0.755042375 / 0.757255084 / 0.757629291 | 0.668086625 / 0.670668291 / 0.674102209 | -7.8312300571211368 | 10 |
| si8 | 3.718965375 / 3.725411875 / 3.784349875 | 3.618803500 / 3.634331333 / 3.652253708 | -31.321803601031377 | 15 |
| synthetic-skew | 0.004853708 / 0.005038209 / 0.007560167 | 0.004880459 / 0.004967958 / 0.005261416 | -1.4978409056181083 | 15 |
| synthetic-slab | 0.004391750 / 0.004410417 / 0.006599292 | 0.004145750 / 0.004789917 / 0.008129250 | -1.6712320453488845 | 12 |

| Case | Median engine s: main → fused | Median density phase s: main → fused | Final eigen residual (Ha) | Certificate density residual (Bohr^-3/2) |
| --- | --- | --- | --- | --- |
| gamma | 0.118449417 → 0.107002292 | 0.010643623 → 0.003618750 | 5.402074254e-09 | 7.042873109e-08 |
| k222s | 0.397525292 → 0.353284042 | 0.042536415 → 0.014616916 | 1.473995449e-09 | 8.643077422e-10 |
| k222 | 0.749887000 → 0.665082916 | 0.084886625 → 0.029029209 | 1.484683820e-09 | 5.399170572e-10 |
| si8 | 3.710476500 → 3.622145209 | 0.104090459 → 0.052569165 | 5.299856401e-09 | 2.793067528e-08 |
| synthetic-skew | 0.001908541 → 0.001837708 | 0.000122999 → 0.000043124 | 1.606908186e-09 | 2.519428291e-09 |
| synthetic-slab | 0.001407541 → 0.001323625 | 0.000104334 → 0.000031293 | 1.062641555e-09 | 1.555704306e-10 |

The tiny synthetic slab's process wall median increased, while its measured
engine and density phases decreased. No wall-time speedup is claimed for that
case and no size-specific tuning was introduced. Si8 now has real converged
before/after SCF evidence, superseding the old-base failure above.

The first resumed synthetic run was rejected by P8G's newly required element
identity (CLI exit 1); before collecting the final samples, both binaries'
inputs were given the same explicit synthetic-test declaration:
`Xx = { path = "/tmp/p10c-evidence/synthetic.upf", element = "He",
allow_element_mismatch = true }`. This declares the element metadata override
rather than bypassing validation; the Gaussian fixture remains nonphysical.
No pseudopotential data or numerical solver parameter changed. The final
paired driver exited 0 after all six cases.

Final kernel sweep: all five release integration tests passed, including the
full support/shape sweep and timings. Three unit tests passed for scratch
stability, permanently zero inactive columns and the generic backend adapter.
Ordered output still matches every oracle bit; the symmetric maximum errors
remain 6.12098052203517064e-16 (local) and 1.15185970029673752e-15 (density).
The compatibility adapter's local and density outputs are independently
compared against the Rust oracle; its buffers survive calls and errors.

<details><summary>Final integration raw samples and binary provenance</summary>

```json
{
  "gamma": {
    "main": {
      "wall_s": [
        0.135435417,
        0.125023,
        0.125359875
      ],
      "engine_s": [
        0.121553291,
        0.118449417,
        0.118139501
      ],
      "density_s": [
        0.010643623,
        0.010611167,
        0.010661583
      ],
      "energy_ha": -7.295530468210769,
      "density_sha256": "74963159f8eb058d58220f3bbde70a806a3622ca3c568b7d7d1ef068207610b0"
    },
    "fused": {
      "wall_s": [
        0.12363816700000002,
        0.113186167,
        0.11377158300000001
      ],
      "engine_s": [
        0.109592291,
        0.107002292,
        0.106897167
      ],
      "density_s": [
        0.003594125,
        0.00361875,
        0.003632625
      ],
      "energy_ha": -7.295530468210769,
      "density_sha256": "74963159f8eb058d58220f3bbde70a806a3622ca3c568b7d7d1ef068207610b0"
    }
  },
  "k222s": {
    "main": {
      "wall_s": [
        0.4066160830000001,
        0.4048847499999997,
        0.4006446669999999
      ],
      "engine_s": [
        0.399143666,
        0.397525292,
        0.395562793
      ],
      "density_s": [
        0.042784918,
        0.042285793,
        0.042536415
      ],
      "energy_ha": -7.919423412469083,
      "density_sha256": "c1de1d09cb8dc2c53500b10ce90e7ee8dc1bf692e1a87e2f2fec59cc64315617"
    },
    "fused": {
      "wall_s": [
        0.3593167909999999,
        0.357913833,
        0.3618220409999999
      ],
      "engine_s": [
        0.353284042,
        0.352257666,
        0.354483875
      ],
      "density_s": [
        0.014736668,
        0.014543125,
        0.014616916
      ],
      "energy_ha": -7.919423412469083,
      "density_sha256": "c1de1d09cb8dc2c53500b10ce90e7ee8dc1bf692e1a87e2f2fec59cc64315617"
    }
  },
  "k222": {
    "main": {
      "wall_s": [
        0.7550423749999999,
        0.7576292909999998,
        0.7572550840000005
      ],
      "engine_s": [
        0.747134124,
        0.749887,
        0.751801124
      ],
      "density_s": [
        0.08456675,
        0.084886625,
        0.085714666
      ],
      "energy_ha": -7.831230057121137,
      "density_sha256": "a38c265d381a9b918f667fbb14d43e210d3513387c9d3c41b1533930cc2515b3"
    },
    "fused": {
      "wall_s": [
        0.6741022090000004,
        0.6680866249999999,
        0.6706682909999993
      ],
      "engine_s": [
        0.666599291,
        0.662151917,
        0.665082916
      ],
      "density_s": [
        0.031678666,
        0.028813918,
        0.029029209
      ],
      "energy_ha": -7.831230057121137,
      "density_sha256": "a38c265d381a9b918f667fbb14d43e210d3513387c9d3c41b1533930cc2515b3"
    }
  },
  "si8": {
    "main": {
      "wall_s": [
        3.7254118750000007,
        3.718965375,
        3.784349875
      ],
      "engine_s": [
        3.7104765,
        3.708267958,
        3.7691032509999998
      ],
      "density_s": [
        0.10402254,
        0.104090459,
        0.104500208
      ],
      "energy_ha": -31.321803601031377,
      "density_sha256": "484f4f279ea968d7b0466f85130c7e33552cc37595becf09d29e306359c26698"
    },
    "fused": {
      "wall_s": [
        3.6188035000000003,
        3.6343313330000004,
        3.652253708
      ],
      "engine_s": [
        3.604504125,
        3.622145209,
        3.638588958
      ],
      "density_s": [
        0.052018002,
        0.053247082,
        0.052569165
      ],
      "energy_ha": -31.321803601031377,
      "density_sha256": "484f4f279ea968d7b0466f85130c7e33552cc37595becf09d29e306359c26698"
    }
  },
  "synthetic-skew": {
    "main": {
      "wall_s": [
        0.007560167000001172,
        0.004853708000002399,
        0.00503820900000207
      ],
      "engine_s": [
        0.002146333,
        0.001908541,
        0.001894375
      ],
      "density_s": [
        0.000122999,
        0.000121295,
        0.000123249
      ],
      "energy_ha": -1.4978409056181083,
      "density_sha256": "0f7e57524f3abdc8ed78a9ab96e45b5c0970e03539cb2e6ecf709caf9477c04f"
    },
    "fused": {
      "wall_s": [
        0.005261415999999741,
        0.004880458999998893,
        0.004967958000001715
      ],
      "engine_s": [
        0.002038791,
        0.001837708,
        0.001819916
      ],
      "density_s": [
        4.3124e-05,
        4.1329e-05,
        4.3211e-05
      ],
      "energy_ha": -1.4978409056181083,
      "density_sha256": "0f7e57524f3abdc8ed78a9ab96e45b5c0970e03539cb2e6ecf709caf9477c04f"
    }
  },
  "synthetic-slab": {
    "main": {
      "wall_s": [
        0.006599292000000645,
        0.004410416999999001,
        0.0043917499999999166
      ],
      "engine_s": [
        0.001446833,
        0.001399417,
        0.001407541
      ],
      "density_s": [
        0.000105418,
        0.000104334,
        0.00010425
      ],
      "energy_ha": -1.6712320453488845,
      "density_sha256": "73ad78a4fb85b2579e1955ab0795520e745853310d6360c7a6db0e18fca794bd"
    },
    "fused": {
      "wall_s": [
        0.008129249999999644,
        0.004789917000000088,
        0.004145750000002835
      ],
      "engine_s": [
        0.001330333,
        0.001323625,
        0.001270376
      ],
      "density_s": [
        3.1543e-05,
        2.8751e-05,
        3.1293e-05
      ],
      "energy_ha": -1.6712320453488845,
      "density_sha256": "73ad78a4fb85b2579e1955ab0795520e745853310d6360c7a6db0e18fca794bd"
    }
  },
  "binaries": {
    "main-cli": "e525d5d52d962bc3b5b0849df6f614af9768489647aa5521d446917c8042a674",
    "fused-cli": "fcb073065fda3372d1f63ca2d703458b3dd4b1665a7979bafcb7f3732ddb4533"
  }
}
```

</details>

Commands for this integration run (the main binary was built before changing
its constructor, with task files temporarily restored as described above):

```sh
scripts/check.sh
CARGO_PROFILE_RELEASE_DEBUG=true cargo test --release -p dft-numerics --test fused -- --include-ignored --nocapture --test-threads=1
CARGO_PROFILE_RELEASE_DEBUG=true cargo test --release -p dft-numerics fused::tests -- --nocapture
python3 -I scripts/bench_scf.py --binary /tmp/p10c-resume/main-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,si8 --threads 1 --repeats 3 --baseline docs/performance/baseline-after-o4c.json --require-identical --record /tmp/p10c-resume/main-bench.json
python3 -I scripts/bench_scf.py --binary /tmp/p10c-resume/fused-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,si8 --threads 1 --repeats 3 --baseline /tmp/p10c-resume/main-bench.json --require-identical --record /tmp/p10c-resume/fused-bench.json
python3 -I scripts/bench_scf.py --binary /tmp/p10c-resume/fused-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 3 --baseline docs/performance/baseline-after-o4c.json --record /tmp/p10c-resume/historical-gate.json
python3 -I /tmp/p10c-resume/paired.py
python3 -I scripts/bench_scf.py --binary /tmp/p10c-resume/fused-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,si8 --threads 16 --repeats 1 --baseline /tmp/p10c-resume/main-bench.json --require-identical --record /tmp/p10c-resume/threads16.json
```

The historical strict main comparison exits 5 for the pre-existing hash drift;
the final historical energy and current-main strict comparisons exit 0. The
16-thread run is a numerical determinism check; its timing is not used because
it overlapped the final workspace check. Full workspace build/test/fmt/clippy
pass. The separate MKL-feature check exits 101 for the upstream import error
already described; it is not reported as a passing gate.
