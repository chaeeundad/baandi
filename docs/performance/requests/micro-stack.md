# dev08 request: micro-stack of small threads-only speedups (2026-11-10)

Status: implemented, opt-in, NOT yet measured on dev08. Adoption rule (AGENTS.md, user 2026-11-10): every switch that improves the dev08 median wall without breaking a correctness clause is adopted and stacked; there is no percentage threshold. Mac numbers below are observations on a shared machine (other agents' builds and tests ran concurrently; wall noise +-20 percent), never acceptance evidence. Ledger rows: MICRO-* (attempts-ledger.md, section "Micro-stack").

## Switches (all default off; every one is bit-identical to the default when off)

| Switch | What it does | Exactness when ON | Where |
|---|---|---|---|
| `BAANDI_DENSITY_REDUCE=on` | band sum of the density over fixed 4096-point grid blocks, pool-parallel, every point still summed in ascending band order; zero-occupation contribution grids (exact +0) are cleared once and skipped | **bit-identical** (unit test `pair_and_reduce_switches`, Mac si16 T=1/T=8 seeds 601-602: energy and density SHA equal) | dft-engine band_density.rs |
| `BAANDI_APPLY_LEAN=on` | real Gamma H apply: kinetic term column by column instead of an integer modulo per element; the nonlocal add skips its repeated input scans (X is scanned once at entry, Y once at the end, same error) | **bit-identical** (test `apply_lean_is_bit_identical`, including NaN rejection) | dft-physics real_gamma.rs |
| `BAANDI_B6=on` | in-solve return certificate of the real Davidson: P^T X from the cached products of the basis columns, (P^T V) C, instead of a fresh projector product (cache kept through appends, contractions and rebuilds) | eigenpairs, density and work counts **bit-identical**; only the certificate residual values change at rounding (relative 4e-7 of a 1e-10 residual on si16); a residual exactly at the threshold could flip | dft-solvers real_davidson.rs (hook), lib.rs (trait defaults), dft-engine gamma_real.rs (forwarding) |
| `BAANDI_NL_ENERGY_REAL=on` | energy audit nonlocal term on the real Gamma path: per-band <x|V_nl|x> from one real projector GEMM of the returned ORC vectors instead of complex projector sums over the full sphere | rounding class: energy differs by 3e-14 Ha (si32 seed 601), density SHA and all counts identical (energy is not fed back) | dft-physics real_gamma.rs, dft-engine gamma_real.rs, energy.rs, scf.rs |
| `BAANDI_ATOMIC_REAL=on` | atomic-start Ritz step: H Q through the real Gamma Hamiltonian (pair FFTs, real projector GEMM) instead of the complex apply; the real state built there is the one the first solve reuses | rounding class: start vectors differ at rounding, so the trajectory changes at rounding level (si32 seed 601: energy 1.8e-13 Ha, SCF/Davidson/H counts identical; density SHA differs) | dft-engine gamma_real.rs, atomic_start.rs, scf.rs |
| `BAANDI_DENSITY_PAIR=on` | real Gamma density: two occupied bands share one inverse FFT of c_a + i c_b; each band keeps its own contribution grid (so CERT_SCREEN's per-band integrals still work); pairs fixed as (0,1), (2,3), ... and worker ranges even-aligned | rounding class (density within 1e-14 relative, thread-count independent; si16 seed 601 energy 7e-15 Ha) | dft-engine band_density.rs, dft-numerics fused.rs / fused_pair.rs |
| `BAANDI_WIDE_DOT=on` | in-block MGS dot products of the real Davidson with eight lane accumulators | rounding class | dft-solvers real_davidson.rs |

Bit-identical when ON (the architect can switch the default at once if dev08 wall is not worse): DENSITY_REDUCE, APPLY_LEAN. B6 keeps eigenpairs, densities, energies and counts bit-identical; only reported residual digits move.

## Mac observations (release, pure Rust, T=1 unless stated; shared machine)

- NL_ENERGY_REAL: `total_energy` phase si32 seed 601 0.586 -> 0.078 s (-87 percent, about 2 percent of the 23 s wall), si16 0.092 -> 0.027 s.
- ATOMIC_REAL: `potential` phase (contains the atomic start) si32 0.677 -> 0.493 s, si16 0.160 -> 0.102 s; per-step potential build is only about 5 ms on si32, so the phase is the atomic start.
- B6: certificate residuals change at rounding; counts identical; the saved product is about 2 count n bands flops per strict solve (wall effect below Mac noise).
- DENSITY_REDUCE: density phase at T=8 si16 0.041 -> 0.035 to 0.039 s (inside noise); serial before, now parallel.
- DENSITY_PAIR: halves the density inverse FFTs on occupied bands; Mac density timer inside noise.
- APPLY_LEAN, WIDE_DOT: no Mac-visible change (WIDE_DOT: the scalar MGS stays bound by L2 traffic on the M-series core; the switch is left for the x86/AVX2 case where the four-lane loop compiles to a latency-bound scalar chain).

## Cells and protocol

Same binary, MKL build (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) and pure build; environment-only A/B. Cells `benchmarks/cases/si16-gamma` and `si32-gamma`, profiles strict and standard, perturbed seeds 601-610 (+-1e-7 Bohr as MPI-M1-GATE3), T=1 (taskset core 0) and T=8 (cores 0-7). Interleave OFF/ON per seed; node idle.

Arms (each against OFF = all seven unset, with the current defaults CERT_SCREEN, NL_ENERGY=fast, MKL_NONLOCAL on):
1. Each switch alone (si32 strict and standard, T=1 and T=8, seeds 601-605) for attribution.
2. STACK = `BAANDI_DENSITY_REDUCE=on BAANDI_APPLY_LEAN=on BAANDI_B6=on BAANDI_NL_ENERGY_REAL=on BAANDI_ATOMIC_REAL=on BAANDI_DENSITY_PAIR=on BAANDI_WIDE_DOT=on` on all cells and seeds 601-610.
3. If a switch alone is slower, STACK without it.

## Correctness clauses (decisive; unchanged by the adoption rule)

1. Bit-identical switches (DENSITY_REDUCE, APPLY_LEAN): energy bits, density SHA, SCF/Davidson/H counts equal to OFF in every run, T=1 and T=8.
2. B6: density SHA, energy bits and counts equal to OFF; certificate verdicts equal (report any `residuals_certified` or SCF count difference).
3. Rounding-class switches (NL_ENERGY_REAL, ATOMIC_REAL, DENSITY_PAIR, WIDE_DOT, STACK): every run converges and certifies; |E_on - E_off| within the class (strict 1e-10 Ha, standard 5e-9 Ha); density SHA reproducible run to run and equal between T=1 and T=8 for each arm (thread-count independence).
4. Report the Davidson and H counts per seed (a rounding change may move a trajectory; that is not a failure by itself).

## What to send back

Per run: wall, phase timers (`density_s`, `potential_s`, `total_energy_s`, `diagonalization_s`), `convergence.work` (Davidson and H vectors recorded plus certificate), energy (17 digits), density SHA, exit code; binary sha256, commit, build command, uptime before and after, interleave order.
