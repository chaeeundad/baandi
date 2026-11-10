# dev08 request: real Gamma H application (2026-11-11)

Status: implemented, opt-in, NOT yet measured on dev08. Adoption rule (AGENTS.md, user 2026-11-10): every switch that improves the dev08 median wall without breaking a correctness clause is adopted and stacked; there is no percentage threshold. Mac numbers below are observations on a shared machine (load average 4 to 8 from other agents during the runs), never acceptance evidence. Ledger rows: HAPPLY-* (attempts-ledger.md, section "H application").

Trigger: dev08 stage profile (MKL, si32 strict, new defaults): `h_apply` is the largest Davidson stage (38 percent at T=1, 33 percent at T=8); trace MPI-M1-CEILING: 66 percent of apply calls carry ONE column, calls with 65 to 68 columns carry 69 percent of the column work.

## Switches (all default off; every one is bit-identical to the default when off)

| Switch | What it does | Exactness when ON | Where |
|---|---|---|---|
| `BAANDI_HAPPLY_OVERLAP=on` | the local part (pair FFTs) and the nonlocal image (P^T X, coupling, P W into a private buffer) run as the two branches of one pool join; then Y += kinetic and Y += image in the unchanged serial order. Also for the B6 cached-projection apply. At T=1 the join runs the branches one after the other | **bit-identical** for every width and thread budget (test `apply_overlap_is_bit_identical`, widths 1/2/5/8, T=1/2/4, plain and cached-projection apply, non-finite rejection; Mac si8/si16 SCF: energy bits, density SHA and Davidson/H counts equal) | dft-physics real_gamma.rs |
| `BAANDI_HAPPLY_FOLD=on` | the odd tail column of a block (every single-column call, and the last column of odd blocks) goes through a half-length complex 3D transform: the real field is folded along z (`Z = psi(2m) + i psi(2m+1)`, grid nx x ny x nz/2), multiplied by the potential split into even and odd z samples, and unfolded with the two-real-halves identity. Even blocks and the leading pairs of odd blocks are unchanged. Inactive (pair path kept) when the last grid axis is odd | rounding class: the tail column agrees with the zero-partner pair to relative 1.0e-16 (test `apply_fold_matches_pairs`, 4 cells, widths 1/2/3/6/7, T=1/2/4, thread-count independent, fold plus overlap equal to fold alone bit for bit; unit test `fold_matches_pair_with_zero_partner`, 5 grids incl. odd x/y) | dft-numerics fused_fold.rs (new), dft-physics real_gamma.rs |
| `BAANDI_PAIR_FFT=pruned` (MKL builds only) | the pair local transform (and the pair density of direct minimisation) uses the pure-Rust mirror-slab layout `PairFusedPlan` with MKL DFTI 1D batch descriptors (one per lane buffer, contiguous lanes, unscaled) instead of the full-grid `GridFusedPlan` (MKL 3D, full zero fill). This is backlog B33 step 2 | rounding class (MKL 1D kernels instead of the MKL 3D transform); one plan per pair piece, so results stay thread-count independent. The pure build ignores the switch; its `PairFusedPlan` refactor (pluggable lanes) is bit-identical | dft-mkl adapter.rs (`Plan1dBatch`), dft-numerics fused_pair_plan.rs (`PairLanes`, `MklLanes`, `MklPairPlan`) |

The fold uses the build's 3D FFT (`FftPlan3d`: MKL DFTI in MKL builds, rustfft otherwise). On the pure build the default pair plan is already pruned, so the unpruned half-grid rustfft fold is only neutral there (Mac observation below); the target is the MKL build, where the pair path is a full-grid 3D transform with an exact zero partner.

**MKL code was only type-checked on the Mac** (`cargo clippy --workspace --features dft-numerics/mkl --all-targets -- -D warnings` exit 0; no MKL library here). Before any timing run on dev08 run the MKL unit test:

```sh
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo test --release --features mkl -p dft-numerics --lib -- mkl_lanes_match_rust_lanes fused_fold
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo test --release --features dft-numerics/mkl -p dft-physics --test real_gamma_hamiltonian
```

`mkl_lanes_match_rust_lanes` checks the MKL lanes against rustfft lanes (relative 1e-14, four grids incl. 15x20x24 and 25x9x36), repeatability and that the default `MklPairPlan` is the grid plan. If it fails, skip the `BAANDI_PAIR_FFT` arm and report the message.

## Mac observations (release, pure Rust; shared machine)

H-apply microbenchmark (`happly_widths`, cubic 16 Bohr, 32 atoms with 8 projectors each = 256 columns, 12053 rows, grid 60^3; median of 30 rounds for widths 1 and 2, 12 rounds for 16 and 68; ratio to base):

| T | width | base us | overlap | fold | fold+overlap |
|---|---|---|---|---|---|
| 1 | 1 | 2780 | 1.011 | 1.034 | 1.029 |
| 1 | 2 | 3403 | 1.000 | 1.001 | 1.000 |
| 1 | 16 | 19948 | 0.995 | 0.997 | 0.995 |
| 1 | 68 | 82482 | 1.002 | 0.995 | 1.000 |
| 4 | 1 | 2025 | **0.724** | 1.005 | **0.743** |
| 4 | 2 | 2113 | **0.695** | 1.003 | **0.697** |
| 4 | 16 | 6422 | 1.037 | 0.994 | 1.002 |
| 4 | 68 | 26023 | 1.006 | 0.999 | 1.007 |

Full-grid proxy for the MKL path (`fold_full_grid_proxy`, the SAME 3D FFT implementation, rustfft, in both arms: `GridFusedPlan` pair with a zero partner against the fold; median of 15): 48^3 0.381, 60^3 0.459, 72^3 0.520 of the pair time. On dev08 the single-column local part should drop to about half (MKL half-grid descriptor against the full-grid one).

SCF (si8 and si16, strict and standard, T=1 and T=4, two seeds, interleaved, one run each): OVERLAP bit-identical in all 32 comparisons (energy bits, density SHA, Davidson and H counts); wall si8 strict T=4 0.954 and 0.987, others within +-1 percent. FOLD: |dE| 3.2e-14 (si8 strict), 2.8e-14 (si8 standard), 4.9e-13 (si16 strict), 5.5e-13 Ha (si16 standard); SCF, Davidson and H counts identical; density SHA equal between T=1 and T=4 for every arm; wall neutral (pure build, see above).

## Cells and protocol

Same binary per build, environment-only A/B. MKL build (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) is the priority; pure build as a second arm for OVERLAP only (FOLD and PAIR_FFT are not expected to help the pure build). Cells `benchmarks/cases/si16-gamma` and `si32-gamma`, profiles strict and standard, perturbed seeds 601-610 (+-1e-7 Bohr as MPI-M1-GATE3), T=1 (taskset core 0) and T=8 (cores 0-7), paired per seed, interleaved OFF/ON, node idle (load below 0.30), uptime before and after. All current defaults (APPLY_LEAN, B6, MKL_NONLOCAL, C1 ladder, ...) stay at their defaults in both arms.

Arms (each against OFF = the three switches unset):
1. `BAANDI_HAPPLY_OVERLAP=on` alone (MKL and pure).
2. `BAANDI_HAPPLY_FOLD=on` alone (MKL).
3. `BAANDI_PAIR_FFT=pruned` alone (MKL; only if the MKL unit test passes).
4. STACK = all three on (MKL).
5. If an arm is slower alone, STACK without it.

Optional stage profile (`BAANDI_DAVIDSON_STAGES=on`) for si32 strict seed 601, T=1 and T=8, OFF and STACK: report the `h_apply` stage time.

## Correctness clauses (decisive; unchanged by the adoption rule)

1. OVERLAP (bit-identical): energy bits, density SHA, SCF/Davidson/H counts and certificate verdicts equal to OFF in every run, T=1 and T=8, both builds.
2. FOLD, PAIR_FFT, STACK (rounding class): every run converges and certifies; |E_on - E_off| within the class (strict 1e-10 Ha, standard 5e-9 Ha); density SHA reproducible run to run and equal between T=1 and T=8 for each arm (thread-count independence).
3. Report the Davidson and H counts per seed (a rounding change may move a trajectory; that is not a failure by itself).

## What to send back

Per run: wall, phase timers (`density_s`, `potential_s`, `total_energy_s`, `diagonalization_s`), `convergence.work` (Davidson and H vectors recorded plus certificate), energy (17 digits), density SHA, exit code; binary sha256, commit, build command, MKL unit test output, uptime before and after, interleave order.
