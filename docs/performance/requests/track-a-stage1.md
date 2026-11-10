# dev08 request: Track A stage 1, direct minimization pre-converger (pre-registered 2026-11-10)

Status when written: no dev08 data exists. The Mac stage-0 counts (docs/performance/track-a-stage0.md) already trigger the pre-registered Track A kill rule for variant A0, so this request is only to be sent if the orchestrator (or the user) decides that the wall-time observation in that file justifies a dev08 check despite the count verdict. The acceptance rules below were taken from the proposal (section 3.6) and the task brief, not derived from the Mac numbers, and must not be moved after dev08 data.

## What is measured
Branch: the Track A worktree branch (commit in the ledger row TRACK-A0). There is no CLI switch: `[scf] method = "direct_then_scf"` deliberately does not exist yet. The engine option is `dft_engine::scf::direct_min::ScfMethod::DirectThenScf` (default `Scf`, which is exactly the ordinary `run_outcome`), driven by the release-only bench test:

```
cd <worktree>
export DFT_SI_UPF=<Si.pz-vbc.UPF>  DFT_PSEUDODOJO_DIR=<dir with Mg.upf O.upf>
# pure build
cargo test --release -p dft-engine --test direct_min -- --ignored --nocapture --test-threads=1
# MKL build
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo test --release -p dft-engine --features dft-numerics/mkl --test direct_min -- --ignored --nocapture --test-threads=1
```
Environment of the bench: `TRACK_A_CASES` (comma list), `TRACK_A_SEEDS` (perturbed starts, every atom moved by a seeded uniform +-1e-7 Bohr per axis as in MPI-M1-GATE3; 0 = unperturbed), `TRACK_A_THREADS` (worker budget), `TRACK_A_GTOL` (minimizer gradient tolerance, default 1e-6 Ha; do not change for acceptance), `TRACK_A_TRACE=1` (per-iteration strict and hand-off records). Each output line gives, per case and seed: minimizer iterations, energy evaluations, H columns and H-equivalents per occupied band, final projected gradient, stop reason, hand-off strict iterations and H columns (recorded + certificate), strict-only iterations and H columns (recorded + certificate, plus atomic-start columns separately), energy differences (Ha) of the minimizer and of the certified hand-off result to the strict-only result, certification of both paths, and wall times; a second line gives FNV-64 hashes of the final densities (hand-off, strict, minimizer).

## Cells
- si16 and si32 strict (`TRACK_A_CASES=si16,si32`), perturbed seeds 601-610.
- Non-silicon: `syn-skew2`, `syn-fcc`, `syn-lowsym6` (synthetic multi-projector fixtures from crates/dft-engine/tests/common) and `mgo8` (rocksalt MgO conventional cell, PseudoDojo PBE with NLCC on O, 35 Ha), seeds 601-610.
- T=1 (taskset core 0) and T=8 (cores 0-7), `MKL_NUM_THREADS=1`, node idle (load below 0.30), pure and MKL builds separately.

## Pre-registered acceptance (all required)
1. Counts first (the kill rule, unchanged): Track A stops if the median total H columns (minimizer + hand-off recorded + hand-off certificate) exceeds the median strict-only total H columns (recorded + certificate) on si32 and on at least two non-silicon fixtures. Atomic-start columns are reported for both paths but not used in the rule.
2. Accuracy: on unique-solution fixtures the certified energy of the hand-off path is within 1e-12 Ha of the strict-only (SCF-only) path at every seed; the certificate passes on both paths at every seed; no new convergence failure. A fixture known to have several fixed points (docs/validation/scf-fixed-point-multiplicity.md) is reported separately and is not a waiver for any other.
3. Determinism: the hand-off density hash and the minimizer energy bits are identical between T=1 and T=8 within each backend.
4. Wall (second, only if 1 to 3 pass): median wall of minimizer + hand-off below the strict-only median at T=1 and T=8 on si32 strict; no non-silicon fixture slower by more than 2 percent at the median.

A failure of 1 stops Track A for variant A0 (record in the ledger, no further stage). A failure of 2 or 3 rejects the stage. Variants A1 and A2 (proposal 3.4 point 7) need their own registration.

## What to send back
The complete bench output per build and thread budget, binary or test-executable sha256, commit, build commands, uptime before and after, and the TRACE output for si32 seed 601 at T=1.
