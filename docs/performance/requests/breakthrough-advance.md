# dev08 request: BREAKTHROUGH-3 bounded subspace SCF (stage 1, one case, one core)

Status: pre-registered 2026-11-17 after the Mac stage-1 screen REJECTED the idea
as a speed change (ledger BT3-ADVANCE-CONTRACT, BT3-GROWTH-TRIGGER,
BT3-BUDGET-VARIANTS, BT3-COARSE). This request is the dev08 confirmation of the
negative result on the measurement machine (so the ledger row carries a dev08
number and the idea is not re-screened on the Mac), not a candidate for
adoption. Every job here finishes in well under 10 minutes (si16 T=1 runs take
2 to 3 s on the Mac). BAANDI only, one model (threads, T=1): no QE column.
Accuracy classes unchanged; the switch is off by default and bit-identical to
main when off.

## What changes
`BAANDI_SCF_ADVANCE=on|off` (default off; `crates/dft-engine/src/advance.rs`,
`session.rs`, `gamma_real.rs`; solver `RealDavidson::advance_reusing`,
`SolveIntent::Advance(SubspaceBudget)` in `crates/dft-solvers/src/real_davidson.rs`
and `lean_davidson.rs`). Intermediate standard-path maps of the serial real
Gamma lean path run a bounded number of Davidson correction blocks from the
kept start (`BAANDI_SCF_ADVANCE_BLOCKS`, default 2) and return the approximate
occupied subspace; the density, the energy audit and the Pulay pair are built
from it; the map is completed by the exact solve at the same density when its
inner estimate B = sum f |r|^2 exceeds 0.1 x the Coulomb SCF residual E, when
E is below the threshold (convergence candidate), when E grew
(`BAANDI_SCF_ADVANCE_GROWTH`, on) or when the occupied and empty Ritz values
cross. The first map, its redo, every completion, the stop step and the
verification are exact (ADR 0043 unchanged). Research sub-switches:
`SCF_ADVANCE_GROWTH=on|off`, `SCF_ADVANCE_ADAPT=on|off`,
`SCF_ADVANCE_SCOPE=all|fine`. Design: `docs/design/bounded-subspace-scf.md`.

## Build
Branch HEAD (commit given by the orchestrator) after merge with main; MKL build
(`--features mkl`, `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`),
one binary for every arm, sha256 and build exit code recorded.

## Job (stage 1; ONE case, ONE core)
si16-gamma f1e-3 seed 601, T=1 (`taskset -c 0`), the harness input of
`scripts/fast_screen.py --case benchmarks/cases/si16-gamma.toml --seed 601
--threads 1 --xc pbe --upf <PseudoDojo Si.upf>` (the Mac screen used the
PseudoDojo PBE file; `--upf <Si.pz-vbc.UPF>` without `--xc` is the LDA
alternative, state which). Three alternating pairs, same binary, environment:
- A: default.
- B: `BAANDI_SCF_ADVANCE=on`.
- C: `BAANDI_SCF_ADVANCE=on BAANDI_SCF_ADVANCE_GROWTH=off BAANDI_SCF_ADVANCE_ADAPT=off`
  (the best Mac variant on this seed, blocks 2).
Command: `python3 -I scripts/fast_screen.py --binary <dft-cli> --upf <UPF> --xc pbe
--case benchmarks/cases/si16-gamma.toml --seed 601 --threads 1 --pairs 3
--env-b BAANDI_SCF_ADVANCE=on --json <out>.json` (and the C arm as a second
screen with `--env-b` of C). Report per arm: converged, verified, energy,
density SHA, SCF steps, verification attempts, H vectors (recorded +
certificate), `convergence.work.scf_advance` (maps advanced, completions by
trigger, advance and completion H columns), median wall and spread.

## Decisive clauses (fixed now)
1. Identity: arm A energy bits and density SHA equal the main build's (the
   switch off is bit-identical).
2. Class: |E(B) - E(A)| and |E(C) - E(A)| <= 5e-9 Ha; B and C converged and
   verified.
3. Work: the Mac rows are confirmed when B's total H columns exceed A's by more
   than 10 percent and C's are not below A's by 5 percent or more. If C is
   below A by 5 percent or more with the wall inside the spread or better, the
   idea goes back to stage 1 on three seeds (602, 603, 605) before any
   stage-2 widening; otherwise the ledger rows stand as REJECT with the dev08
   numbers and no further job is queued.

## Expansion list (only after clause 3 reopens the idea)
Seeds 602/603/605 on si16; then si32 seed 601; then si-skew. One job per line,
each under 10 minutes.
