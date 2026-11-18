# Performance plan digest, 2026-11-18 (external planning agent; digest by the orchestrator)

Source: a planning-agent report written against published main e902853 (planning only, nothing executed). The full report was handed over in the conversation; this file keeps the decisions that change our work. MPI against MPI only; no cross-model comparison anywhere.

## Verdict of the plan
- Finish the pending distributed potential/energy/mixing work (M2_GRID, M2_MIX), then remove remaining replicated construction and full-grid materialization, re-measure, and only then consider fewer SCF/verification maps (issues #21, #23) through bounded experiments.
- No new dense-eigensolver campaign (REL-SLICES closed), no new Advance/bounded-subspace, no restart-keep, no generic async alltoall, no generic MKL FFT switch, no ladder removal, no direct PR-CG repeat (TRACK-A0 already tried: total H 1.5-1.6x worse, non-Aufbau failures in small-gap fixtures; issue #22 must acknowledge it).
- Honest feasibility: the required removal is 56-69 % of the current MPI wall (si32 np8: 65 %, si64 np8: 69 %). Even a 50 % total reduction leaves 1.3-1.4. With W1-W8 alone the plan puts the probability of 0.90 on every MPI row below 10 %; the architecture that could reach it is an owned grid/orbital state through the whole SCF map plus lower total map work (response strategy or a genuinely cheaper occupied-subspace solver).

## Work packages (priority order)
W0 provenance and timers; W1 M2_GRID; W2 M2_MIX; W3 init-2 / owned state (no more light-state switches); W4 remove remaining replicated glue in stop/verification/output; W5 verification cost (valid-product reuse; delay probe, trajectory class); W6 FFT plan reuse at new shapes only; W7 rank contention/binding diagnostics; W8 overlap only after a positive exposed-wait budget; W9 response preconditioning (#21, research); W10 handoff error budget (#23, research); W11 occupied-subspace minimization (#22, stopped); W12 k-point pools (separate scope); W13 hybrid MPI+threads (separate model).

## Rules worth keeping
- Noise rule: A/A floor first (4 alternating pairs, first discarded); effect must be >= 3x floor with a moving mechanism counter; each job one case and one model, load < 0.30, <= 600 s except final- jobs.
- Ten-step screens are diagnostic only; any trajectory change needs full SCF; every adopted change passes the strict 1e-10 / standard 5e-9 Ha gates, certificate, and a generality matrix (cubic, skewed, elongated, slab, MgO real PP, synthetic PPs, small and large band counts).
- Never put a QE-only result in the headline as "best open-source" until the competitor inventory (R00) closes; freeze diago_full_acc and keep one row per configuration before taking a minimum.

## Missing data the plan asks for (cheap dev08 jobs)
- R00: provenance manifest (SHA, binary SHA256, features, effective switches, binding); the repo layout of the dev08 checkout is not known to the planner.
- R01: current absolute exclusive phase tree for si64 np8, then si32 np8, then np4/np2, then si8/si16; refresh one-core and threads tables per cell and thread count (old ranges 1.84-2.25 and 1.48-1.86 have no cell mapping).
- Suspicious records to re-check: wrong-binary and truncated-stdout incidents, zero MPI timers (now fixed by M2-STAGE-TIMERS), old QE energy offset 1.273e-6 Ha/atom (resolution not established), fastest-over-seeds bias in the QE target table.

## Status of the packages on 2026-11-18 (orchestrator)
- W1/W2: implemented as BAANDI_M2_GRID and BAANDI_M2_MIX (default off, rounding class, fixed-order sums, same bits for every rank count); dev08 A/B pending.
- W3: init-2 landed (shared projector table, preflight reuse, density split, Ewald gate np >= 4); owned-row construction of the complex table not done (conflicts with bit identity of the split atomic-start Gram).
- W4: BETWEEN-MATERIALIZE landed (import check split over ranks, fused lazy assembly, lazy delivered orbitals); dev08 A/B pending.
- MPI stage timers landed (BAANDI_DAVIDSON_STAGES now fills on the M2 path).
