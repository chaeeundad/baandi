TASK_ID: R-SCF-START
GOAL: |
  Part of the new SCF path (GitHub issue #25; read docs/design/reference-scf.md). User decisions 2026-11-23: the legacy SCF machinery is replaced, old code need not be preserved, bold rewrites allowed where they promise speed. This ticket builds the START of the reference SCF run: everything from process start to the first Davidson call, at the FINAL cutoff (no coarse ladder), lean in time and memory.
  Measured motivation (dev08, np8, decomposition total = per-map cost x map count x non-SCF overhead): the non-SCF overhead ratio is 1.49-1.66 at si32/si64; QE's init_run is 0.19 s at si64 np8 (0.06 s si32) while BAANDI spends about 0.98 s in init + atomic start (atomic_start alone 0.60 s at si64 np8, inside the coarse ladder 0.614 s si64, 2.778 s si128, 10.721 s si216) and 0.25-0.54 s in setup; outside-process (MPI init, binary load, output) 0.41-0.51 s at small cells; memory BAANDI/QE 9.9x at si128, 13.4x at si216 (1.93 GB vs 0.19 GB at si128). QE starts from the atomic superposition density and an atomic-orbital start block (startingwfc atomic) at the final cutoff; see docs/research/qe-study/qe-scf-cycle.md (our own words); never look for or use QE source.
TASK:
  (1) Audit the current setup path for a real-Gamma integer-occupation run under MPI (problem construction, basis and grid plans, projector tables, structure factors, local ionic potential, Ewald, atomic superposition density, atomic start block, first potential, distributed layout setup), with per-rank bytes and calling frequency, using the existing INIT/LAUNCHIO/OWNED lane audits (docs/performance/lanes/INIT, LAUNCHIO, OWNED) and the merged default-off M2_OWNED_INPUT. List what the reference path needs and what it can skip (anything for the coarse stage, replicated serial operator, full-grid materialization, legacy diagnostics).
  (2) Implement a lean start function for the reference driver in a NEW file crates/dft-engine/src/refstart.rs (do not edit crates/dft-engine/src/refscf/, owned by ticket R-SCF-1): inputs = the parsed problem; outputs = owned (distributed) operator state, the atomic superposition density at the final cutoff (normalized to N_el, G=0 = N_el/Omega), the atomic-orbital start block (owned rows, orthonormal; randomization rule documented and deterministic), the first potential; built directly in owned form (no full-length temporaries where avoidable), memory accounted per rank. Reuse existing kernels where they are good; rewrite where they waste (state why with numbers).
  (3) Validation on the Mac: the start density and start block equal (rounding class) what the existing code produces for the same final-cutoff problem on si8/si16/si32/skew and a synthetic PP, at 1/2/3/8 simulated ranks; peak bytes per rank versus the current path (counters); a short end-to-end check that the legacy driver run from this start converges to the same energy (no change of defaults).
  (4) docs/design/reference-scf.md: append a section "Start" (append only; other tickets edit other sections) and docs/performance/lanes/RSCF/start-dev08-request.md: dev08 measurement of the start path (si64, si128, si216 np8; time to first Davidson call and peak RSS, A/A floor first, crossed A/B against the legacy start, one case per job).
  Defaults unchanged; no Mac timing claims. Scratch inside the worktree.
OWNED_FILES:
  - crates/dft-engine/src/refstart.rs
  - crates/dft-engine/src/lib.rs   # module declaration only
  - crates/dft-engine/src/atomic_start.rs   # only if a reusable owned-form builder is needed; behavior unchanged for the legacy path
  - crates/dft-engine/src/gamma_real.rs   # only owned-form constructors, legacy behavior unchanged
  - crates/dft-physics/src/distributed_real_gamma.rs   # only owned-form constructors, legacy behavior unchanged
  - crates/dft-engine/tests/   # new test files
  - docs/design/reference-scf.md   # "Start" section only
  - docs/performance/lanes/RSCF/start-dev08-request.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/design/reference-scf.md, docs/research/qe-study/qe-scf-cycle.md, docs/performance/lanes/INIT/, docs/performance/lanes/LAUNCHIO/, docs/performance/lanes/OWNED/
TEST_COMMANDS:
  - scripts/check.sh full
ACCEPTANCE_CRITERIA: audit with bytes and frequencies; refstart function producing owned start state equal to the existing start within rounding at 1/2/3/8 ranks; per-rank memory counters; dev08 request; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
