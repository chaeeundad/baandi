TASK_ID: R-SCF-2
GOAL: |
  Issue #25 patch 2, part A (read docs/design/reference-scf.md and the merged patch 1 code in crates/dft-engine/src/scf/refscf* or wherever R-SCF-1 placed it, crates/dft-engine/src/refstart.rs, crates/dft-solvers/src/ref_davidson.rs): turn the reference driver into a complete, runnable SCF path that uses ONLY the new components, so that its trajectory can be compared with QE 7.6 step by step and later timed.
  dev08 status of patch 1 (main 29beda5): the density update reproduces QE's next density for si128 (np1 and np8) and si8 at every k = 1..6 within 4.9e-18 e/bohr^3 max (RMS ~1e-20), G=0 and total charge exact, dr2 equal to QE's printed 8 digits; the design note's very strict tolerance 1e-19 + 1e-13 x amplitude fails only at si128 k = 5 on one lowest-shell coefficient (1.6e-18 absolute, 3.8e-13 relative, history matrix condition 6689), which is BELOW QE's own rank-count floor (QE si8 np8 vs np1 differ by 1.7e-17..6.8e-17 in the output density and 7.3e-18..1.7e-17 in the mixed density). Treat patch 1 as passed; set the documented tolerance from these measured floors (state it), not the strict 1e-19.
TASK:
  (1) Driver integration: the reference driver (`[scf] driver = "reference"`, `BAANDI_REFERENCE_SCF=on`) uses refstart (final-cutoff start: atomic superposition density and atomic-orbital start block, owned rows), the new inner solver ref_davidson (QE-like block Davidson with the note's preconditioner and the note's threshold rule: atomic start 1e-2, update min(ethr, 0.1 dr2/max(1,N_el)), floor 1e-13, the first-iteration redo when dr2 < ethr N_el, empty-band threshold as in the note, all bands treated as occupied at the first step) and the patch 1 Broyden update; no coarse ladder, no legacy session, no legacy switches. Warm start between maps: previous orbitals as the next start block (as QE does). Stop when dr2 < the configured threshold; map the standard accuracy class to conv_thr explicitly (document; the dev08 si128 QE runs used conv_thr 1e-8 Ry). Then BAANDI's independent fixed-point verification (expose the existing verification from session.rs as a callable function with no behavior change for the legacy path) and the certificate in results.json; both charged to the run.
  (2) Input keys: accept `mixing = "broyden"` (and history length, beta) for the reference driver (the general validation rejected it on dev08, they used the "pulay" alias); keep the legacy behavior for the legacy driver.
  (3) Records per iteration (results.json `reference_scf`): dr2 (Ry and Ha), ethr requested and used, first-step redo flag, Davidson iterations, active bands per iteration (sum), H columns, reduced-problem orders, restarts, mixing history occupancy and condition, energies (total, band, Hartree, xc, double counting), wall per phase (start, solve, density, potential, mix, verification, output), peak bytes per rank. Comparable one to one with QE's printed lines (estimated scf accuracy, ethr, avg # of iterations, total energy).
  (4) Validation on the Mac: si8, si16, si32, skew, a slab or elongated cell and a synthetic PP: converges, passes the certificate, energy versus the strict legacy reference within the standard class (report the differences), rank counts 1/2/3/8 (simulated MPI) give the same trajectory within rounding; legacy default bit-identical. Compare the si8 trajectory (dr2 per step, ethr, iterations) against the committed QE si8 data (tests/reference/qe/i25-si8/, k1..k3 stdout lines) and report where it first differs and why.
  (5) docs/design/reference-scf.md: update the Driver section; docs/performance/lanes/RSCF/run-dev08-request.md: si128, si64, si216 np8 and si32/si64 one rank full runs with the reference driver (record everything above), plus the matching QE full runs already on dev08 for side-by-side dr2/ethr/iteration tables (no timing claims yet; one case per job; final- jobs 7200 s).
  Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh.
OWNED_FILES:
  - crates/dft-engine/src/   # reference driver files, refstart.rs, minimal accessors in scf.rs/session.rs (verification exposure only, legacy behavior unchanged), lib.rs
  - crates/dft-solvers/src/ref_davidson.rs   # fixes needed for integration only
  - crates/dft-cli/src/run.rs
  - crates/dft-cli/src/config.rs
  - crates/dft-cli/src/refscf_check.rs   # if shared helpers move
  - crates/dft-engine/tests/
  - crates/dft-cli/tests/
  - docs/design/reference-scf.md
  - docs/performance/lanes/RSCF/
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/research/qe-study/qe-density-mixing-spec.md, docs/research/qe-study/davidson-iteration-structure-gamma.md, docs/research/qe-study/qe-scf-cycle.md
TEST_COMMANDS:
  - scripts/check.sh full
ACCEPTANCE_CRITERIA: the reference driver runs end to end with only new components and passes the certificate on the listed systems; per-iteration records match the QE line format; rank independence; legacy bit-identical; dev08 request ready.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
