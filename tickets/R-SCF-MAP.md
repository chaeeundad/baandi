TASK_ID: R-SCF-MAP
GOAL: |
  Issue #25 patch 2, part B: the "same input, compare output density" test. From a common input density (QE's rho_in(k) of its own trajectory), compute BAANDI's output density rho_out with its existing Hamiltonian/density kernels and the inner-solve threshold QE used at that step, and compare with QE's rho_out(k). A mismatch points at potential, occupations or the inner solve; agreement followed by a different next input points at mixing (patch 1 already showed mixing agrees to <5e-18). dev08 has QE trajectories for si128 (np8), si8 (np8 and np1) at k = 1..6: /data/sig-dft/i25/{si128,si8,si8-np1}/k{k}-{mix,out}/charge-density.dat (rho_in(k) = file k-1 mix, rho_in(1) recovered as in patch 1), stdout with ethr per step; BAANDI inputs matching the QE runs: /data/sig-dft/i25/{si128,si8,si8-np1}/input.toml. The Mac has the si8 one-rank fixture k1..k3 in tests/reference/qe/i25-si8/.
TASK: (1) A CLI subcommand (for example `dft-cli refscf-map <input.toml> --qe <dir> --steps K [--ethr-from-stdout|--ethr X] [--solver legacy|reference]`) that, for each k, loads QE's rho_in(k) by Miller index into BAANDI's density representation, builds the potential with BAANDI's kernels, solves the bands (legacy Davidson by default; the new ref_davidson as an option once available on main, crates/dft-solvers/src/ref_davidson.rs) to the given threshold, forms rho_out and compares with QE's rho_out(k): max/RMS difference, G=0 and total charge, lowest 10 shells, the energy terms if QE prints them, the band energies (if derivable), plus the Hartree-like metric distance (same metric as the mixer). Report the first failing term and a split of the potential (Hartree, xc, local) differences where possible. Write the reusable QE-reader parts as shared helpers with refscf_check.rs only if needed. (2) Tests on the Mac with the si8 fixture (k1..k3), documenting the measured differences and a tolerance derived from QE's own floor (QE np8 vs np1 differs by 1.7e-17..6.8e-17 e/bohr^3 in the output density). (3) docs/performance/lanes/RSCF/map-dev08-request.md: si128 np1 and np8, si8 np1 commands (with the QE ethr per step), what to return. Defaults unchanged. Scratch inside the worktree. Do not edit scripts/delegate.sh.
OWNED_FILES:
  - crates/dft-cli/src/refscf_map.rs   # new
  - crates/dft-cli/src/main.rs   # subcommand registration only
  - crates/dft-cli/src/refscf_check.rs   # shared helpers only
  - crates/dft-engine/src/   # only small public accessors needed to build a potential from a given density and run one fixed-potential solve; no behavior change
  - crates/dft-cli/tests/
  - docs/performance/lanes/RSCF/map-dev08-request.md
  - docs/performance/attempts-ledger.md
READ_ONLY_CONTEXT: AGENTS.md, docs/design/reference-scf.md, docs/research/qe-study/qe-density-mixing-spec.md, docs/research/qe-study/qe-scf-cycle.md
TEST_COMMANDS:
  - scripts/check.sh full
ACCEPTANCE_CRITERIA: subcommand works on the si8 fixture with a documented comparison table; dev08 request ready; defaults bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json
