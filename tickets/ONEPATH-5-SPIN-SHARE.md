TASK_ID: ONEPATH-5-SPIN-SHARE
GOAL: |
  A collinear spin run with zero total moment and equal channel occupations (fixed sector n_up = n_down, initial moment seed 0, no constraint that breaks the symmetry) solves the same Hamiltonian twice per SCF map (spin up and spin down), because both channels see the same potential. Make the solver detect the exactly degenerate case and solve ONE channel, reuse its orbitals/density for the other channel, so the spin wall time falls toward the unpolarized cost. Keep everything bit-identical in meaning: the spin results (energy, magnetization zero, per-channel orbitals, eigenvalues) must equal those of the two-channel run to the rounding class, and any symmetry breaking (non-zero seed, different channel counts, spin-dependent potential difference above a stated threshold, fractional occupation differences) falls back to the existing two-channel solve with a recorded reason.
OWNED_FILES:
  - crates/dft-engine/src/collinear.rs
  - crates/dft-engine/src/scf.rs   # only if the channel solve loop lives there
  - crates/dft-engine/tests/spin_share.rs   # new
  - crates/dft-numerics/src/switches.rs      # new switch BAANDI_SPIN_SHARE=on|off (default OFF)
  - crates/dft-numerics/tests/switches.rs    # default assertion + TABLE.len() +1
  - crates/dft-cli/tests/run.rs              # switches.len() +1
  - docs/performance/attempts-ledger.md      # append a row
READ_ONLY_CONTEXT:
  - AGENTS.md
  - crates/dft-engine/AGENTS.md
  - docs/design/onepath3-spin-metals.md
  - docs/decisions/0043-single-standard-path.md
  - docs/performance/fast-loop-protocol.md
  - crates/dft-engine/src/collinear.rs (CollinearProblem::map, the per-channel solve, atomic start sharing added in ONEPATH-3B)
INPUT_CONTRACT: |
  Existing collinear inputs. With BAANDI_SPIN_SHARE=off (default) behaviour is bit-identical to main. With on, the sharing applies only when ALL hold: fixed sector with n_up == n_down, both channels have the same band count and the same occupations, the two channel potentials (V_eff,up and V_eff,down) are bit-identical at this map (this holds when the spin density is exactly zero magnetization; test it by comparing the potentials, not by assuming), no spin-dependent external field/constraint, standard path. Otherwise fall back and record a reason in results.json (execution/spin block) and in a counter.
PHYSICS_CONVENTIONS: docs/physics-conventions.md (Hartree atomic units; spin-unpolarized limit: rho_up = rho_down = rho/2)
IMPLEMENTATION_STEPS:
  1. Read collinear.rs map(); locate where each spin channel's Davidson solve and density accumulation happen.
  2. Add the switch (default off) and the degeneracy test (bitwise potential comparison + occupation/band-count equality).
  3. When shared: run the solve for channel up only, copy the orbitals/eigenvalues/residuals to channel down (same data), build both channel densities from them, keep every downstream quantity (energies, magnetization, estimator terms E_H[d_n], E_H[d_m], E_0[d_m], residual guard, verification map) computed exactly as before from the two (identical) channel states. The verification map uses the same rule.
  4. Counters: number of shared maps, number of fallbacks with reason; expose in results.json (spin block) count-only.
  5. Tests (tests/spin_share.rs, bounds fixed BEFORE running and never weakened): (a) zero-moment spin run with SPIN_SHARE on equals the same run with it off within |dE| <= 1e-12 Ha, identical SCF/Davidson counts, per-channel eigenvalues equal; (b) the shared run issues exactly half of the Hamiltonian applications of the off run (count via the existing work counters) on a synthetic fixture (see crates/dft-engine/tests/collinear.rs fixtures); (c) fallback: non-zero seed, n_up != n_down, finite-temperature fractional occupations each take the two-channel path (counter and reason); (d) results with the switch off are bit-identical to main (compare against a fixture energy bit pattern recorded from main BEFORE your change).
TEST_COMMANDS:
  - scripts/check.sh dft-engine
  - cargo test -p dft-engine --test spin_share --offline
  - cargo test -p dft-numerics --test switches --offline
ACCEPTANCE_CRITERIA:
  - Switch off: bit-identical to main (fixture energy bits).
  - Switch on, zero-moment fixed sector: |dE| <= 1e-12 Ha vs off, same SCF and Davidson counts, H applications halved (+- bookkeeping), ranks/threads independent (T=1 and T=4 equal bits).
  - All fallback cases recorded.
  - scripts/check.sh all four gates PASS (the final full check is run by the orchestrator after review).
NON_GOALS: sharing for fractional occupations/finite temperature, k-point symmetry reduction, MPI paths, changing any tolerance.
DEPENDENCIES: none (main 23c938c or newer)
GIT: the worktree sandbox cannot commit; leave changes uncommitted and say so in the report; Claude reviews and commits.
RETURN_FORMAT: docs/schemas/report.schema.json
