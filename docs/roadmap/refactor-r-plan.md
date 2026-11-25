# Stage R plan: make the reference driver the only production SCF path

Status: draft by a Claude planning agent (2026-11-25), adopted by the orchestrator as the plan for stage R. Decisions D1-D6 (section 6) are pending with the user. ADR 0046 is taken by G4E (XC small/negative density), so the single-path ADR is 0047.

## 0. What I found that changes the plan

These facts were checked in the code. Each one limits what the tickets can do or the order they go in.

1. **The reference driver is not independent of legacy code at module level.**
   - `crates/dft-engine/src/scf.rs:3-4` declares `#[path = "refscf/mod.rs"] pub mod refscf;`, so `refscf` is a submodule of `scf`.
   - `refscf/driver.rs` and `refscf/k_driver.rs` do `use crate::scf::*`.
   - The reference path uses these legacy items:
     - `ScfOptions`, `ScfResult`, `ScfOutcome`, `PhaseTimings`, `AtomicStartReport`, `Verification*` (from `scf.rs` and `convergence.rs`)
     - `assess_fixed_point` (`session.rs:3132`)
     - `cutoff_ladder::prolong_{indices,orbitals,density,frequency}` (`cutoff_ladder.rs:1050,1099`)
     - `gamma_real::{GammaReal, LazyOrbitals, set_delivered_orbitals, lazy_result_allowed}`
     - `collinear::SpinSector`, `finite_temperature::{ThermodynamicLedger, validate_options}`
     - `scf::pool_fixture` (tests)
   - So shared types must be extracted before any legacy file can move or be deleted.
2. **The strict oracle is the legacy driver.**
   - `BAANDI_REFERENCE_STRICT=on` runs `ScfProblem` with `standard = None` (`config.rs:455`, `ScfOptions::strict_reference`).
   - `ScfProblem::run_impl` runs through `ScfSession` (`scf.rs:4202`).
   - So the oracle needs at least `scf.rs`, `session.rs`, `convergence.rs`, `collinear.rs`, `finite_temperature.rs`, `band_density.rs`, `atomic_start.rs`, `events.rs`, and the legacy solvers. If it must keep running at MPI8 on dev08, it also needs `gamma_real.rs`, `m2.rs` and `mix_dist.rs`.
3. **Precedence bug that would make oracle checks pass falsely.**
   - Today `run::options()` (run.rs:123) and `execute_mode_inner` (run.rs:1502) check `reference_requested` first.
   - With `REFERENCE_SCF=on` and `REFERENCE_STRICT=on`, the reference driver runs and `STRICT` is silently ignored.
   - After the default flips, every "oracle" run would actually be a reference run, and `E_ref - E_oracle = 0` would pass falsely.
   - R1 must make STRICT win, and every gate script must check the recorded driver.
4. **Config validation reads the raw switch.**
   - `config.rs:854` and `config.rs:1227` call `switches::on(REFERENCE_SCF)` and reject `mixing = "linear" | "kerker-pulay"` and `mixing_q0_inverse_bohr`.
   - If only the default flips, those legacy inputs fail at parse time instead of falling back.
   - Validation must use the routing decision instead.
5. **The README validation block (generated) says more than the "Works today" table.** `scripts/make_readme_figures.py:645-665` prints checkpoint/restart bit-identity and "paused runtime/audit restoration and cross-process CLI restart also pass", reading `tests/reference/checkpoint/results.json` and `tests/reference/runtime/generated/results.json`. These are `ScfSession` features. They cannot be removed while that text stands, unless the evidence is reproduced on the reference path or the user approves a generator change.
6. **The bench baseline is a legacy baseline.** `docs/performance/baseline.json` was recorded on the default path, which is legacy today. After the flip, `scripts/bench_scf.py --require-identical` against it is meaningless. A reference baseline must be recorded before the flip (R0).
7. **The reference start rejects some inputs** (`refstart.rs:398-422`, `driver.rs:244`, `k_driver.rs:109,135`, `electronic.rs:293-318`):
   - `strain_reference` (variable cell)
   - `pw92_prefactor != Modified` (affects `compat.mode = "reference"`, ADR 0038)
   - UPFs without `PP_RHOATOM`, without `PP_PSWFC`, or with l>3 orbitals (likely many synthetic fixtures)
   - density floors, uniform start
   - RMM or ladder on k lists, polarized LDA-PZ (legacy rejects this too)

   Today these only fail when the reference path is opted into. After the flip they need declared fallbacks.

---

## 1. Callers of the legacy SCF entry points, and what replaces each

Abbreviations:
- **Legacy entry points**: `ScfProblem::{new, with_kpoints}` + `run`, `run_outcome` or `step_kpoints`; `ScfSession` and `Runtime`; `cutoff_ladder::run_ladder`; `collinear::CollinearProblem`; `finite_temperature::FiniteTProblem`; `dft_engine::nscf`.
- **Reference entry points**: `refscf::driver::run[_with]`, `k_driver::run[_prepared]`, `electronic::run[_prepared]`, `delivered::{run, force_state[_prepared], ionic_seed, gamma_seed}`.

### 1a. CLI (`crates/dft-cli/src`)

| Caller | Legacy call | Reference replacement | Missing / action |
|---|---|---|---|
| `run.rs` `options()` :123 | `options_with_reference` → `ScfOptions` for every command | `reference_options` :246 | Routing must be centralised (R1); STRICT must take precedence |
| `run.rs` `execute_mode_inner` :1757-1993 | `ScfProblem::new` / `with_kpoints` + `enable_prepared_symmetry`, `set_density_floor`, `run_ladder` or `run_outcome` | `k_driver::run` / `run_prepared` :1886 and `delivered::force_state[_prepared]` :2046 (already wired) | Becomes the fallback and oracle branch only |
| `run.rs` :1720 → `finite_run::execute` | `CollinearProblem::new_fermi_dirac`, `FiniteTProblem::run` (finite_run.rs:53,70,146) | `reference_electronic_run` (run.rs:3004) → `electronic::run[_prepared]` | G4D free-moment PBE (no port; research). Moment FD stays legacy unless all three G4B switches are on (`reference_spin_fd_fallback_reason` :149) |
| `run.rs` :1745 → `spin_run::execute` | `CollinearProblem::run` (spin_run.rs:105) | `reference_electronic_run` | Nothing for fixed sectors (G4C qualified); used as fallback and oracle only |
| `telemetry.rs` `drive` :180 (`run --observe/--agent`) | `ScfSession` + `Runtime` (RetainLast) | none | The reference driver is monolithic (run.rs:1558 rejects it). Needs the reference stepper (R8a/R8c) |
| `control_channel.rs` :210 (test module only) | `ScfProblem` + `ScfSession` + `Runtime` | stepper | Moves with telemetry |
| `heal.rs` :130-183 (self-healing demo) | `ScfSession` | stepper | README validation block cites runtime evidence (finding 5) |
| `checkpoint_io.rs` `create` / `restart` :1853-2044 | `ScfProblem`, `ScfSession::{new, warm_start, restore_exact}`, Pulay/RMM state in `dft_engine::checkpoint` | none | Needs checkpoint v2 for Broyden history, density, orbitals and the ethr ratchet (R8b) |
| `checkpoint_inspect.rs` | reads checkpoint v1 files | – | Keep; extend for v2 |
| `relax.rs` `evaluate_with` :398-529 | `ScfSession` (Strict completion) with `initialize_density_transfer` / `initialize_cell_density_transfer`, C3 predictor | `delivered::run` + `force_state` + `ionic_seed` (:421-643, ions mode) | Variable-cell mode (`relax.rs:761` rejects it; refstart rejects `strain_reference`) and density extrapolation (`:767`) |
| `relax/cell.rs` :492 | `run::options` for grid / `strain_reference` | – | Ported in R7 |
| `bands.rs` :96-165 and `dos.rs` (uses `bands::Context`) | `ScfProblem::run` then `dft_engine::nscf` (nscf.rs:56, which calls `ScfProblem::with_kpoints` + `step_kpoints`) | SCF: `k_driver::run` | NSCF at fixed density on a k list does not exist on the reference path (R4) |
| `refscf_map.rs` :205-305 (diagnostic `refscf-map`) | `ScfProblem` for `Solver::Legacy` | its own `Solver::Reference` arm | Legacy arm moves under the oracle feature or is dropped |
| `refscf_check.rs` | none (reference only) | – | – |
| `inspect.rs` :181 (`explain`) | `run::options` | routing (R1) | Report the selected driver and reason |
| `config.rs` :551,643 | `collinear::SpinSector` type | – | Type moves to a shared module (R9a) |
| `examples/davidson_bench.rs` :627 | `ScfProblem` | – | Bench of legacy Davidsons: move under the oracle feature or delete in R10 |

### 1b. Python and MCP

| Caller | Legacy call | Replacement | Missing |
|---|---|---|---|
| `baandi-py/src/lib.rs` `Session::new` :216, `initialize` / `step` / `run` / `status` / `iteration` :253-286 | `ScfProblem::new` + `ScfSession` (in a worker thread) | `run()` → engine facade (R3) | `step()`, `status()`, `iteration()` need the stepper |
| `baandi-py` `forces` / `stress` :329,382 | `forces::forces` / `stress::stress` with `ElectronicState::Fixed(&outcome.result)` | `delivered::force_state` (same kernels) | none |
| `baandi-mcp/src/lib.rs` `register_input(ScfProblem)` :74, tools `create_calculation`, `run`, `step`, `pause`, `resume`, `cancel`, `propose_policy_change` :147-219 | `ScfSession` + `Runtime` + `policy::PolicyChange` | stepper + runtime adapter | Policy changes at safe points need a reference equivalent. Many legacy policy knobs (mixing beta, locking and similar) have no reference counterpart; the policy surface must be redefined |
| `baandi-mcp/tests/server.rs`, `examples/tool_schemas.rs` | same | – | Follow the MCP port |

### 1c. Engine tests and benches (legacy-using files; counts are `#[test]`/`#[ignore]` from grep)

| Group | Files (tests/ignored) | Disposition |
|---|---|---|
| Physics oracles on the strict path | `all/scf.rs` 8/4, `kpoints` 5/2, `nlcc` 6/2, `pbe` 2/1, `collinear` 6/1, `forces` 7/3, `stress` 6/6, `symmetry` 7/4, `symmetry_full_group` 2/1, `real_nc` 2/2, `electrostatic` 1/1, `force_identity` 2/1, `frozen_cell_metric` 1/0, `threading` 1/0, `spin_atomic_start` 4/1, `nscf` 3/1, `spin_share.rs` 5/0, `distributed_scf.rs` 11/0 | Keep, gated by `required-features = ["strict-oracle"]` (R9c). Physics oracle tests are never deleted |
| Legacy standard-path experiments | `cutoff_ladder` 17/2, `direct_min` 5/2, `scf_advance` 1, `standard_{stall,robust}_{contract,trajectory}` 8/2, `verify_schedule_alloc` 1, `first_block_diag` 2, `snapshot_identity` 3, `fp32_snapshots` 1/1, `instrument` 2, `trajectory_divergence` 3/2, `mixing_survey` 2/1, `preconditioner_{contract,study}` 3/1, `capacity_{contract,trajectory}` 4/1, `band_*_scan` / `_spectrum` 6/4, `atomic_start_{contract,option,prototype}` 9/2, `density_transfer` 3/1, `hybrid_budget` 1/1, `davidson_profile` 1, `stage_timers` 1, `init_scale` 1, `other_init` 3, `run_config` 1, `diagnostics` 9/3, `gamma_real_*` 15/8, `gamma_workspace_snapshot` 1/1 | Delete with their code in R10. Evidence-generator tests that read `tests/reference` (the `band_*`, `atomic_start_*`, `mixing_survey`, `trajectory_divergence`, `standard_*_trajectory` files) are deleted, but their stored JSON stays |
| Legacy MPI M2 | `m2_*.rs` (12 files, about 24 tests), `owned_handoff` 5/1, `init_owned_states` 6 | Kept with the full oracle. Deleted only if the oracle is trimmed (decision D2) |
| Runtime and checkpoint | `all/runtime_policy` 18/1, `all/checkpoint` 6/2, `all/stepwise_equivalence` 7/2, `common/checkpoint.rs` | Retarget to the stepper (R8) or keep while the legacy runtime is kept |
| Reference tests | `refscf.rs` 100/24, `refscf_forces.rs` 8/2, `refstart.rs` 6/2, `reference_driver.rs` 2/2 | Production tests. They use legacy only as an in-test oracle (37 hits in refscf.rs), so they need the feature too |

### 1d. CLI tests (`crates/dft-cli/tests`)

- **Run the binary on the default path; results change at the flip:** `all/run.rs` 21/4, `all/single_path.rs`, `all/onepath_*` (15/8), `all/volumetric.rs`, `all/stress.rs`, `all/relaxation.rs`, `all/self_healing.rs`, `spin.rs`, `gamma_real_config.rs`, `pseudo_convention_config.rs`, `pw92_prefactor_config.rs`, `preflight` tests. Each is pinned or re-baselined in R2.
- **Legacy-only:** `checkpoint_io.rs` 8/4, `checkpoint_inspect.rs`, `atomic_start_checkpoint.rs`, `reference_mode_restart.rs` (ADR 0038 compat mode), `all/cell_relaxation.rs` 7/7.
- **Reference:** `reference_driver.rs`, `reference_gaps.rs`, `refscf_map.rs`, `stop_to_verified.rs`.

---

## 2. Features only available on the legacy path

| # | Feature (legacy home) | User-visible? | Plan | Stage |
|---|---|---|---|---|
| L1 | Free-moment collinear FD with PBE (G4D) | yes | Keep the legacy fallback with reason `G4D_FALLBACK` until G4D qualifies (a research ticket outside R) | fallback |
| L2 | Moment FD without the G4B switches | yes | Turn the three G4B switches on by default, with dev08 gen6 evidence re-run at default; then delete the predicate | R6 |
| L3 | NSCF bands/DOS at fixed density (`nscf.rs`, `step_kpoints`) | yes (Works today) | Port: `refscf::nscf` | R4 |
| L4 | Variable-cell relax (`strain_reference`, frozen Miller basis; `initialize_cell_density_transfer`) | yes (Works today, README cell table) | Port: refstart frozen-strain support plus relax wiring | R7a/R7b |
| L5 | Relax density extrapolation C3 (`RELAX_EXTRAPOLATION` default off) | opt-in | Decision D5: port into `ionic_seed`, or retire with a migration error | R7c |
| L6 | Checkpoint/restart v1 (`checkpoint.rs`, `checkpoint_io.rs`) | yes (README validation block) | Port: stepper plus checkpoint v2 (D1) | R8b |
| L7 | observe/agent runtime, control socket, events, policy changes, self-healing (`session.rs`, `runtime.rs`, `policy.rs`, `events.rs`, `telemetry.rs`, `heal.rs`) | yes (README validation block) | Port through the stepper (D1) | R8a/R8c |
| L8 | Python `Session.step/status/iteration`, MCP tools | yes (APIs) | `run()` ported in R5; step and MCP in R8d | R5, R8d |
| L9 | Legacy-only SCF input keys: `mixing = linear / kerker-pulay`, `mixing_q0_inverse_bohr`, `hartree_mixing`, `spectral_guard`, `locking`, `preconditioner`, `davidson_capacity_per_band`, `start_density = uniform`, `start_wavefunctions` (random), `energy_tolerance_ha`, `density_residual_tolerance`, `eigen_tolerance` (strict only), `standard_stall_guard`, `standard_density_forcing`, legacy `cutoff_ladder`, `xc.density_floor` | yes | During the transition: route to legacy with a reason (R1). In R10c: input error with a migration message such as "applies only to the strict oracle (BAANDI_REFERENCE_STRICT=on)" (D4) | R1, R10c |
| L10 | `compat.mode = "reference"` (ADR 0038: PW92 original prefactor, QE local mesh) | yes | refstart rejects `pw92 != Modified`. Port: pass the prefactor, local mesh and radial rule through refstart. These are physics constants, so the cost is small | R6b |
| L12 | Legacy cutoff ladder, direct minimisation, SCF_ADVANCE, fixed map, spectral SCF, first block, replay capture | experiment | Delete. The reference path has its own `refscf/ladder.rs` and RMM, both opt-in | R10a |
| L13 | Legacy standard-path stop machinery (stall guard, density forcing, stop predictor, verify schedule, cert screen, stop metric) | default on the legacy standard path | Delete once the fallback list is empty; the oracle keeps only the strict path | R10b |
| L14 | Legacy MPI M2 distributed Gamma, `gamma_real` complex↔real bridge, lean/distributed Davidsons | oracle at MPI | Decision D2 | R10d |
| L15 | `refscf-map --solver legacy` | diagnostic | Under the oracle feature | R9c |

L11 (optional atomic charge/orbitals): **done in R6c**. The reference start uses
neutral uniform density and deterministic kinetic completion when needed; it no
longer routes these UPFs to legacy.

Reference-only features need nothing: RMM, the reference ladder, kept HX, and the verify gate are already reference.

---

## 3. Proposed final module layout

```
crates/dft-engine/src/
  lib.rs                 re-exports; `pub mod refscf` at top level (no longer nested in scf.rs)
  model/                 NEW (R9a): pure data and contracts shared by production and oracle
    options.rs           ScfOptions, Mixing, CompatMode, GammaReal, StartDensity/Wavefunctions, OccupationPolicy
    result.rs            ScfResult, ScfOutcome, KPointSolution, IterationRecord, PhaseTimings, FixedPointCertificate
    verify.rs            Verification*, StandardProfile, assess_fixed_point (from session.rs:3132)
    spin.rs              SpinSector, ThermodynamicLedger, finite-T validate_options
    error.rs             ScfError
  refstart.rs            (2075) production start
  refscf/                (10112 incl. refstart) production SCF: driver, k_driver, electronic, delivered,
    prolong.rs           NEW: prolong_* moved from cutoff_ladder.rs (~150 lines)
    lazy.rs              NEW: LazyOrbitals, set_delivered_orbitals from gamma_real.rs (~350 lines)
    nscf.rs              NEW (R4)
    stepper.rs           NEW (R8a, if D1 = port)
    api.rs               NEW (R3): one facade for CLI / Python / MCP
  energy.rs forces.rs stress.rs force_identity.rs functional.rs grid_dist.rs parallel.rs
  profile.rs phase_clock.rs        shared (6386 with nscf/lib)
  oracle/                NEW (R9c), `#[cfg(feature = "strict-oracle")]`, default feature ON
    scf.rs session.rs convergence.rs collinear.rs finite_temperature.rs band_density.rs
    atomic_start.rs start_density.rs spectral.rs diagnostics.rs nscf.rs events.rs
    gamma_real.rs m2.rs mix_dist.rs          (kept unless D2 = trim)
```

**Where the oracle lives: a Cargo feature module, not a separate crate.**
- A separate `dft-oracle` crate needs the frozen root `Cargo.toml` and `ownership.toml` to change.
- It would also need dozens of engine-private items (band density, energy audit internals, the gamma-real bridge) made public.
- A feature-gated `oracle/` module gives the same isolation. `cargo check -p dft-engine --no-default-features` must compile, and that build contains no oracle code; it becomes a new `scripts/check.sh` gate.
- The oracle stays reachable from the production binary only through `BAANDI_REFERENCE_STRICT=on`. All dev08 gate scripts already use that switch.
- `dft-cli` enables the feature by default. A later `--no-default-features` production build is possible once the fallback list is empty.

**Line counts (`wc -l`, today):**

| Group | Files | Lines |
|---|---|---|
| dft-engine src total | | 46,066 |
| dft-engine tests total | | 36,151 |
| Reference (refscf + refstart) | | 10,112 |
| Oracle core | scf 4868, session 3479, convergence 734, collinear 1675, finite_temperature 551 + tests 376, band_density 807, atomic_start 1196, start_density 170, gamma_real 3401, m2 595, mix_dist 1379, verification_tests 983, spectral 172, diagnostics 473 | 20,859 |
| Legacy experiments (delete) | cutoff_ladder 3395, direct_min 873, advance 344, fixed_map 265, spectral_scf 68, replay_capture 56, first_block 361, verification_profile 237, floor_tests 87 | 5,686 |
| Runtime (replace) | runtime 1169, policy 329, events 209, checkpoint 386, snapshot 906, state_view 24 | 3,023 |
| CLI legacy-only | finite_run 538, spin_run 315, telemetry 288, heal 297, control_channel 299, checkpoint_io 2252, checkpoint_inspect 128 | 4,117 |
| Experiment tests (engine) | | 5,830 |
| M2 tests | m2_*, owned_handoff, init_owned_states, distributed_scf | 4,369 |
| Runtime tests | | 2,205 |
| Legacy-only solvers (dft-solvers) | real_davidson 3587, lean_davidson 2448, davidson_stages 859, distributed_davidson 2292 (+ subdirectories with mixing and ref_davidson, 6935 total), filtered_subspace 782, bordered 719, rmm 649, locking 452, mixing 1345, legacy 131, spectral 327 | 13,591 + part of 6,935 |
| `dft-physics/src/real_gamma.rs` (legacy only) | | 1,012 |
| Switch table `dft-numerics/src/switches.rs` | 186 switches | 2,194 |

**Estimated removal:**

| Scope | What | Lines removed |
|---|---|---|
| Certain (R10a + R10b + R10c + R11) | Experiments src ~5.5k (minus ~150 moved prolong lines); experiment tests ~5.8k; legacy standard-path branches in `session.rs` / `scf.rs` ~2-3k; CLI legacy branches ~1.2k (run.rs ladder/legacy failure path ~700, relax legacy branch ~130, bands ~60, refscf_map legacy arm ~150, all estimates); switch entries and dead arms ~0.5-1k | ≈ 14-16k |
| With D1 = port | Runtime and checkpoint are rewritten, not deleted; legacy session-runtime glue | net ≈ -2k |
| With D2 = trim the oracle to serial/threaded complex + real-Gamma | m2, mix_dist, most of gamma_real, distributed/lean Davidsons, davidson_stages, real_gamma.rs, M2 tests | ≈ +18-20k more |

Final dft-engine plus dft-solvers is roughly 35-40% smaller with D2, or about 15% without it.

---

## 4. Ordered stages (each stage is one ticket)

### Gate definitions (referenced by every stage)

**Local (Mac):**
- `L-CHECK`: `scripts/check.sh` exits 0. Record passed, failed and ignored counts plus the binary count, and compare with the previous merge.
- `L-README`: `python3 -I scripts/make_readme_figures.py` leaves no diff, and `python3 -I -m unittest scripts/tests/test_readme_figures.py` passes.
- `L-ID`: release `scripts/bench_scf.py --require-identical --baseline docs/performance/baseline-reference.json` (cases gamma and k222s) passes, plus the new debug golden test `reference_identity` (density sha256 and energy bits, keyed by target triple).
- `L-ORC`: the same with `BAANDI_REFERENCE_STRICT=on` against `baseline-oracle.json`.
- `L-ROUTE`: the routing-matrix test (R1) passes.

**dev08 (through the master session):**
- `RG-ID`, hot-path identity: si32, si64 and si128 Gamma LDA at the default settings, MPI8 pinned and one rank with 8 threads. Must match the R0 capture: density sha256, energy hex bits, map count, H columns. Wall median of 3 crossed A/B repeats must be within the A/A floor.
- `RG-ORC`, oracle identity: `REFERENCE_STRICT=on` on si8-gamma, Si2 k222s, Si2 FD θ=0.01 k444, the G4C spin fixture and the G5 si8 displaced forces/stress case. Energy bits, density hash and force bits must equal R0.
- `RG-ACC`, accuracy matrix (a "gen7" rerun of the roadmap evidence table at the production default, no switches set):
  - Cases: Mg NLCC; si32 PBE; MgO PBE NLCC; Si2 k444 unshifted; Si2 k222 shifted; Si2 FD θ=0.01 k444; G4C fixed-sector spin PBE; G6 symmetry case; si8 displaced forces and stress; 3-step ion relax.
  - Gates: |E − E_oracle| ≤ 5e-9 Ha (standard class); label `verified`; forces ≤ 1e-5 Ha/bohr (G5 saw 8.2e-6); stress ≤ 1e-8 Ha/bohr³; |μ − μ_oracle| ≤ 1e-6 Ha; relax at threshold 5e-13 within dE 5e-9, dF 1e-6, dr 1e-6.
  - `results.json` must show `scf.driver = "reference"`, and oracle runs must show `"strict_oracle"`. A gate script that finds the wrong driver fails.
- `RG-QE`: the six README QE cases (Si 30 Ry: LDA Gamma, k222, k222s, k444s; PBE Gamma, k222s) on the default binary. They must be within the README tolerances (1e-4 Ha/atom, 1e-3 Ha). Report the change against the stored differences; this is for information and decides decision D6.
- `RG-ROUTE`: the routing corpus run on dev08 with MPI8, checking the recorded driver and reason.

A stage marked "no numerical change" needs `RG-ID` and `RG-ORC` to be bit-identical. A stage that ports a feature adds that feature's own oracle comparison.

---

### R0: identity baselines and routing corpus (no production code change)

- **Goal:** freeze what must not change before anything moves.
- **OWNED_FILES:**
  - `docs/performance/baseline-reference.json` (new; `bench_scf.py --record` with `BAANDI_REFERENCE_SCF=on`)
  - `docs/performance/baseline-oracle.json` (new; with `BAANDI_REFERENCE_STRICT=on`)
  - `crates/dft-cli/tests/reference_identity.rs` (new): in-process tiny Si2 Gamma LDA (synthetic UPF from `scripts/gen_synthetic_upf.py` with `PP_RHOATOM`/`PP_PSWFC`), under 10 s in debug. It runs both paths and asserts density sha256 and energy bits per `target_os`/`target_arch`, skipping with a note on other targets.
  - `crates/dft-cli/tests/fixtures/routing/*.toml` (new): one input per L-row of section 2, each with an `expected_driver` and `expected_reason` comment.
  - `scripts/bench_scf.py` (records `scf.driver`; refuses a baseline whose driver differs)
- **Tests:** `L-CHECK` (count +2 tests), `L-README`.
- **dev08:** capture the R0 reference for `RG-ID`, `RG-ORC` and `RG-ACC`, stored as `/data/sig-dft/r0/RESULT.md`, with the commit and binary sha.

### R1: one driver selector, with fallback reasons recorded everywhere (default still legacy)

- **Goal:** add `driver_select::select(config, mode, command, upfs) -> Selection { driver: Reference | LegacyStandard | StrictOracle, reason: Option<&'static str>, source: Default | Env | Input }`.
  - Precedence: `REFERENCE_STRICT=on` → StrictOracle (this fixes finding 3); then explicit `scf.driver`; then env `REFERENCE_SCF`; then the default.
  - The legacy-only predicates (L1, L2, L4-L11) are checked after the UPF parse and before engine construction. They are declared predicates; nothing is caught and retried after a refstart error, because that would hide bugs.
  - Every command (`run`, `relax`, `bands`, `dos`, `checkpoint`, `restart`, `explain`, `refscf-map`) and `config.rs` mixing validation (finding 4) use the selection instead of `switches::on(REFERENCE_SCF)`.
  - `results.json` gains `scf.driver`, `scf.driver_reason` and `scf.driver_source` for every command, generalizing today's `reference_driver_fallback_reason`. The routing record is written to stderr for the user.
- **OWNED_FILES:** `crates/dft-cli/src/driver_select.rs` (new), `crates/dft-cli/src/lib.rs`, `run.rs` (`options`, `reference_requested`, `reference_spin_fd_fallback*`, the `execute_mode_inner` header), `config.rs` (:854, :1227, :1312), `relax.rs` (:398, :760), `bands.rs`, `checkpoint_io.rs` (option resolution only), `inspect.rs`, `refscf_map.rs`, `crates/dft-cli/tests/driver_routing.rs` (new, drives the R0 corpus).
- **Tests:**
  - `L-ROUTE` covers all corpus inputs under default off, env on, input reference, input legacy, and STRICT combined with each.
  - `L-ID` and `L-ORC` stay bit-identical.
  - `L-CHECK` count is the previous count plus the routing tests; no existing test changes except added JSON fields.
- **dev08:** `RG-ROUTE`, and `RG-ORC` with `REFERENCE_SCF=on` + `STRICT=on`, which must now report `strict_oracle`.

### R2: flip the CLI default to the reference driver

- **Goal:** set the `REFERENCE_SCF` default to `"on"` and rewrite its meaning text ("production driver; off = legacy standard path, deprecated").
  - Any existing CLI test whose expected output is legacy-specific is either pinned with `[scf] driver = "legacy"` (when it tests a legacy feature) or re-baselined (when it tests generic behavior). The ticket lists every test in each category.
  - `bench_scf.py` defaults to `baseline-reference.json`, and `baseline.json` is renamed `baseline-legacy.json` with a note.
- **OWNED_FILES:** `crates/dft-numerics/src/switches.rs` (one entry), `crates/dft-numerics/tests/all/switches.rs` (`legacy_defaults_are_unchanged`), `crates/dft-cli/tests/**` (pins and re-baselines), `scripts/bench_scf.py`, `docs/performance/baseline*.json` (rename), `docs/performance/README.md` (baseline paragraph).
- **Tests:**
  - `L-ID`: the default run is now bit-identical to R0's `REFERENCE_SCF=on` baseline, which proves the flip changes only the selection.
  - `L-ORC` is unchanged.
  - `L-ROUTE`.
  - `L-CHECK`: the test count stays the same; the list of modified expectations goes in the report.
- **dev08:** full `RG-ACC` at default (no switches), `RG-ID` (default now equals R0 with the switch on), `RG-ROUTE`, `RG-QE`.

### R3: engine facade

- **Goal:** add `refscf::api::solve(StartInput, Request) -> Delivered`.
  - `Request` holds: profile, points, `Option<symmetry::Prepared>`, `Option<electronic::Policy>`, `Option<IonicSeed>`, and the derivatives wanted (forces, stress).
  - It dispatches to `driver`, `k_driver` or `electronic` and `delivered` exactly as `run.rs:1886-2050` and `relax.rs:421-443` do today.
  - The CLI `run` and `relax` reference branches call it. It is the single entry point later used by Python, MCP and bands.
- **OWNED_FILES:** `crates/dft-engine/src/refscf/api.rs` (new), `refscf/mod.rs`, `crates/dft-cli/src/run.rs` (reference branches only), `relax.rs` (reference branch).
- **Tests:** `L-ID`; `results.json` byte-identical except wall fields (add a JSON-diff helper test with timing keys masked); `L-CHECK` unchanged.
- **dev08:** `RG-ID`.

### R4: reference NSCF (bands and DOS) — two tickets

- **R4a (engine):** `refscf::nscf(input, density, points, tolerance)`.
  - Builds the complex k start from a supplied density. This needs a new `refstart::start_kpoints_from_density` sharing the atomic-orbital start vectors.
  - Solves all bands, plus at least 4 empty bands, with `ref_davidson_complex` to a residual of `eigen_tolerance` for every band, with no mixing.
  - The density grid is fixed, matching legacy `nscf.rs` semantics (charge check, at least four empty bands).
  - OWNED: `crates/dft-engine/src/refscf/nscf.rs` (new), `refstart.rs` (new entry only), `refscf/mod.rs`, `crates/dft-engine/tests/refscf_nscf.rs` (new).
  - Test: on the same input density, reference NSCF eigenvalues match legacy `nscf` within 1e-8 Ha (both residual-certified at 1e-8 Ha) on Si2 k222s, a skew two-species cell and MgO Gamma PBE NLCC.
- **R4b (CLI):** `bands.rs` and `dos.rs` use `api::solve` for the SCF and `refscf::nscf` for the path. The `driver_select` reason L3 is removed.
  - OWNED: `crates/dft-cli/src/bands.rs`, `dos.rs`, `driver_select.rs` (remove the L3 predicate), `crates/dft-cli/tests/all/bands_dos.rs`.
  - Test: the existing `bands_dos` format tests, plus an output comparison with the legacy pin within 1e-7 Ha (standard-class SCF difference propagated).
- **dev08:** Si2 band path, Si2 k444 DOS and MgO PBE bands, each against the oracle SCF density followed by legacy NSCF: max |Δε| ≤ 1e-6 Ha (this includes the standard-class density difference); `RG-ID`.

### R5: Python `run()`, forces and stress on the reference path

- **Goal:** `Session.run()`, `forces()`, `stress()` and `result()` go through `api::solve`.
  - Before the stepper exists (R8), `initialize()`, `step()`, `status()` and `iteration()` either raise `NotImplementedError("stepwise control requires driver='legacy' until R8")` or run on an explicit `Options.driver = "legacy"`.
  - Note: the `Options` mixing defaults change (reference beta 0.7, history 8). This must be documented in `VALIDATION.md`.
- **OWNED_FILES:** `crates/baandi-py/src/lib.rs`, `crates/baandi-py/VALIDATION.md`, `crates/baandi-py/AGENTS.md`.
- **Tests:** a native build test of energy, forces and stress against the CLI on the same input (bit-identical, same engine call); ASE evidence regenerated as a new file (old evidence kept).
- **dev08:** none beyond `RG-ID`.

### R6: remove the cheap fallbacks

- **R6a G4B default:** turn `REFSCF_SPIN_FD_MOMENT`, `REFSCF_SPIN_STIFFNESS` and `REFSCF_HISTORY_CONDITIONING` on by default; remove the L2 predicate.
  - OWNED: `switches.rs` (3 entries), `crates/dft-cli/src/driver_select.rs`, `crates/dft-cli/tests/spin.rs`, `crates/dft-numerics/tests/all/switches.rs`.
  - dev08: the gen6 G4 matrix at default (moment FD Si2, the G4C fixture) against the oracle, standard class; `RG-ID` (the Gamma LDA path does not read these switches).
- **R6b compat reference mode (L10):** refstart accepts `Pw92Prefactor::Original`, the QE local mesh and the radial rule (`refstart.rs:403`).
  - OWNED: `crates/dft-engine/src/refstart.rs` (validation and pass-through), `refscf/map.rs` or `electronic_map.rs` if the XC prefactor is threaded there, `driver_select.rs`, `crates/dft-cli/tests/reference_mode_restart.rs` and `all/reference_mode_*`.
  - Test: compat-mode Si2 reference vs oracle within 5e-9 Ha.
  - dev08: `RG-ID` (the default prefactor path must stay bit-identical).

### R7: variable-cell relaxation on the reference path

- **R7a refstart frozen strain:** honour `ScfOptions::strain_reference`. Plane waves are selected on the reference cell's Miller set and evaluated at the strained cell (frozen basis, as in legacy). This applies to the Gamma real and complex k starts, and to `delivered::stress_state` at the strained cell.
  - OWNED: `crates/dft-engine/src/refstart.rs`, `refscf/delivered.rs`, `crates/dft-engine/tests/refscf_forces.rs`.
  - Test: si2 strained cell, reference stress vs oracle within 1e-8 Ha/bohr³; the frozen-basis energy derivative checked by finite difference.
  - Size: the largest single ticket. If the implementer overruns, split it into Gamma and k parts.
- **R7b relax wiring:** `relax.rs` and `relax/cell.rs` cell modes use `api::solve` with the cell seed (density transfer between cells reuses `ionic_seed` with the frozen Miller map). Remove the `relax.rs:761` rejection.
  - OWNED: `crates/dft-cli/src/relax.rs`, `relax/cell.rs`, `driver_select.rs`, `crates/dft-cli/tests/all/cell_relaxation.rs`.
  - dev08: the stored cell-relax cases (`tests/reference/relaxation/cell/validation.json` protocol) against the oracle: final lattice ≤ 1e-5 bohr, pressure ≤ 1e-6 Ha/bohr³, dE ≤ 5e-9 Ha per step; `RG-ID`.
- **R7c extrapolation (D5):** either port the C3 predictor into seed construction or retire `RELAX_EXTRAPOLATION` and the `start_extrapolation` key with a migration error. Look up its row in the attempts ledger first.

### R8: stepper, checkpoint and runtime (only if D1 = port) — four tickets

- **R8a stepper:** split the loops of `refscf::driver`, `k_driver` and `electronic` into `Stepper { initialize(), step() -> StepReport, result() }`. `run()` becomes `Stepper::run_to_end()` with the same arithmetic order.
  - OWNED: `refscf/driver.rs`, `k_driver.rs`, `electronic.rs`, `refscf/stepper.rs` (new), `refscf/record.rs`.
  - Gate: this is the riskiest hot-path ticket. `L-ID`, `RG-ID` and a dev08 wall check (no regression beyond the A/A floor at si128 MPI8) are all mandatory.
- **R8b checkpoint v2:** serialize the density in planes, the Broyden history, kept orbitals, the ethr ratchet and the verification state. Restart must be bit-identical to the uninterrupted run.
  - OWNED: `crates/dft-engine/src/refscf/checkpoint.rs` (new), `crates/dft-cli/src/checkpoint_io.rs`, `checkpoint_inspect.rs`, `docs/checkpoint-format.md`, `crates/dft-cli/tests/checkpoint_io.rs`.
  - Gate: memory and file restart identity on Gamma, k4 and skew two-species fixtures (the same protocol as `tests/reference/checkpoint/results.json`); the new evidence is written to a new file.
- **R8c runtime:** `telemetry.rs`, `control_channel.rs` and `heal.rs` drive the stepper. Policy changes are limited to reference-meaningful ones (beta, history, threshold, verification, pause/resume/cancel). Self-healing evidence is regenerated as a new file.
- **R8d:** MCP tools and Python `step()` on the stepper.
- **README:** the generator gains the reference evidence files. The old checkpoint and runtime paragraphs are re-pointed only after the user approves (D6).

### R9: move the strict oracle (pure moves, bit-identical by construction)

- **R9a extract shared surface** into `crate::model` (section 3), keeping `pub use` at the old paths so no caller changes. Also move `prolong_*` into `refscf/prolong.rs` and `LazyOrbitals` into `refscf/lazy.rs`.
  - OWNED: `crates/dft-engine/src/{scf.rs, session.rs, convergence.rs, collinear.rs, finite_temperature.rs, cutoff_ladder.rs, gamma_real.rs, lib.rs}`, `model/*` (new), `refscf/{prolong.rs, lazy.rs, mod.rs, ladder.rs, driver.rs}`.
  - This can be pulled forward to right after R2 if OWNED sets allow.
- **R9b de-nest:** `refscf` and `atomic_start` become top-level modules in `lib.rs` (remove `scf.rs:1-4` `#[path]`). `refscf` imports `crate::model` instead of `crate::scf::*`.
- **R9c feature gate:**
  - `git mv` the legacy files into `src/oracle/`; `[features] default = ["strict-oracle"]` in `crates/dft-engine/Cargo.toml`; `#[cfg(feature = "strict-oracle")]` on the module and re-exports.
  - Tests that use the oracle get `required-features`.
  - `dft-cli` declares `dft-engine = { ..., features = ["strict-oracle"] }`.
  - New `check.sh` gate: `cargo check -p dft-engine --no-default-features` and `cargo clippy -p dft-engine --no-default-features -D warnings`.
  - OWNED: those files, `crates/dft-engine/Cargo.toml`, `crates/dft-cli/Cargo.toml`, `scripts/check.sh`, `crates/dft-engine/AGENTS.md`, test manifests.
  - Not allowed: the root `Cargo.toml` (frozen; not needed).
- **Tests for all of R9:** `L-ID`, `L-ORC` bit-identical, `L-CHECK` with identical counts (moves only).
- **dev08:** `RG-ID` and `RG-ORC`, both bit-identical.

### R10: delete dead code

- **R10a legacy experiments (L12):**
  - Delete `cutoff_ladder.rs` (except the already-moved prolong code), `direct_min.rs`, `advance.rs`, `fixed_map.rs`, `spectral_scf.rs`, `replay_capture.rs`, `first_block.rs`, `verification_profile.rs`, `floor_tests.rs`, and their tests (section 1c group 2: cutoff_ladder, direct_min, scf_advance, first_block_diag, snapshot/fp32/instrument where only those modules use them).
  - Delete the CLI ladder JSON in `run.rs`.
  - Remove the dependent switches: `CUTOFF_LADDER*`, `LADDER_*` (9), `SCF_ADVANCE*` (5), `FIRST_BLOCK`, `SNAPSHOT_DIR` if orphaned.
  - Removal: about 5.5k src lines plus about 4-5k test lines; about 30 integration tests. The ticket reports the exact before and after counts.
  - Tests: `L-ORC` (the oracle does not use these) and `L-ID`.
  - dev08: `RG-ID`, `RG-ORC`.
- **R10b legacy standard path (only once the fallback list is just G4D, or empty):**
  - Remove the standard-profile branches from `oracle/session.rs` and `oracle/scf.rs`: stall guard, density forcing, stop predictor, verify schedule, cert screen, stop metric, legacy verification attempts. The oracle becomes strict-only.
  - Remove the standard-mode `spin_run` and `finite_run` CLI routes, except a G4D fallback if D3 keeps it.
  - Delete the matching `standard_*` tests and `verify_schedule_alloc`.
  - Gate: `L-ORC` and `RG-ORC` bit-identical (strict arms untouched).
- **R10c CLI switches and keys:**
  - Retire the env switch `REFERENCE_SCF` and the input value `scf.driver = "legacy"` (both become migration errors).
  - Turn the L9 and L11 input keys into errors that name the oracle.
  - Reduce `driver_select` to Reference | StrictOracle.
  - Update `L-ROUTE`.
  - dev08: `RG-ROUTE`, `RG-ACC`.
- **R10d (only if D2 = trim):**
  - Delete oracle MPI M2 (`m2.rs`, `mix_dist.rs`, most of `gamma_real.rs`), `dft-physics/real_gamma.rs`, the legacy lean, distributed and staged Davidsons, and the M2 tests.
  - Precondition: the trimmed oracle agrees with the R0 oracle within 1e-12 Ha (decision of 2026-11-11 #1) on the generality set.
  - Large-cell oracle runs then use single-rank threads.

### R11: switch retirement — three tickets

**Criteria:**
- (a) **Orphaned**: read only by deleted code. Removed in the same ticket that deletes the code.
- (b) **Oracle-internal**: read only under `oracle/`. Freeze the oracle at its defaults by deleting the non-default arm and the switch, so that only `REFERENCE_STRICT` selects it. The default arm is untouched, so `RG-ORC` must be bit-identical. A tunable oracle is not an oracle.
- (c) **Reference experiment, default off**: its attempts-ledger status is rejected, superseded or gated, with no later premise-change row. It is not an isolated algorithm module protected by the Generality rule (keep `REFSCF_RMM*`, `REFSCF_LADDER*` if listed as opt-in, `MKL_*` and Kerker/TF preconditioning). It is not diagnostic. Delete the switch and the off-default code. `RG-ID` must be bit-identical, since the default arm is untouched.
- (d) **Default-on rollback switches on the reference path** (REF_DAV_*, REFSTART_*, REFSCF_KEPT_HX, NONLOCAL_*, FFT_*): keep them. AGENTS.md says "defaults stay switchable for A/B and rollback". Retiring them needs a sunset rule from the user (D3). Suggested rule: adopted at least two dev08 evidence generations ago, with no rollback use in any gate since.
- (e) **Always keep**: diagnostics and profiling (`PHASE_PROFILE`, `STAGE_TIMERS`, `DAVIDSON_STAGES`, `FFT_PROFILE`, `COLLECTIVE_TRACE`, `INIT_PROFILE`, `GLUE_PROFILE`, `PHASE_VALIDATE`, `PHASE_WAIT` if still read), infrastructure (`PSEUDO_DIR`, `BINARY_HASH_CACHE`, `DGEMM_TIMERS`, `MPI_TREE`), and `REFERENCE_STRICT`.

**Classification from grepping `Id::X` in each crate:**

| Class | Count | Examples |
|---|---|---|
| Read only in legacy engine files | ~63 | CERT_SCREEN, POTENTIAL_REUSE, EIGEN_FACTOR, EMPTY_LOOSE, BUFFER_TOL, LADDER_* (9), SCF_ADVANCE* (5), STANDARD_STOP, DENSITY_PAIR/REDUCE, ATOMIC_*, MPI_M2, M2_{BATCH, DENSITY, REUSE, LAZY, SCREEN, ATOMIC, GRID, MIX, LOCAL_*, OWNED_*, LIGHT_STATE, LAZY_SERIAL, LAZY_RESULT}, STOP_METRIC, STOP_PREDICTOR, VERIFY_SCHEDULE, VERIFY_REUSE, MIX_COULOMB_GUARD, MIX_GROWTH_RESET, MIX_CELL_KERKER, SPIN_SHARE, EWALD_CACHE, IMPORT_CHECK, VERIFY_PROFILE, FIRST_BLOCK, SNAPSHOT_DIR |
| Read only in legacy solver and physics code | ~37 | B1*, B5, B6, RITZ_FUSED, NARROW_*, SECOND_VECTOR, ORTHO, MGS, ORTHO_BLOCK, WIDE_DOT, PROJECT_WIDE, DAVIDSON, LEAN_*, DAVIDSON_{PACKED_REDUCE, ROTATE_COLUMNS, UPPER_BORDER, PROFILE}, RESTART_KEEP, M2_{ORTHO, LEAN, EIGEN, B5, WORKSPACE, SELECTED_ROOTS, PANEL_GS, START, LEAN_ALGEBRA, LAYOUT_CACHE}, MIX_INCREMENTAL, MKL_NONLOCAL, APPLY_LEAN, HAPPLY_{OVERLAP, FOLD} |
| Reference and shared kernels | ~70 | REF_DAV_* 14, REFSCF_* ~20, REFSTART_* ~12, FFT_* 8, NONLOCAL_* 4, numerics kernel 6, EWALD_FAST, IONIC_FAST, M2_FFT, M2_STICKS, M2_DGEMM (shared with `distributed_real_gamma`) |
| CLI and infrastructure | ~16 | – |

The legacy engine and solver groups (about 100 switches) are classes (a) and (b). The reference group goes through (c), (d) and (e). Default-off candidates for (c) include `REF_DAV_STACK_ROTATE`, `REF_DAV_SETUP_FUSED`, `REFSCF_KEPT_RECOMPUTE`, `REFSTART_LOCAL_SHELLS`, `REFSCF_VERIFY_REUSE`, `REFSCF_VERIFY_GATE_FIRST`, `REFSCF_PERMAP`, `FFT_SPARSE_CLEAR`, `FFT_DIRECT_PACK`, `FFT_A2A_BALANCE`, `FFT_KEPT_FUSE`, `HAPPLY_*` and `FORCE_BAND_FIRST` (in `forces.rs`, so shared). Each needs its ledger row cited; e.g. the ledger line 3735 says "STACK_ROTATE and SETUP_FUSED stay off".

**Tickets:**
- **R11a:** add a `scope: Reference | Oracle | Infra | Diagnostic` field to the `table!` macro, plus a test that an Oracle-scoped `Id` is read only under `oracle/`. Add `RETIRED: &[(&str, &str /*former default*/)]`: setting a retired variable to its former default prints a warning; any other value is an error, so no behavior changes silently. This protects dev08 scripts that still export `BAANDI_REFERENCE_SCF=on`, which is the default value. Update `RUN_COUNTERS`/labels and `legacy_defaults_are_unchanged`.
  - OWNED: `crates/dft-numerics/src/switches.rs`, `crates/dft-numerics/tests/all/switches.rs`, `scripts/bench_scf.py` (manifest check).
- **R11b:** oracle-internal retirement, (b). Gate `RG-ORC` bit-identical.
- **R11c:** reference default-off retirement, (c), citing the ledger row per switch. Gate `RG-ID` bit-identical. Add a ledger row "retired in R11c" per switch.

### R12: documentation

- New ADR `docs/decisions/0047-reference-driver-single-path.md`: the routing, the oracle feature, the fallback list and its sunset, the switch rules. It supersedes the sentence in `docs/design/reference-scf.md` "Rules" that says the old driver and switches are removed in one step.
- Also update:
  - `docs/roadmap/reference-generalization.md` (row R evidence)
  - `docs/design/reference-scf.md` (new section "R: production routing")
  - `docs/INDEX.md`
  - `crates/dft-engine/AGENTS.md`, `crates/dft-cli/AGENTS.md`
  - `docs/performance/README.md` (baselines)
  - `docs/checkpoint-format.md`
  - `docs/progress.md`, `artifacts/*.json`
- README: only the generator-owned block changes (the stored evidence is unchanged, so `L-README` passes). Factual corrections outside the block, to be reported to the user: line 94 "Pulay and linear density mixing" (the production mixer becomes Broyden/plain), line 101 "Buffer bands" (legacy `EMPTY_LOOSE`/`BUFFER_TOL`), line 104 "finite-T restart unavailable" (still true), and line 112 (cell relax) after R7.

**Dependencies:**
- R0 → R1 → R2 → R3.
- {R4, R5, R6, R7} run in parallel after R3 (their OWNED sets are disjoint except `driver_select.rs`; serialise edits to it).
- R8 after R3 and only if D1 = port.
- R9 after R2 (R9a may go earlier).
- R10a after R9c.
- R10b after R4, R6, R7 and R8 (or D1 = retire).
- R10c after R10b. R10d is optional.
- R11 after R10a (R11b after R9c).
- R12 last, with ADR text drafted in R2.

---

## 5. Risks

1. **Gamma-LDA hot path identity.**
   - The production hot path is `refscf/driver.rs`, `map.rs` and `kept_hx.rs`, plus `dft-solvers` `ref_davidson*` and `distributed_fft`, `dft-physics` `distributed_real_gamma`, and `dft-numerics` `real.rs`/`fused_pair_plan`.
   - R1-R7, R9 and R10 do not touch these files except R8a (stepper) and R11 (switch arms).
   - Moving code between modules cannot change floating-point results, but it can change inlining and code layout. So every move stage still needs `RG-ID` wall within the A/A floor at si128 MPI8, not only hash identity.
   - R8a is the highest risk: the stepper must keep the arithmetic order and must not add allocations per map.
2. **False oracle passes (finding 3).** Every gate script must assert `scf.driver`. R1 is a prerequisite of the flip.
3. **Bench baseline mismatch (finding 6).** R0 must record `baseline-reference.json` before R2. After R2, the old `baseline.json` becomes the legacy baseline only.
4. **Test-count changes.**
   - R2 changes expectations, not counts.
   - R10a removes about 30 integration tests plus unit tests.
   - R10b removes about 10 more.
   - R10d removes about 30 M2 tests.
   - Each ticket reports exact before and after passed/ignored/binary counts from `run_test_binaries.py`, and lists every deleted test with a reason ("tests deleted code"). Physics oracle tests (analytic, synthetic, finite-difference) are never deleted, only moved or retargeted.
   - `artifacts/tests.json` is updated by the orchestrator.
5. **README generator.**
   - `scripts/make_readme_figures.py` reads only stored `tests/reference/**` files, so no stage changes the block as long as those files are never deleted. R10 must not delete any `tests/reference` file, even when the test that produced it is deleted.
   - The truth risk is semantic: checkpoint/restart and runtime paragraphs (finding 5), and "stored runs strict; default runs verify every step on the standard path" (still true after the flip, since ADR 0043 defines the standard path by its stop rule, not by the driver).
   - R8 must not remove v1 checkpoint and runtime code until v2 evidence exists and the user approves the paragraph change.
   - `RG-QE` checks whether the QE table numbers, printed to three significant figures, would change on the reference path. If they would, the generator states which driver produced the stored rows (D6).
6. **Legacy inputs break at the flip** (findings 4 and 7). The R1 corpus must cover every key in L9-L11 and every synthetic UPF fixture used by tests. A missing predicate turns into a refstart `InvalidInput` error at the user, not a silent wrong result.
7. **Mixing defaults change at the flip** (beta 0.3/0.5, history 6 → 0.7/8). Results differ within the standard class. Users who set `mixing_beta` explicitly get the reference meaning of beta, which is a Broyden step. Record this in ADR 0047.
8. **MPI and threads semantics.** Legacy spin/FD replicated the whole calculation across ranks; the reference path distributes it. Threads-only k parallelism (`[scf] threads`) must stay result-independent on the reference k path. Add a threads 1/4 identity test in R2 if one is missing.
9. **The oracle needs MPI for large cells.** If D2 trims the oracle, si128-class oracle runs become slow. Large-cell accuracy gates would then use tightened-threshold reference runs (as G5 did at 5e-13) cross-checked against the oracle on smaller cells. This must be stated in ADR 0047.
10. **Frozen files.** `docs/schemas/` is frozen. If `results.json` has a schema with `additionalProperties: false`, the R1 fields need a schema-owner exception in the ticket. The root `Cargo.toml` and `ownership.toml` are not touched by this plan.

---

## 6. Decisions needed from the user (recommendation in brackets)

- **D1 — runtime control** (checkpoint/restart, observe/agent, self-healing, MCP, Python step): port through a reference stepper and checkpoint v2; keep the legacy session for runtime only (this blocks R10b); or retire, which needs README block changes. [Port, R8; legacy stays as the routed fallback until R8 lands.]
- **D2 — oracle scope:** move it whole, or trim it to serial/threaded strict (−18-20k lines). [Move whole first (R9). Trim later as R10d after a 1e-12 Ha oracle-versus-oracle check.]
- **D3 — sunset rule** for default-on rollback switches on the reference path. [Two dev08 evidence generations without rollback use.]
- **D4 — legacy-only input keys after R10:** error with a migration message, or permanent fallback. [Error.]
- **D5 — `RELAX_EXTRAPOLATION`:** port or retire. [Decide from the ledger; it is default off.]
- **D6 — README factual corrections and driver provenance** in the generated block after `RG-QE`. [Report first, then edit with approval.]
- **G4D free-moment PBE** stays a research ticket. Its fallback is the last legacy standard route. [Keep the fallback; R10b deletes everything else and keeps `collinear`/`finite_temperature` under the oracle.]

### Critical files for implementation
- `crates/dft-cli/src/run.rs` (routing :123-176, the `execute_mode_inner` legacy and reference branches :1479-2060, `reference_electronic_run` :3004)
- `crates/dft-cli/src/config.rs` (raw switch reads :854, :1227; `scf_path_with` :442)
- `crates/dft-engine/src/scf.rs` (nested `refscf` :1-4, shared types, `run_impl` → `ScfSession` :4195)
- `crates/dft-engine/src/refscf/delivered.rs` and `crates/dft-engine/src/refstart.rs` (facade, seeds, `strain_reference` and pw92 rejection :398-422)
- `crates/dft-numerics/src/switches.rs` (186-entry `table!` :194-740, `REFERENCE_SCF`/`STRICT` :450-453)
- Also: `crates/dft-cli/src/relax.rs`, `crates/dft-cli/src/bands.rs`, `crates/dft-engine/src/nscf.rs`, `crates/baandi-py/src/lib.rs`, `crates/baandi-mcp/src/lib.rs`, `scripts/make_readme_figures.py`, `scripts/bench_scf.py`, `docs/performance/baseline.json`


## 7. User decisions (2026-11-25)

- D1: approved as recommended (port checkpoint/restart, runtime control, Python/MCP stepping to a reference stepper, R8).
- D2: approved as recommended: move the strict oracle whole (R9); trimming its MPI part is decided separately at the end (R10d) after a 1e-12 Ha oracle-versus-oracle check.
- D3: removal approved, more broadly than the recommendation: default-on rollback switches on the reference path are removed together with their pre-adoption code (no two-generation wait; "they can be recreated quickly if needed"), and rejected default-off experiment switches are removed. Kept: diagnostics/profiling, infrastructure, REFERENCE_STRICT, and isolated algorithm modules protected by the Generality rule (e.g. RMM, reference ladder). This supersedes, for stage R, the AGENTS.md wording that adopted defaults stay switchable for rollback.
- D4: legacy-only input keys are removed outright after R10 (no migration messages: there are no existing users). They become ordinary unknown-key errors.
- D5: approved as recommended (decide port or retire of RELAX_EXTRAPOLATION from the attempts ledger).
- D6: approved as recommended (report README factual corrections first, edit after approval).

## 8. Evidence log

- gen14 (d81ef8af, R2 default flip): RG-ACC at the default, all 8 cases within 5e-9 Ha of the strict oracle and VERIFIED; si8 forces 8.2e-6 Ha/bohr, stress 1.1e-9; 3-step relax at 5e-13 within gates. RG-QE: 5 of 6 README cases unchanged (<= 4e-10 Ha/atom); PBE Gamma failed on a transient negative density (fixed by R1B screening default).
- gen19 (22531ab1: R1B, R4b, R7b, R8a, R8b): R8a si128 MPI8 energy bits, density hash, maps, iterations and H columns identical to 5d787b70 (verify off and on), wall 1.021x / 1.005x single samples, inside the precise-campaign range; bands/dos default both stages reference with all writers; PBE Gamma README case VERIFIED, E - E_QE +1.6983e-6 Ha/atom (README +1.6982e-6); O atom triplet PBE+NLCC converged and VERIFIED in 9 iterations with screening (legacy/oracle did not converge in 200); Si variable-cell relax default vs legacy final lattice 6.3e-8 bohr, pressure 2.3e-9 Ha/bohr^3; si32 checkpoint at map 5 + restart bit-identical to the uninterrupted run. Note for R12 docs: `restart` writes an export container, not a results directory.
- gen15 bcc Fe: reference VERIFIED in 75.8 s; legacy and strict oracle did not finish in 45 min (no agreement claim).
Evidence: dev08 /data/sig-dft/gen14, gen15, gen19 RESULT.md.
- gen20/gen21: R10A changed the strict si8 Gamma oracle (explicit cutoff_ladder = "off" was ignored; dev08 bisect /data/sig-dft/bisect20); fixed at ab77ec10, RG-ORC five cases bit-identical to gen19 again (MPI8 and np1, forces/stress). R8c observe run bit-identical to plain. Moment FD defaults: bcc Fe bit-identical to gen15; O atom FD fixed moment (PBE+NLCC, vacuum) converged and VERIFIED after G4F (legacy/oracle stop at the growth limit). Evidence: /data/sig-dft/gen20, gen21 RESULT.md.
