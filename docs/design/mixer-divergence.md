# G-MIXER-DIVERGENCE: SCF control diagnosis and bounded experiments

Source baseline: `6a060639934b60859c1a70f6bd8a918b9aac315f`, 2026-11-23.
Mac results below are **work counts and correctness checks, never timing evidence**.
The ticket's si128/si216 histories are relayed observations, not reproduced dev08
measurements. No QE reference was generated in this task.

## Scope and current implementation

Implemented: `BAANDI_MIX_CELL_KERKER=on|off`, default **off**, a geometry-derived
Kerker residual filter that reuses the existing FFT and Pulay implementation.
The default arithmetic, standard stop rule and fixed-point verification remain
unchanged. Configured beta was already applied correctly; its record semantics
are documented and tested, without multiplying the proposal by beta a second time.

The session-level recovery and eigen-schedule experiments described below are
**designs, not implemented switches**. In this checkout the iteration driver is
`crates/dft-engine/src/session.rs`, not `scf.rs`. The ticket owns `scf.rs`'s
mixing/control sections but excludes `session.rs`. A scope-extension question
was submitted; no response had arrived when these independent experiments were
prepared. No edits to that file are included. Likewise the existing Kerker
wrapper's diagnostic forwarding lives in unowned `mixing/kerker.rs`.

## Answers to the five control questions

1. **Beta is applied.** `mixing.rs`, `PulayMixer::mix`, initializes
   `next = x + beta R`; both its ordinary and incremental updates subtract
   `(Delta x + beta Delta R) gamma`. `mix_dist.rs::mix_block_with` implements
   exactly the same expression on the owned block. `session.rs::make_mixer`
   forwards `Mixing::Pulay { beta, history }` to both constructors. There is no
   missing beta on the M2 path. `IterationRecord::mixing_fraction` is a separate
   **feasibility fraction**, initialized to 1 in the session and changed only
   when negative samples in a proposed density require convex interpolation.
   The configured beta is in `scf.mixing.beta`; there is no unique scalar
   effective beta for an Anderson step that is not parallel to the raw residual.
   The new scalar one-secant oracle tests beta 0.2/0.5/0.8 in both update modes;
   distributed first-step tests cover chunk/plane layouts and ranks 1/2/3.

2. **The reset counters are not an energy watchdog.** `session.rs` resets for
   spectral rescue, the opt-in `MIX_GROWTH_RESET`, or positivity feasibility.
   A positive, finite density may have a very large Hartree error and still
   accept the full proposal. Neither `delta_energy_ha` nor `standard_error_ha`
   triggers the ordinary reset. The default fatal divergence branch returns an
   error instead of resetting/recovering. Thus +5.74 Ha alone cannot fire it.

3. **`predicted_ratio` is diagnostic only.** For an unweighted fit it is
   `sqrt(max(0, R.R - gamma.rhs)/(R.R))`, the norm of the least-squares
   **projected residual** relative to the current residual. It is not a measured
   next-map contraction or a condition number. A value near zero means the
   residual was well explained by the history and is normal, for example at
   exact termination of a linear problem. No production control reads this
   field. The separate default-off Coulomb guard computes
   `coulomb_predicted_ratio` of the L2 projection and refits *this step's*
   coefficients only if its Coulomb norm exceeds the raw norm. It neither
   clears the history nor verifies the nonlinear next density. Under a refit,
   the generic `predicted_ratio` describes the final Coulomb fit. Counting a
   small L2 ratio alone as divergence would incorrectly reset successful fits.

4. **The standard eigen schedule is deliberately monotone.**
   `StandardProfile::ratchet` returns
   `clamp(min(previous_tolerance, 0.1 * previous_standard_error / max(Ne,1)),
   5e-14, 5e-3) Ha`. `session.rs` supplies the last record; the first tolerance
   can be seeded by the cutoff ladder. When error increases, the new error term
   is larger than the historical minimum and the previous tolerance wins.
   At si128 map 1, `0.1 * 0.03284 / 512 = 6.4140625e-6 Ha` (rounded relay),
   explaining map 2 and the plateau through map 9. This is not a beta bug or a
   hidden freeze counter. The optional stall/density-forcing rules are off;
   they have prior negative generalization evidence. Tightening the inner
   solve is an ablation, not an established cure for charge sloshing.

5. **The default fatal test watches a different norm.** `density_distance`
   is `sqrt((Omega/Ngrid) sum_r (rho_out-rho_in)^2) / Ne`, reported in
   Bohr^-3/2. At every ordinary map the driver tests
   `r_k > divergence_factor * max(best_previous_r, density_tolerance)`, counts
   consecutive successes, updates the best residual, and aborts after
   `divergence_patience` successes. With the defaults this requires **three
   successive residuals above 20 times the best L2 norm**, not a 20-fold
   increase in the Coulomb energy. The relayed map-4 ratio is only
   `4.489e-4 / 1.575e-4 = 2.85016`; it cannot trigger the default branch.
   `standard_error_ha = (Omega/2) sum_(G != 0) (4*pi/G^2)|R_G|^2` is an
   extensive, strongly low-G-weighted quadratic quantity. It can grow 272x
   while the normalized L2 norm grows less than 3x.

Exact code references in this working tree:

| Question | Source anchors |
|---|---|
| beta / reported fraction | `session.rs:379` constructor forwarding, `mixing.rs:626` and `mixing.rs:716` beta terms, `mix_dist.rs:609` distributed update, `scf.rs:372` record contract, `session.rs:837` initialization and `session.rs:1521`, `session.rs:1596` damping assignments |
| reset causes | `session.rs:1428` spectral rescue, `session.rs:1441` optional growth reset, `session.rs:1516` and `session.rs:1591` feasibility resets |
| projected ratio | `mixing.rs:759`, `mix_dist.rs:702`; distinct Coulomb check `mixing.rs:502`, `mix_dist.rs:197` |
| eigen tolerance | `session.rs:582` history-driven schedule; `convergence.rs:40` ratchet formula |
| fatal divergence norm | `session.rs:921` comparison, `session.rs:930` patience/abort; `scf.rs:4287` normalized L2; `convergence.rs:295` Coulomb weights |

All engine filenames above are under `crates/dft-engine/src/`; `mixing.rs` is
under `crates/dft-solvers/src/`. The machine-searchable anchors are also in
`../performance/lanes/MIXER/source-anchors.txt`.

## Existing experiments and the changed premise

- Ledger `MIX-GROWTH-RESET`: factor 2.5/3/5 resets did not cure the old
  elongated/slab plateau and often worsened the best residual. We do not
  relabel that failed mechanism as a new result.
- `MIX-COULOMB-GUARD` and `MIX-COULOMB-GUARD-DIST`: early trajectory gains were
  diluted or reversed by extra verification attempts. The guard is already
  available on replicated, distributed chunk and owned-plane Pulay. The
  repository records candidates; it does **not** contain a successful dev08
  promotion that would justify enabling it. Its off/parked state is consistent
  with case-dependent full-run evidence, not proof that no port exists.
- `SCF-TAIL-PREDICTOR` (`STOP_PREDICTOR`): the proposed cross estimate is not
  conservative (observed estimate/true error as low as 0.654); a safety factor
  erased most of the gain. `on` can change acceptance and is excluded here.
  All experiments explicitly pin it **off**.
- P40A mixing survey: no universal beta/history/Hartree-metric winner; old
  elongated/slab cases failed across all arms. A31's fixed-q Kerker experiments
  also failed as universal defaults. New premise: a deterministic cell-length
  rule, exact cubic inactivity, and a normally convergent 1x1x4 supercell that
  reproduces the severe low-G divergence. No scalar q0 sweep was performed.
- [Issue #21](https://github.com/chaeeundad/baandi/issues/21), read through the
  GitHub connector: a response model requires a projected Sternheimer action,
  finite-difference validation, refresh/invalidation policy and a measured
  response H/FFT/reduction budget. A geometry filter is **not** that response
  primitive. This bounded ticket does not claim to implement LRDM or chi0.

## Independent candidate designs

### C: cell-derived residual filter (implemented)

Let `gmin^2` be the smallest positive squared norm in the actual density basis,
including the frozen-basis convention where applicable. Let
`gcube = 2*pi/Omega^(1/3)`. Set

`q0^2 = max(0, gcube^2 - gmin^2)`,
`P(G) = |G|^2/(|G|^2 + q0^2)` for nonzero G, and `P(0)=0`.

The denominator of the longest represented mode is brought up to the
frequency scale of the equal-volume cube. There is no material, atom-count,
seed, band-count or fitted threshold in this rule. It scales as an inverse
length and uses actual reciprocal modes rather than reciprocal-column norms,
which can be wrong under a skew re-basing. A 64-epsilon relative comparison
only suppresses roundoff at the cubic equality; it is not a physical threshold.
A cubic cell takes the original Pulay branch without an FFT round trip.

The existing `Kerker` adapter handles signed residuals, FFT normalization,
zero mode and skew-grid Nyquist conjugacy. `PreconditionedPulay` fits the
filtered map; the raw density residual and Coulomb stop/verification still
measure the physical SCF map. Beta/history are retained. Explicit Linear or
Kerker input is not overridden. The switch resolves once during problem
construction and its effective q0 is recorded in the run path note. The CLI
`scf.mixing` block still describes the requested input mixer; the switch path
note distinguishes the effective Kerker mixer from that request.

**MPI limitation:** C uses the existing replicated Kerker-Pulay implementation.
M2 Davidson and grid stages remain enabled, but the plain distributed mixer
and owned-density loop become ineligible. This can spend the saved maps on
FFTs, replicated mixing and density gathers. Count reductions do not establish
a speed gain. A distributed filter must be a separately tested follow-up;
this patch does not silently change `M2_MIX`'s default.

**Diagnostics limitation:** the existing Kerker wrapper does not forward the
inner Pulay diagnostic, so `history[].pulay` is null under C. Reset, feasibility,
raw residual, standard error and tolerance remain observable. Null is not a
zero predicted ratio. Adding that forwarding requires the unowned wrapper file.

### G: observed model-failure recovery (design; pending session ownership)

Use the last accepted best pair `(x_best, F(x_best))`, the raw L2 norm and
Coulomb error, plus the *previous* proposal's projected ratio `p`. Trigger
above the unchanged class floor when actual error increases and the relative
loss exceeds the promised relative squared-norm gain `1-p^2`. This compares
measured loss to model prediction; there is no Si-specific growth constant.
A small p by itself never triggers. For strict maps without a Coulomb error,
use squared raw residual ratios and label the different metric.

On a trigger discard the rejected secants, retain the best pair, and take
one simple step from that pair. A conservative candidate damping is
`beta_safe = min(beta, 1/(1+Lsec))`, where `Lsec = ||Delta R||/||Delta x||`
in the chosen metric on the actual rejected step. This is a local measured
response estimate, not a global contraction proof; zero or nonfinite secants
must fail closed to a valid saved pair or the existing fatal exit. Preserve
charge, apply the existing positivity safeguard and consume the ordinary map
budget. Next map rebuilds history. Do not return the saved state as converged;
only the existing final certificate may accept a state.

Proposed switch `BAANDI_MIX_DIVERGENCE=off|trace|on` (not registered yet).
Trace must preserve all density/energy bits. On must update `mixer_reset` and
`mixing_resets` in both owned and replicated branches, retain the origin of a
rollback in diagnostics, and reject exact checkpoints until controller state
is serialized. Allocate the two best-pair buffers at initialization, never in
Hamiltonian/eigensolver loops. MPI decisions use canonical global scalar
reductions and best pairs stay on owned planes. Tests: stable linear maps with
p approaching zero, hidden low-G growth, nonlinear false predictions, repeated
rejections, zero residual, positivity, rank agreement and checkpoint rejection.

### E: growth-reactive eigen cap (design; pending session ownership)

Proposed `BAANDI_MIX_EIGEN_REACT=off|on` (not registered yet): after two
completed standard maps set

`a = max(1, r_k/r_(k-1), sqrt(E_k/E_(k-1)))`,
`tau_next = max(FLOOR_HA, min(ordinary_ratchet, tau_k/a))`.

Finite positive denominators are required. It is dimensionally consistent,
never loosens the existing tolerance, and acts only on observed growth; it
uses no atom/material parameter. It is deliberately different from P29B's
always-on dimensional `r^2` cap. It may oversolve genuine nonlinear sloshing,
so total H is a rejection counter even if maps decrease. Verification uses
its existing rule. Required tests include monotone/no-growth identity,
first-map identity, floor handling, missing/zero/nonfinite observations and
one-step vs repeated growth. Test G, E, C independently before G+C, E+C and
G+E+C; the implemented matrix below covers C and the existing Coulomb guard,
not unimplemented G/E combinations. C plus the existing Coulomb guard was also run explicitly.

## Measurement and reproducibility

Nine systems, seeds 601/602/603, paired fractional displacement +/-1e-3 from
Python `random.Random(seed)` in atom/component order. One process and one
thread; original benchmark cutoffs; Pulay beta 0.5/history 6; production ladder,
standard 5e-9 Ha and fixed-point verification; 100 fine-map limit. `si-long4`
is generated from four conventional Si8 cells along z (32 atoms, 30 Ry).
Synthetic Gaussian-fcc (s/s/p) and p3 (three p projectors) are nonphysical
adapter/solver fixtures, not Si material evidence. UPFs remain external/scratch;
records contain their SHA-256, never contents.

Strict references use the same displaced inputs, `BAANDI_REFERENCE_STRICT=on`,
density tolerance 1e-10 Bohr^-3/2, energy 1e-12 Ha and eigen residual tolerance
1e-11 Ha. Failed references have **no accuracy oracle**, even if an energy was
printed. No tolerance is relaxed to turn a failure into a reference. The two
old stress cases must remain failures in the denominator.

Run `docs/performance/lanes/MIXER/reproduce.py --help` for the collector.
Exact invocation, binary/input/UPF hashes, exit codes, full fine histories and
certificate attempts are in the per-run JSON files. The first reference
collection used the invalid CLI key `eigen_tolerance_ha` (exit 1), corrected to
`eigen_tolerance` and rerun; these input errors are retained separately and are
not numerical failures. Initial bad reference records are never oracles.

H work is reported by stage. Verification work is already in each history's
`fixed_point_certificate`; never add `verification.attempts[].h_vectors` again.
The full known H total includes ordinary fine H, certificate H, coarse H and
coarse atomic-start H (fine atomic-start H separately when present). A result
file absent after timeout/error means unavailable counters, not zero work.
All failed attempts and unqualified energy comparisons are retained.

Generated tables and numerical verification summary:
[Mac evidence](../performance/lanes/MIXER/evidence.md).
Remote protocol: [dev08 request](../performance/lanes/MIXER/dev08-request.md).

## Outcome and verification limits

C alone is the supported Mac candidate. On the normally convergent 1x1x4
fixture, fine maps across three seeds sum to 142 -> 26 and total recorded H to
40,534 -> 7,585. Si8 is bit-identical, Si16/Si32 improve, and the skew and two
synthetic fixtures are unchanged. Across the 21 qualified standard runs, fine
maps total 296 -> 164 and H 59,397 -> 22,640; these totals are conditional on
this fixture set and exclude the six failed stress cases. Maximum standard
energy difference from the qualified tight reference is 1.673591e-9 Ha;
strict-C vs strict-A is at most 6.252776e-13 Ha. C+the existing guard adds no
benefit where C is active and worsens the skew fixture; it is not recommended.

The old elongated/slab fixtures still do not converge in 100 maps. Their
stagnation is distinct from the recoverable early sloshing in a normal
supercell; this task does not resolve their physical-map/occupation issue or
claim that one cause has been proved. No si128 Mac run was necessary: its MPI8
performance and exact trajectory remain the remote experiment.

Final `scripts/check.sh full`: exit **1**, build/fmt/clippy PASS; **798 tests
pass, 3 fail, 181 ignored**. All three failures are pre-existing pseudopotential
cache tests rejecting the required in-repository scratch path. The permission
profile permits no external scratch on the 990 volume, and no safety check was
bypassed. An earlier full run also failed fmt because tests were added while it
was running; the final run's fmt gate passes. README figure tests exit **1**:
8/9 pass, with the existing structural-scope assertion still expecting the
obsolete phrase `standard profile is unsupported`. README and its generator
are outside this ticket; no accepted/default result was published there.

New geometry, beta and distributed beta unit tests pass in the full engine/
solver tests. The release geometry test exits 0; the M2 si16 simulated P2/P3
smoke exits 0 with 8 fine maps/36 Davidson iterations/619 fine H vectors each,
energy within 5.471e-13 Ha of serial. The two rank counts differ by one energy
ulp and different density hashes: no rank-count bit-identity claim. Thread
budgets 1/4 do give identical numerical traces in all 27 C pairs. Final-binary
OFF checks match the original binary on all 9 systems (seed 601); the earlier
full 27-pair off comparison is also exact. Runner contract tests (fake local
launchers, no remote action) pass 3/3. See `lanes/MIXER/validation.json`.

The live session registry/dashboard, README evidence generator and wrapper
integration are outside OWNED_FILES. The master owns those updates and the
real MPI8 run. No recurring session timer or dashboard upload was started by
this implementer. Changes remain uncommitted.
