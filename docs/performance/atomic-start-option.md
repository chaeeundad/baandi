# P46A opt-in atomic starting wavefunctions

Status: R2 rejects recommendation on si-k222; mandatory stop after the third
new system. The architect-authorized opt-in implementation remains default off.

## Preregistered R2 (2026-11-08, before new data)

The architect explicitly authorizes a default-off production option despite
P45A's rejection under its different frozen rule. P45A remains a negative
result under that rule; it is not retroactively relabelled a pass. The new
premise is an opt-in option and an independently tested rule, not tuning the
existing default or the P45A data.

New systems: si32-gamma (real), si-k444s, si-k222, material-mgo-40,
material-si-20, and si-skew from benchmarks/cases (distinct from P45A's
externally extracted 16-atom skewed cell). Spin and Fermi-Dirac fixtures test
fallback only. Si8 and Si16 are consistency controls only.

Paired free-electron/atomic standard starts use mean-centered relative density
noise 1e-10 and new LCG seeds 301 through 310. Two unperturbed strict reference
levels use (density, energy, eigenvalue) tolerances (1e-9,1e-11,1e-10) and
(1e-10,1e-12,1e-11), energies in Ha. References must certify and their energy
plateau must be <=1e-12 Ha. The ticket also specifies one reference/five starts
for si32; the architect clarification below resolves this before measurement.

Recommend the option only if every condition passes:

- Median SCF does not rise on any system.
- Where default spread (max minus min) is >=4, atomic spread is <=half.
- Median whole-SCF single-thread wall time, including atomic initialization,
  rises no more than 5% on any system and falls on the median over systems
  (architect clarification below). H applications and MACs are diagnostics.
- No default-convergent start fails; strict still certifies.
- Maximum delivered standard energy error is <=5e-9 Ha and <=3 times default
  on every system. This replaces P45A's 1.5 ratio before new data.
- Delivered strict accuracy remains within 1.5 times default.

Stop at the first failed condition, identify the system, and do not tune.
Missing coverage or undefined accuracy ratios cannot pass. H applications and
complex MACs have different units: record both separately, together with MGS,
radial work and EVD order. No conversion of MACs to H applications is specified
in the original ticket. The architect subsequently selected whole-SCF wall
time before any new data, as recorded below. H-only improvement is not
full-cost acceptance. The original clarification request is retained in history.

## Implementation contract

`[scf] start_wavefunctions = "free_electron" | "atomic"`, default free_electron.
Use the P45A Fourier coefficients, all supplied orbitals including zero
occupation, two-pass MGS with the existing relative rank threshold, and a
Rayleigh-Ritz rotation under the first input-density Hamiltonian. Starting
density, mixing and tolerances remain unchanged. Supplied initial vectors win.
Missing orbitals, spin and FD require explicit fallback metadata. Old
checkpoints load as free_electron; exact restart compares the option.

At Gamma, real spherical harmonics have parity (-1)^l; multiplying by (-i)^l
and the structure factor gives c(-G)=conj(c(G)). A real reduced Hamiltonian
preserves this under Ritz rotation. Plane-wave completion must use real sine
and cosine pairs (and the real G=0 column), not individual complex waves.
The existing HalfSphere/ORC bridge must verify Hermitian completion.

Evidence belongs in P46A-evidence/, never target/. Use one target directory.
No UPFs, artifacts updates, remote commands, default changes or QE source.
Mac counts are observations, not dev08 speed acceptance.

## Implemented representation and bookkeeping

`atomic_start.rs` is a child module of `scf`, re-exporting the option and work
report. Parsed orbitals are immutable full-mesh samples; missing occupation or
label remains `None` rather than invented data. Angular momentum above 3 is
retained by the parser but the atomic option explicitly falls back. Every
supplied supported orbital, including zero occupation, enters the span.

The Fourier coefficients are
`(-i)^l exp(-i (G+k).tau) Y_lm((G+k)hat) / sqrt(Omega)` times the existing
projector transform `4*pi*integral r*chi(r)*j_l(|G+k|r) dr`. Exact squared
momentum bits key the per-orbital radial cache. Species/orbital/position/magnetic
component ordering matches P45A. Rank-deficient atomic columns reject without
retuning the library's relative 1e-12 MGS rank test. Plane-wave completion is
ordered by kinetic energy with stable ties. At Gamma it uses real cosine/sine
pairs, including an unpaired final component for odd band counts. The existing
Gamma bridge packs supplied full-sphere coefficients with HalfSphere's 1e-14
Hermitian check and passes ORC vectors to real Davidson; no bridge tolerance
or free-electron solver path was changed.

P45A never exercised plane-wave completion. Its append call lacked reserved
column capacity; the new completion fixture detected that incompatibility.
Production reserves the requested capacity before appending. This is a setup
storage correction, not a changed normalization or rank threshold.

Atomic radial input is copied at problem construction only when supported and
requested; span construction and Ritz rotation are deferred until a cold map
has assembled the actual input-density Hamiltonian. Warm and restored vectors
therefore take precedence even for a rank-deficient candidate atomic span.
Gamma reduced imaginary matrix entries must be below 1e-12 Ha before zeroing,
exactly as in P45A. Full-rank P45A arithmetic is retained.

`atomic_start` in scalar results separately reports initial H vectors,
Gram/rotation complex MACs, MGS projection complex MACs, radial samples, EVD
orders and the fallback reason. Counts sum over k points and cold maps; the
orbital count is therefore a sum over k, not the number of distinct radial
functions. These initialization H applications are additional to the existing
history counters. Norms/divisions and EVD work are not represented as zero H
work or silently assigned a conversion factor. Preflight reports available
orbitals per k; for missing PSWFC, spin and FD it records zero used and a reason
in result metadata, including spin/FD result profiles. Supplied vectors are
identified by scalar runtime metadata. The option and accumulated startup work
survive exact continuation; old wire records omit/default to free_electron.

No edits to gamma_real.rs were necessary: its existing HalfSphere packing and
initial-vector path already satisfy the required bridge contract. No defaults,
starting densities, mixers, solver stopping criteria or tolerances changed.

## Reproduction and current limits

```sh
PATH=/opt/homebrew/opt/rustup/bin:$PATH scripts/check.sh
cargo test --release -p dft-engine --test atomic_start_contract
cargo test --release -p dft-cli --test atomic_start --test atomic_start_checkpoint
cargo test -p dft-upf --test atomic_wavefunctions
P46A_SI_UPF=/absolute/Si.pz-vbc.UPF P46A_CASES=si8-gamma,si16-gamma P46A_ARM=atomic P46A_GAMMA_REAL=off P46A_SEED_BASE=201 cargo test --release -p dft-engine --test atomic_start_option sweep -- --ignored --nocapture --test-threads=1
# Repeat with P46A_GAMMA_REAL=on for consistency controls only.
python3 P46A-evidence/threads.py
python3 P46A-evidence/identity.py after
python3 P46A-evidence/rejections.py
```

The ensemble test is ignored in debug. Its copied diagnostic loop is checked
against a native ScfSession for every atomic and free-electron trajectory,
including energy/residual/tolerance bits, Davidson/H counts and terminal
convergence. Standard results are not labelled certified. Thread runs use both
strict reference levels and standard seed 201, T=1/2/4/8/12, Si8/Si16 and both
paths. FNV-1a of every final density's little-endian f64 bytes is recorded;
small synthetic tests compare all density bits directly. These are consistency
checks, not R2 new-system evidence or remote performance acceptance.

README generation and its nine freshness tests are run; no accepted result
dataset was added, and generated README/figures remain unchanged. Artifacts
are untouched per ticket ownership. No extra target directory or UPF copy was
created. No remote or QE source access occurred.

## Verified local results

Final `scripts/check.sh` exits 0: build/test/fmt/clippy all PASS. Synthetic parser,
CLI fallback and checkpoint codec tests pass; four engine contract tests pass,
including rank rejection, warm/exact restart, missing PSWFC and 15 completion
map/path/thread combinations. Existing parser tests are unchanged. Earlier
compile/fixture/format failures and their fixes are retained in `commands.json`
and individual logs; they are not numerical R2 failures.

The 48 primary consistency SCFs comprise 40 standard starts and eight tight
references, all converged. Si8 standard is 13/13/13; Si16 is 9/9/9, separately
on complex and real Gamma. Maximum within-arm reference energy plateau is
3.979039320256561e-13 Ha. Native-session histories match the diagnostic driver.

Default pre-hook versus option-build identity: 70 comparisons across five cases,
T=1/2/4/8/12, strict/standard, and all applicable real off/on paths; every energy
difference is 0 Ha and every density SHA-256 matches. The baseline executable
was saved before engine/CLI hooks, after the additive parser accessor (which
changes no radial numerical data). Thirty forced-real multi-k runs reject as
expected; these are not counted as successful SCFs. Default densities also
match across thread budgets within each case/profile/path.

Atomic complete-SCF density fingerprints match for all 48 T>1 versus T=1
comparisons: Si8/Si16, both paths, two tight references and standard seed 201.
The test-only FNV-1a fingerprint covers all f64 density bits; it is not a
cryptographic hash. Synthetic tests additionally compare density bits directly.
The evidence verifier exits 0 and produces `verified-summary.json`. This does
not establish R2 generality on the new systems or a total-cost recommendation.

Historical pause: implementation verification completed before new-system data.
The sampling and cost-rule questions were resolved by the architect below;
there are no remaining clarification requests.

## Architect clarification before R2 data (2026-11-08)

The architect resolved the two pending points without changing R2 thresholds:
si32 uses five standard seeds 301–305 and one strict reference per variant;
all other systems use ten standard seeds 301–310 and two strict references.
For si32 the single reference is the tighter already specified level (run 1).
There is no within-arm plateau estimate from a single reference; report that
limitation and the cross-arm reference difference, never fabricate a plateau.

The cost gate is the ratio of median **single-thread whole-SCF wall times**
(atomic/free_electron), bounded by 1.05 per system and decreasing on the median
across systems. H vectors and MACs remain separate diagnostics. Use one build,
one machine, sequential processes and alternating AB/BA order for paired seeds.
No build or other task-owned measurement overlaps timed SCFs. Time the native
ScfSession from problem construction through its terminal step, including all
atomic preparation and initialization, certification and mixing. Input parsing,
noise generation and post-run printing are outside this interval. An untimed
source problem supplies the same density grid to both arms; no atomic Fourier
span or Ritz rotation occurs in it. No duplicate diagnostic-driver SCF runs
inside or next to the measured native SCF.

System order is si32-gamma (real), si-k444s, si-k222, material-mgo-40,
material-si-20, si-skew (real). Collect reference pairs before standard pairs;
check terminal failure/certification immediately and per-system distribution,
wall and maximum-energy-error criteria after its planned starts. Stop at the
first failed condition, without retuning or running subsequent systems.
Standard errors use each arm's tight reference as in the original protocol;
retain cross-arm reference differences. With two levels, report strict energy
plateau errors and their ratio without an additive floor. The si32 single-level
exception cannot independently estimate delivered strict error; certification
and cross-arm agreement are reported without claiming that unavailable ratio.
Local Mac timings are observations; dev08 confirmation remains outstanding.

## R2 outcome: stop on si-k222 (2026-11-08)

**Do not recommend the option under R2.** On si-k222 the maximum atomic standard
energy error is 8.8816420884541e-10 Ha, versus 5.6632920575339085e-11 Ha for
free_electron: **15.682825463 times default**, exceeding the frozen factor 3.
It remains below the separate 5e-9 Ha absolute bound. Lower SCF count and wall
time do not override the independent accuracy condition. No tolerance, seed,
band count, mixing or implementation parameter was retuned.

Each completed system was evaluated after its planned paired ensemble; the
collector wrote `r2/STOP.json` on the first failed system and refuses to start
another system while that marker exists. Material MgO, PseudoDojo Si and the
second low-symmetry Gamma case were **not measured**, not counted as passes.
No aggregate recommendation across the complete system set is available.

| System / path | Free SCF min/median/max | Atomic SCF | Free median wall (s) | Atomic median wall (s) | Wall change | Maximum free / atomic standard energy error (Ha) | Decision |
|---|---|---|---:|---:|---:|---|---|
| si32-gamma / real | 19/29/47 | 11/11/11 | 14.786158 | 6.205336042 | -58.0328% | 1.7722641132422723e-9 / 1.681144112808397e-11 | No failed evaluable condition; single-reference strict-error ratio unavailable |
| si-k444s / complex | 10/10/10 | 9/9/9 | 1.233611104 | 1.255154646 | +1.7464% | 2.227231732376822e-10 / 2.2710011648996442e-10 | Pass observed per-system conditions |
| si-k222 / complex | 10/10/10 | 9/9/9 | 0.3089391875 | 0.2814508955 | -8.8976% | 5.6632920575339085e-11 / 8.8816420884541e-10 | Fail relative standard energy accuracy; stop |

There are 60 completed native SCFs: 50 standard (uncertified) and 10 tight
strict references (all certified), with no numerical convergence failure.
si32 uses five paired standard seeds 301–305 and one tight reference per arm;
the other two use ten seeds 301–310 and two references per arm. Strict
cross-arm tight energy differences are 4.973799150320701e-13,
3.552713678800501e-15 and 2.220446049250313e-14 Ha respectively. For k444s,
within-arm plateaus are 2.1316282072803006e-14 / 3.197442310920451e-14 Ha;
for k222, 7.105427357601002e-15 / 1.0658141036401503e-14 Ha. Both plateau
ratios equal 1.5. These roundoff-scale two-reference checks are internal
consistency observations, not independent external error bounds. si32's
single-reference exception provides no corresponding strict-error ratio.
The negative recommendation does not depend on that unavailable ratio.

### Startup diagnostics (included in the measured wall interval)

| System | Median free / atomic total H vectors | Initial atomic H vectors | Gram/rotation complex MACs | MGS projection complex MACs | Radial samples |
|---|---:|---:|---:|---:|---:|
| si32-gamma | 4665 / 2059 | 128 | 297970176 | 386145024 | 344800 |
| si-k444s | 6293 / 5944 | 256 | 3068928 | 2685312 | 12567960 |
| si-k222 | 1600 / 1449 | 64 | 766080 | 670320 | 732700 |

The Si32 EVD order is 128; the two mesh cases perform an order-8 initial EVD
at each k point (32 and eight points). H totals include ordinary, certificate
and startup applications. MAC counts are diagnostics, never added to H counts
as if they shared units. The k444s startup radial work helps explain why fewer
H vectors need not reduce wall time; no causal timing decomposition is claimed.

### Timing and evidence limits

All paired samples used the same release test binary SHA-256
`de51ad83f163ea575b0bc35ef3bda0a0c631880e002bc450325028577df4ccb2`, T=1, one
Mac, alternating AB/BA order, and no overlapping task-owned builds or runs.
Recorded 1-minute host load reached 4.6411 (si32), 2.9033 (k444s) and 2.6631
(k222). This was not an isolated quiet node; especially the +1.75% k444s
change is a local observation, not a resolved speed claim. dev08 acceptance
was not attempted. The decisive energy ratio is independent of noisy wall time.

Raw records, individual stdout logs, command environments, exit codes and
per-system summaries live in `P46A-evidence/r2/`. `preregistered.md` preserves
the exact protocol before new data; `preregistration.json` hashes it and the
binary. `input-source-hashes.json` identifies the measured harness and source.
`r2.py` executes one system at a time and enforces the stop marker. No new QE
run or reference claim, new default, artifact update or UPF copy was made.

Post-measurement verification: `scripts/check.sh` exits 0 with build/test/fmt/
clippy PASS; `r2/verify.py` independently reconstructs all 60 terminal records,
arm alternation, timing/energy/work statistics and the stop decision (exit 0).
The earlier default-identity/atomic-thread evidence verifier also exits 0.
Production source hashes are unchanged since those implementation gates; this
measurement continuation edits only the test harness, evidence and documents.
README regeneration/freshness remains unchanged (nine tests), artifacts are
untouched, and `git diff --check` exits 0. All changes remain uncommitted.
