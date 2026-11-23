# R-SCF-2 validation: runnable reference SCF, acceptance incomplete

**Status: FAILED acceptance.** The new path runs end to end and all 24 tested
case/rank combinations converge with a passed independent fixed-point check
and energy error below 5e-9 Ha. The Si16 rank-trajectory gate fails at the final
SCF step. No default changes, timing claims or dev08 measurements are made.
Changes remain uncommitted.

Source base: `5c3c9db49fbdb754e7f9bfcc76ddb2cebb8e4b0d`, Mac ARM64, pure Rust.
[Machine-readable evidence](patch2-validation.json) contains source/binary and
UPF hashes, every case/rank result, rank deltas, QE stdout hashes and baseline
identity checks. [Full Si8 reference records](patch2-si8-reference.json) retain
all per-step work, energies, phases and verification attempts. Phase seconds
are observations, not a speed comparison. No UPF is included.

## Implemented boundary

`driver="reference"` or `BAANDI_REFERENCE_SCF=on` branches before constructing
a legacy problem. It uses final-cutoff `refstart`, the owned-row
`ref_davidson`, the Hartree residual metric and patch-1 difference-history
Broyden. Orbitals warm-start the next map. No coarse ladder, session or legacy
solver dispatcher is constructed, including for the independent check.
`mixing="broyden"`, beta and history are accepted through both selectors;
legacy broyden remains rejected. Defaults are unchanged.

Default standard accuracy 5e-9 Ha maps explicitly to `conv_thr=1e-8 Ry`.
Requested/used ethr, first-step redo, each inner call, active widths/sum, H
columns, reduced orders, restarts, unconverged roots, history condition,
output and variational energies, and phase seconds are recorded. Per-band
Ritz changes expose the failed rank branch. The pure fixed-point assessment
is shared with `session.rs`; legacy clipping/predictor decisions are retained.
A failed verification continues the same Broyden trajectory (ADR 0043 V5),
without adopting the verification output or resetting history.

On the displaced Si8 generality fixture the first estimator stop is step 10:
dr2=8.390015129454428e-9 Ry but verification error=3.067234424584117e-8 Ha.
The driver correctly continues and passes at step 12. Verification attempts
and costs are retained separately. The first implementation terminated on
this failed check; that was corrected before the delivery runs.

Memory records distinguish conservative principal-array payload envelopes
(including replicated delivery buffers) from actual process RSS. Linux VmHWM
is sampled once at output preparation and gathered across ranks. It has not
been exercised on Linux here. macOS RSS is null, not a fabricated measurement;
external rank wrappers are requested for whole-process RSS and wall including
final JSON write and MPI finalization.

## Generality and strict oracle

All cases use final-cutoff atomic start and LDA-PZ. Si8/16/32 are respectively
1x1x1, 2x1x1 and 2x2x1 conventional cells, with deterministic 1e-4 fractional
perturbations; elongated is a 3x1x1, 24-atom cell (3:1:1). Skew is the established
2-atom skew fixture from `si-skew.toml`. Silicon cutoff is 15 Ha and band count
is occupied+4. Synthetic is a two-electron Gaussian NC PP with a nonlocal
projector, three bands, cutoff 2 Ha, and an off-origin atom in a skew cell.
The test source constructs every coordinate and synthetic radial sample.

Strict legacy oracle: standard profile disabled, verification off, eigen
residual tolerance 1e-11 Ha, energy change 1e-12 Ha, density tolerance 1e-9,
Pulay beta 0.3/history 8. This is a same-engine strict oracle, **not QE**.
Each new-driver case is run at simulated rank counts 1, 2, 3 and 8.

| Case | SCF steps | Max absolute energy difference vs strict (Ha), all ranks | Max final verification Coulomb error (Ha) | Max density difference vs rank 1 (e/Bohr^3) |
|---|---:|---:|---:|---:|
| Si8 | 12 | 1.43884904e-10 | 1.77650266e-9 | 8.04911693e-16 |
| Si16 | 12 | 1.09466214e-10 | 3.08748482e-9 | **8.17472951e-8** |
| Si32 | 12 | 3.80993015e-11 | 1.02930238e-9 | 2.91433544e-16 |
| Skew Si2 | 8 | 7.59472485e-11 | 9.68938025e-11 | 3.75699472e-13 |
| Elongated Si24 | 13 | 6.06945605e-11 | 4.65785756e-9 | 3.46944695e-16 |
| Synthetic nonlocal | 10 | 1.69247072e-9 | 3.36643347e-10 | 6.90419943e-16 |

All residual guards and final charge checks pass. Maximum charge error is
3.126e-13 electrons. All five other cases pass the rank-trajectory checks.
Exploratory z=4 single-atom synthetic and sheared conventional Si8 fixtures
failed their strict legacy convergence oracle; they are not counted as
successful comparisons. The closed-shell z=2 PP and pre-existing skew Si2
fixture provide valid independent strict references. No tolerance was changed
to accept those exploratory runs.

### Unresolved Si16 rank branch

The test deliberately remains failing (release test exits 101). At SCF step
12, Davidson correction 2, the last occupied band (zero-based index 31) has:

| Simulated ranks | Absolute Ritz change (Ha) | ethr (Ha) | `< ethr` |
|---:|---:|---:|---|
| 1 | 1.9843571230637735e-12 | 1.9834923037939178e-12 | false |
| 3 | 1.9834411890684578e-12 | 1.9834923038151466e-12 | true |

A 9.159339953156e-16 Ha difference straddles the threshold. There are no
rank-repair/drop events. Rank 1 takes 3 corrections/93 H columns; rank 3 takes
2/92. The resulting dr2 differs by -3.15849833477692152e-11 Ry and the final
density differs by 8.174729508592105e-8 e/Bohr^3, beyond the 1e-9 rank test cap.
Steps 1..11 agree within the recorded rounding bounds. Both end states pass
the standard energy class and fixed-point check, but that does not satisfy
the trajectory criterion. [Raw Ritz records](patch2-si16-ritz.txt) are retained.
No threshold quantization, tolerance widening, extra forced correction, or
system-specific branch was introduced. Rank-reproducible numerical reductions
or an approved rounding-aware stopping contract remain follow-up work.

## Genuine QE Si8 comparison

The reference is the committed `tests/reference/qe/i25-si8/si8-np1/k1..k3-mix`
stdout, QE 7.6, command/input/version and UPF identity in the fixture manifest.
No QE execution or source reading occurred on the Mac. Cell, displaced
coordinates, 30 Ry cutoff, 36^3 grid, 20 bands, atomic start, beta 0.7/history 8
match. The adjacent `si8-qe-replay.toml` is used with the external UPF path and
explicit band/mixing fields. These are one-core trajectory observations; no
MPI/threads timing comparison is present.

| Step | QE printed dr2 (Ry) | BAANDI dr2 (Ry) | QE printed ethr (Ry) | BAANDI ethr (Ry) | QE / BAANDI Davidson iterations |
|---:|---:|---:|---:|---:|---:|
| 1 | 0.43855533 | 0.4347506747835862 | 1.00e-2 | 1.00e-2 | 3 / 3 |
| 2 | 0.01956404 | 0.01980782493916525 | 1.37e-3 | 1.358595858698707e-3 | 1 / 1 |
| 3 | 0.00085399 | 0.000881418344485175 | 6.11e-5 | 6.189945293489141e-5 | 4 / 4 |

| Step | QE printed total energy (Ry) | BAANDI variational energy (Ry) |
|---:|---:|---:|
| 1 | -62.61942841 | -62.62000017134319 |
| 2 | -62.64131413 | -62.64154670773554 |
| 3 | -62.64319427 | -62.64316797941382 |

The first reported difference is dr2/energy at **step 1**, before history can
contribute. A separate startup diagnostic recovers QE rho_in(1) from
`(rho_mix(1)-0.7*rho_out(1))/0.3` and compares the actual `refstart` density in
identical Miller order: maximum **8.640452553275752e-7 e/Bohr^3**, RMS
**1.733454203850240e-8 e/Bohr^3**, worst G=(1,-1,-1); G0 difference after this
subtraction is 5.585809592645319e-16 e/Bohr^3. Thus the paths already have
different input densities before solving, and later ethr differences follow
the dr2 ratchet. Equal call counts do not imply equal approximate eigenvectors.

The existing refstart builder directly integrates radial samples and retains
its initial real-space clipping/normalization convention, whereas the supplied
QE startup note describes a radial |G| table and no clipping. This is a known
startup implementation difference; which part causes the observed coefficient
delta has **not** been isolated. Do not attribute it uniquely to Broyden or
claim QE trajectory identity. The dev08 request asks for the startup and first
solve snapshots needed to separate those effects.

Patch 1 remains accepted, as instructed. Its updated same-rank replay limit
is 1e-17 + 1e-12*amplitude e/Bohr^3, based on the relayed Si128 maximum 4.9e-18
and the old strict-limit miss 1.6e-18 at condition 6689. Cross-QE-rank diagnostic
limit is 1e-16 + 1e-12*amplitude, based on output floors 1.7e-17..6.8e-17 and
mixed-density floors 7.3e-18..1.7e-17. These limits do not relax SCF gates.

## Legacy identity and executed checks

A pristine HEAD source archive was compiled into a separate saved binary.
Against the candidate, four default-path runs (`si8-gamma`, `si16-gamma`,
`si32-gamma`, `si-skew`) have bit-identical density bytes, total energies and
eigenvalues. Entire convergence histories and verification records agree after
excluding timing fields; density hashes are stored in the evidence JSON.

| Command/check | Exit | Result |
|---|---:|---|
| `scripts/check.sh full` with worktree TMPDIR and CHECK_SCCACHE=0 | 1 | build/fmt/clippy pass; 825 tests pass, 188 ignored, 3 existing cache tests fail |
| `cargo test --offline --release -p dft-engine --test reference_driver -- --include-ignored --nocapture --test-threads=1` | 101 | synthetic passes; silicon matrix fails only the Si16 rank gate |
| `cargo test --offline --release -p dft-cli --test reference_driver -- --include-ignored --nocapture` | 0 | 2 pass: both selectors/broyden/JSON/verification, and QE startup diagnosis |
| `cargo test --offline --release -p dft-engine --test refscf -- --include-ignored` | 0 | 10 pass, including real QE Broyden replay |
| `cargo test --offline --release -p dft-engine --test refstart -- --include-ignored` | 0 | 6 pass, including Si rank matrix |
| `cargo test --offline --release -p dft-solvers --test ref_davidson` | 0 | 5 pass |
| actual CLI Si8 reference run | 0 | converged, fixed-point check passes; full records retained |
| `bench_scf.py --cases gamma --repeats 1` against pristine same-HEAD baseline, `--require-identical` | 0 | energy and density identical |
| same benchmark against historical `docs/performance/baseline.json` | 3 | rejects standard vs strict profile mismatch |
| historical baseline with `BAANDI_REFERENCE_STRICT=on --require-identical` | 5 | energy difference 3.55e-15 Ha; density differs from historical baseline |
| `python3 -I scripts/make_readme_figures.py` | 0 | generated output unchanged |
| `python3 -I scripts/tests/test_readme_figures.py` | 1 | 8 pass, existing scope-text assertion fails |
| `git diff --check` | 0 | no whitespace errors |

Set `PATH=/opt/homebrew/opt/rustup/bin:$PATH`, `TMPDIR="$PWD/target/scratch/rscf2"`
and `DFT_SI_UPF` to an **absolute** matching UPF path for release tests. Set
`TMPDIR_990="$TMPDIR" CHECK_SCCACHE=0` for the full check. An initial benchmark
invocation with a relative UPF path exited 3 because the script resolves paths
relative to its temporary input directory; absolute-path reruns are tabulated.
Initial integration errors (preflight alias and synthetic chemical label) and
clippy diagnostics were fixed; the final tests above exercise the fixes.

The three full-gate failures are the pre-existing cache tests named in
`patch2-validation.json`: the task requires scratch inside the worktree while
the cache validator requires it outside repositories. No test exclusion or
cache-policy change was used. Sandbox `ps` denial also appears in the load
report; the script continues to execute all four gates. The README assertion
still expects the obsolete phrase `standard profile is unsupported`.

The registry snapshot contains no R-SCF-2 entry and `.delegate/R-SCF-2/report.json`
is absent. Registry/dashboard files, README result-type changes and generator
code are outside OWNED_FILES; the orchestrator must register/upload this report
and handle any subsequent adoption evidence from the master checkout.
No new delegated run, remote action, git metadata write or commit was made.
The requested [dev08 jobs](run-dev08-request.md) are ready for review.
