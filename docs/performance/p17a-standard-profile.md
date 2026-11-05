# P17A: standard convergence and delivered accuracy

Implemented, additive, unpolarized Davidson only. Strict remains the default.
Baseline: `84f06a478549ce39151ba6095f5980bb652d8d26`; changes uncommitted.
No QE source was accessed. Derived from prose notes `qe-scf-cycle.md` and
`qe-solvers-and-subspace.md`, with P15A replay evidence. The architect's follow-up
allows stored QE estimates to calibrate the default without waiting for new
remote outputs. Neither an equal rounded energy nor an error estimate is a
rigorous delivered-error guarantee.

## Contracts and implementation

- `[scf] convergence_profile="strict"|"standard"`, default strict. Standard
  uses `standard_threshold_ha`, positive finite, **5e-9 Ha** by default. Supplying
  the threshold under strict rejects, preventing an ignored convergence control.
- `E_res = 2*pi*Omega*sum_(G!=0, density sphere) |delta_n_G|^2/G^2`, in Ha.
  Delta n is output minus input density. The transform is normalized by grid
  size. There is no electron-count normalization of this stopping estimate.
  Standard stops when `E_res < standard_threshold_ha`.
- The occupied Ritz-change threshold starts at **5e-3 Ha**. Thereafter it is
  `max(5e-14 Ha, min(previous, 0.1*previous_E_res/max(1,N_e)))`. The error and
  requested threshold in history reconstruct the ratchet exactly on restart.
- If iteration one's error is below `threshold*max(1,N_e)`, tighten by the same
  rule and redo the SCF map **once**, from the returned orbitals and the unchanged
  input density, before creating mixing history. Both solves' Davidson steps,
  H-vector applications and phase times are charged. Near-converged warm-start
  regression tests exercise this safeguard and exact restart preserves it.
- Roots are tested by consecutive sorted Ritz-value changes, strictly below
  their threshold. The first solve has no invented previous value. Empty roots
  use `max(5*occupied_threshold,5e-6 Ha)`. All occupied roots **and the first empty**
  must converge, including at the iteration limit. Only higher buffers may remain
  unconverged at that limit. No root locking or residual certificate is implied.
  A full-space first solve repeats its reduced solve once, because it cannot
  expand and still needs two actual Ritz values.
- There is **no fixed-point certificate solve** for standard. Results have
  `convergence_profile="standard"`, `certified=false`, thresholds, criterion,
  raw residuals, redo flags and work counters. CLI and all checkpoint kinds label
  the profile. Exact checkpoint compatibility includes both profile and metric;
  default strict checkpoint extension fields are omitted for legacy compatibility.
  A warm start records its source profile but converges under its new explicit
  configuration. Standard cannot silently enter the collinear path or Hybrid RMM.
- `[scf] mixing_metric="l2"|"hartree"` is independent and defaults to L2. Hartree
  changes only Pulay's Gram matrix and right-hand side to the Coulomb metric.
  Update formula, beta, history, rank regularization, charge projection and
  feasibility handling are unchanged. Kerker-Pulay uses the metric on its
  filtered residual; linear mixing has no least-squares metric to change.
- The Coulomb kernel, density-sphere indices, FFT plans and transform buffers
  are constructed once. Metric Gram scratch is allocated during the existing
  mixer initialization call and reused. No new logging, I/O or per-root scratch
  allocation was introduced. Existing Davidson/Pulay allocations remain.
- The existing optional final spectral audit remains independent. Its H work is
  counted; even a successful audit cannot certify a standard result. With it off,
  neither convergence profile proves lowest-root completeness. Existing NLCC and
  XC-floor checkpoint restrictions remain; this ticket does not extend them.

Strict arithmetic and tolerances are untouched. Standard-only branches perform
no metric FFTs in strict/L2. Existing strict JSON fields and schema version are
unchanged; only opted-in modes add fields. Both `verify_qe.py` and
`compare_materials_qe.py` reject non-strict or explicitly uncertified provenance
before making validation comparisons. The equal-accuracy harness is separate.

## Calibration before SCF and default freeze

`python3 -I tests/p17a_replay_calibrate.py` applies the proposed occupied/empty
thresholds to the existing P15A frozen-H traces before wiring the production
rule. Recorded native versus ratcheted Ritz prefix steps: **si8 283 -> 54**
(1000 H vectors), **k444s 1686 -> 740** (5631 H vectors), no censored solves in
those cases. Several synthetic traces are censored: a native one-step solve has
no second Ritz value. They are not assigned a fictitious stopping point. These
are fixed-trajectory diagnostics, not actual SCF speedups. The implementation
changes correction selection too, so its actual work must be measured afresh.
Evidence: `tests/reference/perf/p17a/replay-calibration.json`.

E1 scans 5e-8, 5e-9 and 5e-10 Ha. Stored relayed QE si8 conv_thr=1e-8 finishes
with estimated accuracy **5.4e-9 Ry = 2.7e-9 Ha**; k444s is approximately
**6e-9 Ry = 3e-9 Ha**. Source: `tests/reference/perf/qe_iteration_data_relayed.md`.
The architect reports the 1e-8/1e-12 energies equal at the printed digits. That
bounds their resolvable difference, not an exact zero error. Eight-decimal Ry
energies have a difference-rounding bound of about 5e-9 Ha.

**Frozen default: 5e-9 Ha**, recorded in `default-freeze.json` before the final
comparison. At this candidate, six reference-backed systems have delivered
energy errors from **1.3328e-10 to 1.0009e-9 Ha**. This lies within the stored QE
estimate scale and is below its ordinary printed-energy resolution. The looser
5e-8 candidate gives Si8 error **1.7910e-8 Ha**, exceeding that scale; 5e-10 adds
work without a requirement from this QE evidence. The decision uses MgO, bulk
Si, an elongated synthetic cell and a 16-atom cell as well as the silicon
performance fixtures. It does not tune the default on silicon alone or promise
that all systems have an error bound of 5e-9 Ha. Eigenvalue errors are visibly
larger than under strict and are reported, not hidden by the energy comparison.

## Actual final comparison

Native Mac, ordinary release, T=1. Times are process-wall observations from
serial runs, not dev08 acceptance or a statistically resolved speedup claim.
Commands, hashes, three thresholds, complete eigenvalue/density differences,
failures and times are in `tests/reference/perf/p17a/final-report.json` and
`final-table.md`. Strict counts below include certificate work; standard counts
include any first-iteration redo and no certificate.

| Case | strict SCF / Davidson / H vectors | standard 5e-9 SCF / Davidson / H vectors | standard dE (Ha) | max occupied+first-empty dEigen (Ha) |
| --- | --- | --- | ---: | ---: |
| Si8 | 15 / 345 / 2265 | 14 / 69 / 742 | 1.00089e-9 | 1.65593e-6 |
| k444s | 10 / 1827 / 11181 | 9 / 866 / 5864 | 1.33280e-10 | 1.56343e-6 |
| P13D Si, 20 Ry | 15 / 2865 / 16526 | 11 / 1042 / 6974 | 2.59334e-10 | 2.17985e-6 |
| P13D MgO, 40 Ry | 18 / 1306 / 12202 | 13 / 609 / 5909 | 1.50877e-10 | 9.93483e-6 |
| elongated 4 | 9 / 31 / 171 | 9 / 28 / 172 | 9.03338e-10 | 1.53106e-6 |
| Gamma 16 | 10 / 69 / 593 | 7 / 24 / 335 | 3.66114e-10 | 6.35428e-7 |

Elongated 4 is a negative work result: **one more H vector**, despite three fewer
Davidson steps. No universal speedup is claimed.

A tight reference requires two strict certified runs with successively tightened
`(energy, density, eigen)` settings `(1e-12,1e-9,1e-10)` and
`(1e-13,1e-10,1e-11)` in existing CLI units, with energy difference <=1e-12 Ha.
The tighter run supplies the reference. This empirical energy plateau is not a
ground-state or completeness proof. All physics inputs and UPFs remain fixed.
Density L2 is `sqrt(integral(delta_n^2))` in Bohr^-3/2 over the full grid;
Hartree norm is `sqrt(E_H[delta_n])` in sqrt(Ha), over the density sphere. Both
are unnormalized by electron count. The harness independently transforms the
stored real-space density using the existing stdlib separable DFT reader.

Three stress fixtures have **no reference**, so all their accuracy fields are
null and comparisons are work-only:

| Case | strict default | standard 5e-9 |
| --- | --- | --- |
| skew 20 | 6 SCF / 15 Davidson / 255 H; tight runs diverge | 100 SCF / 212 Davidson / 4168 H, exhausted |
| slab 4 | 100 / 539 / 2506, exhausted | 13 / 57 / 297, density-growth divergence |
| mesh 16 | 100 / 622 / 5973, exhausted | 100 / 442 / 8182, exhausted |

Skew's old strict convergence is already known to be spectrally incomplete
(P15C). It is not used as a ground-state oracle. Slab diverges sooner under
standard, and mesh16 needs more H vectors despite fewer Davidson steps. These
are negative results, not improvements. The full harness exits **3** because
these references are unavailable; the verification runner explicitly expects
that exit while requiring every validation test to pass.

E4 (strict with only Hartree metric): SCF counts are unchanged on the six usable
reference cases. H vectors change Si8 2265->2190, k444s 11181->11206,
Si 16526->16561, MgO 12202->12382, elongated 171->182, Gamma16 593->565.
MgO energy differs from strict/L2 by **-1.40687e-12 Ha**, slightly outside the
1e-12 Ha equivalence target. Thus **L2 remains the default**. Hartree remains an
explicit experimental switch, not a uniformly accepted optimization. Mesh16
also diverges under strict/Hartree; its error and available counters are retained.

## Reproduction and existing QE files

```
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --release -p dft-cli
python3 -I tests/p17a_replay_calibrate.py
DFT_NC_DIR=/absolute/path/to/pseudodojo/ex python3 -I tests/p17a_verify.py \
  --work-dir /tmp/p17a-final-verification
```

Use a fresh work directory. The verification runner runs the full comparison,
five-case strict identity checks at T=1 and T=8, the standard thread matrix,
CLI/checkpoint tests, Python tests, `scripts/check.sh` and README freshness.
`final-verification.json` records real exits and source/binary hashes. Synthetic
UPFs stay temporary; material UPFs are checked against the P13D manifest via
`DFT_NC_DIR`. Stored input strings reference temporary generated UPF paths;
regenerate them using the runner instead of reusing a stale temporary path.

The strict checks compare the feature build to the stored **pre-feature**
baselines: all five cases at T=1 and T=8 preserve energy exactly (0 Ha difference)
and density hashes. Standard/L2 and standard/Hartree across all nine cases at
T=1/T=8 preserve scientific histories, eigenvalues, counters, output density
hashes and failure statuses exactly: **18 pairs / 36 runs**. Timers are excluded.
Checkpoint tests cover both metrics, exact density/orbital byte identity,
profile mismatch rejection and a warm-start first redo. Solver tests cover the
missing-first-change rule and mandatory root prefix; analytic metric tests
cover zero mode, sphere cutoff, sign, units and the weighted Gram formula.

Explicit stored QE inputs require no new remote run:

```
python3 -I scripts/equal_accuracy.py --cases '' \
  --qe-output tests/reference/kpoints/k222_c/qe.out \
  --qe-output tests/reference/kpoints/k222s_c/qe.out \
  --qe-output tests/reference/si30_auto/qe.out \
  --qe-output tests/reference/si30_ours/qe.out --output-dir /tmp/p17a-qe
```

This parses all four real stored outputs, records energy, eigenvalues, iterations,
error estimates, exact timer counts where present, hashes and internal WALL.
Without a supplied reference it deliberately leaves delivered accuracy null
(exit 3). `--qe-reference` adds a same-model output reference, checking printed
geometry, PP MD5 identity, cutoffs, bands, XC and k order. The existing si30
ours/auto pair tests this interface, but is **not** a 1e-8/1e-12 convergence pair.
A manifest interface additionally binds full inputs, commands and UPF SHA-256.
Optional `--qe-density`/`--qe-reference-density` read supported complex-Gamma
`charge-density.dat`; geometry is taken from rounded stdout, so precision is
limited. Existing Gamma-only half-sphere files remain explicitly unsupported.

QE Davidson steps use initial solver calls plus expansion timer calls where
available, matching the inclusion of BAANDI's initial reduced solve. Raw H block
calls are reported separately; **H-vector counts remain null**, because block
width is not recoverable from these outputs. Rounded average iterations never
manufacture exact work. QE WALL excludes launcher overhead and is rounded;
BAANDI wall includes process start. No cross-code wall ratio or new QE run is
claimed. Error estimates are not substituted for delivered energy differences.

`strict-preparation-*`, `profile-readiness-*`, `verification.json` and
`qe-inventory.json` retain the initial partial run as historical evidence; use
`final-*`, `default-freeze.json` and `thread-identity.json` for current status.
README generation remains fresh with no generated diff; this ticket adds no
new cross-code validation or accepted dev08 performance claim. Dashboard/session
uploads remain the orchestrator's workflow; this worktree's registry had no P17A
running entry and no upload destination, so none was fabricated.
