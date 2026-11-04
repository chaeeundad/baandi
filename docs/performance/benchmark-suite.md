# P10G generality and robustness suite

This suite establishes a broader measurement surface; it changes no Rust code,
physical model, SCF defaults, tolerances or accepted optimization. See
[usage and case schema](../../benchmarks/README.md) and
[the machine-readable baseline](suite-baseline.json).

## Design and limits

- The runner and UPF generator use only Python's standard library, including on
  this Mac's Python 3.9.6. Case files are a documented restricted TOML subset.
- Actual executable SHA-256, repository commit and dirty state, platform, CPU
  count, Python, pinned toolchain and Cargo.lock hash are recorded. Successful
  CLI output additionally reports its own build commit/compiler and threading.
  Repository HEAD alone is not proof of executable provenance. Build with the
  documented release command; the suite does not silently rebuild a binary.
- Timings are subprocess wall seconds, include process setup and output I/O,
  exclude building and generating inputs, and retain every repeat plus median,
  minimum and maximum. Engine and phase times are separately retained. Shared
  Mac timings are indicative, not accepted speedup evidence; quiet dev08 runs
  with repeated measurements are needed before adopting performance changes.
- Operator applications count Hamiltonian vector columns: sum every history
  iteration and its certificate, including unsuccessful certification attempts.
  Never add the duplicated `last_iteration` field. When the CLI fails before
  writing results these counts, energy and iterations are null, not zero.
- All results have explicit status and convergence flags. Unconverged energies
  are diagnostics, never converged references. Timeout kills and reaps the direct
  CLI process; this CLI does not launch descendant processes. Complete stdout
  and stderr remain in the work directory. Timeouts lack final engine counters.
- `--robustness` continues after numerical errors, budget exhaustion, invalid
  output and timeouts and returns zero when recording finishes. Without it the
  first unexpected failure stops further cases with exit 3, preserving results.
  An exact expected-error substring only matches solver errors, not timeouts.
  A repaired expected failure automatically becomes a pass and is identified as
  an improvement. Setup failures are exit 1. Baseline regressions are exit 5,
  including in robustness mode. Exit zero alone is not evidence of convergence.
- Baseline comparisons cover only selected cases. A missing selected entry or
  different input/UPF fingerprint is a comparison failure. The portable input
  fingerprint excludes temporary paths but includes seed, threads, bands,
  tolerances and the potential bytes. Different configurations need their own
  baseline. Baseline files are never updated during comparison.
- For deterministic passing cases, absolute energy change must be <=1e-12 Ha;
  iteration count and convergence/status must agree. `--require-identical`
  additionally checks the density SHA-256. Each run records hashes even without
  that option. Repeats must agree exactly on status, energy, density and iteration
  count. Failure-to-pass is an improvement; changed error text is reported.
  Timing and operator-count differences are evidence, not arbitrary hard gates.
- Inputs include real Si, cubic/skewed/elongated/slab geometries, 1–16 atoms,
  Gamma/even/odd/shifted meshes, low/high cutoffs, s/p/d projectors, valence 1–14,
  and a 100-band stress case. Synthetic fixtures use the in-repository test
  conventions: exponential radial mesh, local/DIJ in Ry, beta=r*projector,
  same-channel symmetric DIJ coupling, zero atomic starting density. Neither
  their NC label nor numerical convergence establishes physical correctness.
- The UPF reader supports **three total projectors**. The suite covers one to
  three in a single channel and one each in s/p/d; nine-projector s3/p3/d3 is
  unsupported and deliberately not presented as covered. Broad physical
  materials validation still needs independently validated real potentials.
- Automatic FFT selection rounds to 2/3/5-smooth dimensions. The supported
  explicit override exercises 7x11x13 instead; the CLI independently validates
  that the grid is large enough. The slab has an outer-layer separation of
  8 Bohr and a remaining periodic gap of 12 Angstrom. It is a numerical probe,
  not a validated surface calculation or dipole-corrected model.
- Cubic cases expose symmetry-degenerate excited multiplets; the suite does not
  assert all eigenvalues remain exactly equal in floating arithmetic. Direct
  and indirect occupied/empty gaps are measured in Ha from reported occupations
  and eigenvalues. A proposed small-gap geometry may be rejected by the engine's
  insulating-gap certificate; the result is retained as a robustness finding.
  No claim of general metallic support or smearing is made.
- Random fractional displacement streams are reproducible per case and seed.
  Every failure includes a standalone exact input and potential hash/parameters.
  The suite does not automatically shrink a structure to a mathematically
  minimal failure or invent missing partial-solve counters.

## Verification

The recorded measurement command, results and findings follow below. No QE run
was made for this ticket and no synthetic energy is described as a QE reference.
Existing repository QE regressions are exercised by the Python test suite.
Only owned files are changed; the ticket explicitly excludes artifacts, tickets,
root README, manifests and crates. No dashboard upload or root README generation
is performed under that restriction. Changes remain uncommitted for review.

## First Mac baseline and findings (2026-11-03)

Host: macOS 27.0.1 arm64, 12 logical CPUs, Python 3.9.6; release Rust 1.99.1,
commit `0142ebf639e90a33c27b945a0e4c61d4c030b538`, requested thread budget 1.
The executable and generator hashes are in JSON. Hardware-model sysctl access
was denied by the sandbox; no exact CPU product name is asserted.

```sh
python3 -I scripts/bench_suite.py \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --robustness --include-long --variants 1 --repeats 2 \
  --output /tmp/p10g-baseline-final.json
```

Exit 0: 25 base cases plus 25 displaced variants, each executed twice (100 CLI
runs). All repeats agree exactly on status, energy, density hash and iterations.
34 inputs pass, 9 exhaust the SCF budget, and 7 return errors. Among the 25 base
cases, 17 pass, 4 exhaust SCF, and 4 error. These failures are explicit findings,
not a claim that all systems converge. The 19 default base cases take a sum of
8.629413 median seconds; the slowest individual default run is 5.812720 seconds.
Six opt-in cases remain conservatively classified long because their intended
workloads are larger, even when a current early failure makes them cheap.
Shared-host variability is visible in the stored repeat samples.

Numerical evidence:

- Si Gamma: -7.295530468210772 Ha, 10 SCF iterations, 422 Hamiltonian vector
  applications including certificates. Last main-solve eigen residual is
  3.904607422511523e-9 Ha; relative charge error is 1.3322676295501878e-15.
- Passing small-gap Si skew: 0.0031446358589312107 Ha occupied/empty gap;
  passing synthetic s/p/d probe: 0.7924507025384524 Ha. These establish different
  numerical gap regimes, not physical validation of the synthetic system.
- Across 34 passing inputs, maximum reported main-solve eigen residual is
  8.574732451721064e-9 Ha and maximum relative charge error is
  2.6645352591003757e-15. Eigenvalue arrays and certificate diagnostics are stored.
  The CLI does not expose orthogonality errors, so none are invented here;
  existing solver/engine acceptance tests cover that invariant.
- Explicit FFT probe ran with exactly [7, 11, 13]; base energy
  -0.7346843876893239 Ha, 12 iterations. This is a synthetic discretization
  experiment, not a cutoff-converged material result.

| Finding | Reproducer IDs | Observed behavior and interpretation |
|---|---|---|
| P7PC symmetry-sensitive Davidson stagnation | si8-gamma | 400 iterations, residual 4.238824601895987e-6 Ha; matches the ticket diagnosis exactly. Default band settings retained. |
| Larger real-cell eigensolver failures | si8-k222; si16-gamma and its variant | Residuals 8.064599029131822e-6, 1.828694229425743e-2 and 1.2279455730020664e-2 Ha. No final energy/counters are available. |
| SCF robustness limits | si-elongated and variant; si-slab and variant; si-skew-variant-0 | All reach 100 SCF iterations without convergence. Slab is retained as a failing experimental case; no successful surface calculation is claimed. |
| Insulating-only occupation limitation | si-k222-variant-0 | `NonInsulating`, gap -0.01776592622951892 Ha versus required 1e-7 Ha. This is a model-boundary rejection, not Davidson stagnation. |
| Multi-projector/high-valence SCF stress | syn-s3 and variant; syn-z14 and variant | Reach 100 iterations. Synthetic energies are diagnostics only. |
| Large-band CLI ceiling | syn-100-bands and variant | `bands/subspace incompatible with k basis`. Engine default Davidson subspace is 48; CLI exposes nbands but no subspace override. The 100-band input is retained so this limitation can become a passing regression after a separately owned fix. No 100-band SCF performance coverage is claimed yet. |

The cubic synthetic probe is nominally symmetric, but finite grids and partial
shell occupation can split its excited levels. The real cubic Si8 case is the
maintained top-multiplet stagnation reproducer. Displacing Si8 breaks the
symmetry and converges in both baseline Gamma and mesh variants; this is
observational evidence of sensitivity, not a new causal proof or a fix.

Additional seed/thread validation is embedded in JSON under
`additional_validation.seed7_threads4`:

```sh
python3 -I scripts/bench_suite.py --cases si-gamma,si8-gamma,syn-spd \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --robustness --seed 7 --threads 4 --variants 2 --repeats 2 \
  --output /tmp/p10g-seed7.json
```

Exit 0: nine inputs, eight pass and the unperturbed Si8 retains its known error;
all two-run status/energy/density/iteration identities hold. Baseline fingerprints
intentionally distinguish seed/thread configurations. This is not a thread
scaling measurement.

Verification completed:

- `cargo build --offline --release -p dft-cli`: exit 0.
- `scripts/check.sh`: exit 0; build/test/fmt/clippy all PASS.
- `python3 -I -m unittest discover -s scripts/tests`: exit 0; 126 tests,
  including 13 new suite tests. Fake CLI tests exercise timeouts, malformed
  output, missing output, nonfinite energy, failures, nonconvergence,
  nondeterminism, portability of fingerprints and regression detection.
  Independent XML/analytic checks validate generated UPF mesh, units, channels,
  same-channel coupling and array lengths; real suite runs exercise Rust parsing.
- `cargo test --offline --release --workspace -- --ignored`: exit 0;
  17 ignored release acceptance tests passed.
- `scripts/bench_scf.py` with gamma,k222s,k222,k444s, one repeat,
  `--baseline docs/performance/baseline-after-o2.json --require-identical`:
  exit 0; all four energy differences are exactly 0 Ha and density hashes match.
  Reported wall-time ratios are not interpreted as speedups.


- Full suite re-comparison, 50 inputs, one fresh repeat, same seed/thread budget:
  `--robustness --include-long --variants 1 --baseline docs/performance/suite-baseline.json --require-identical`
  (same binary/UPF, output `/tmp/p10g-recompare.json`): exit 0, no regressions;
  all 34 passing energies differ by exactly 0 Ha and all density hashes match.
  Failure statuses/messages also match. This additionally verifies comparison
  across different repeat counts (two baseline repeats, one comparison repeat).
  Verification commands, exit codes and comparison records are embedded in JSON.

The complete table and failing single-case inputs below were generated by the
runner from the first baseline. Synthetic generation parameters and hashes are
in each corresponding JSON entry; no UPF bytes are included in either report.

# Benchmark suite results

Synthetic results are non-physical. Timings are indicative on a shared host.

| Case | Class | Status | Seconds (median) | Iterations | Energy (Ha) |
|---|---|---|---:|---:|---:|
| si-elongated | elongated | nonconverged | 0.771 | 100 | -15.087271907813104 |
| si-elongated-variant-0 | elongated | nonconverged | 0.805 | 100 | -14.391129151658008 |
| si-gamma | si-primitive | pass | 0.152 | 10 | -7.295530468210772 |
| si-gamma-variant-0 | si-primitive | pass | 0.167 | 11 | -7.294564148045221 |
| si-high | si-cutoff | pass | 0.434 | 11 | -7.296038607887432 |
| si-high-variant-0 | si-cutoff | pass | 0.450 | 11 | -7.289858372720575 |
| si-k222 | si-primitive | pass | 1.000 | 10 | -7.831230057121145 |
| si-k222-variant-0 | si-primitive | error | 1.062 | None | None |
| si-k222s | si-primitive | pass | 0.542 | 10 | -7.919423412469067 |
| si-k222s-variant-0 | si-primitive | pass | 0.614 | 10 | -7.919082857019964 |
| si-k444s | si-primitive | pass | 3.878 | 10 | -7.925976126648081 |
| si-k444s-variant-0 | si-primitive | pass | 3.868 | 10 | -7.925412930775879 |
| si-low | si-cutoff | pass | 0.021 | 11 | -7.233756452824805 |
| si-low-variant-0 | si-cutoff | pass | 0.022 | 11 | -7.232285425784414 |
| si-skew | skewed | pass | 0.182 | 13 | -7.250555672960804 |
| si-skew-variant-0 | skewed | nonconverged | 0.830 | 100 | -7.170684309735211 |
| si-slab | slab | nonconverged | 1.029 | 100 | -15.23143886907433 |
| si-slab-variant-0 | slab | nonconverged | 1.022 | 100 | -15.235565931334294 |
| si16-gamma | si-supercell | error | 11.352 | None | None |
| si16-gamma-variant-0 | si-supercell | error | 13.836 | None | None |
| si8-gamma | si-conventional | error | 5.792 | None | None |
| si8-gamma-variant-0 | si-conventional | pass | 7.001 | 22 | -31.31645520439176 |
| si8-k222 | si-conventional | error | 10.673 | None | None |
| si8-k222-variant-0 | si-conventional | pass | 23.267 | 12 | -31.663211397321028 |
| syn-100-bands | large-bands | error | 0.230 | None | None |
| syn-100-bands-variant-0 | large-bands | error | 0.230 | None | None |
| syn-cubic-degenerate | synthetic-gap | pass | 0.009 | 10 | -4.261277234135343 |
| syn-cubic-degenerate-variant-0 | synthetic-gap | pass | 0.009 | 10 | -4.261277340138134 |
| syn-d3 | synthetic-projectors | pass | 0.005 | 8 | -0.9465076331236224 |
| syn-d3-variant-0 | synthetic-projectors | pass | 0.005 | 8 | -0.9465084028570097 |
| syn-fft-primes | awkward-fft | pass | 0.007 | 12 | -0.7346843876893239 |
| syn-fft-primes-variant-0 | awkward-fft | pass | 0.008 | 16 | -0.7346843893082529 |
| syn-gaussian-fcc | synthetic-gap | pass | 0.052 | 16 | -4.683023405363695 |
| syn-gaussian-fcc-variant-0 | synthetic-gap | pass | 0.047 | 15 | -4.682616704438532 |
| syn-odd-shift | synthetic-mesh | pass | 0.010 | 9 | -1.5900664019906925 |
| syn-odd-shift-variant-0 | synthetic-mesh | pass | 0.009 | 9 | -1.5900664020054043 |
| syn-p3 | synthetic-projectors | pass | 0.005 | 8 | -0.9465076597173665 |
| syn-p3-variant-0 | synthetic-projectors | pass | 0.005 | 10 | -0.9465086826168907 |
| syn-s1 | synthetic-projectors | pass | 0.005 | 13 | -0.7346754280328645 |
| syn-s1-variant-0 | synthetic-projectors | pass | 0.004 | 13 | -0.7346765225858594 |
| syn-s2 | synthetic-projectors | pass | 0.005 | 17 | 0.63040952833334 |
| syn-s2-variant-0 | synthetic-projectors | pass | 0.005 | 15 | 0.6306518063739004 |
| syn-s3 | synthetic-projectors | nonconverged | 0.009 | 100 | 1.8287903113034467 |
| syn-s3-variant-0 | synthetic-projectors | nonconverged | 0.009 | 100 | 1.7654427206126853 |
| syn-spd | synthetic-projectors | pass | 0.005 | 11 | -0.7346751923893953 |
| syn-spd-variant-0 | synthetic-projectors | pass | 0.005 | 11 | -0.7346742428850579 |
| syn-z1 | synthetic-valence | pass | 0.008 | 7 | -0.9291990923944967 |
| syn-z1-variant-0 | synthetic-valence | pass | 0.008 | 7 | -0.9291645936438467 |
| syn-z14 | synthetic-valence | nonconverged | 0.050 | 100 | -71.1401274918061 |
| syn-z14-variant-0 | synthetic-valence | nonconverged | 0.050 | 100 | -71.14685630604582 |

## By system class

| Class | Cases | Passed | Sum median seconds |
|---|---:|---:|---:|
| awkward-fft | 2 | 2 | 0.015 |
| elongated | 2 | 0 | 1.576 |
| large-bands | 2 | 0 | 0.460 |
| si-conventional | 4 | 2 | 46.734 |
| si-cutoff | 4 | 4 | 0.926 |
| si-primitive | 8 | 7 | 11.283 |
| si-supercell | 2 | 0 | 25.188 |
| skewed | 2 | 1 | 1.012 |
| slab | 2 | 0 | 2.051 |
| synthetic-gap | 4 | 4 | 0.117 |
| synthetic-mesh | 2 | 2 | 0.019 |
| synthetic-projectors | 12 | 10 | 0.065 |
| synthetic-valence | 4 | 2 | 0.117 |

## Baseline comparison

Only selected cases are compared; missing selected entries fail comparison.

## Findings and exact reproducing inputs

Replace SI_UPF_PATH with the read-only validation file. For synthetic cases regenerate synthetic.UPF using the case parameters in JSON. Runtime directories contain executable inputs and logs.

### si-elongated

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-elongated"
lattice_unit = "bohr"
lattice_vectors = [[8, 0, 0], [0, 8, 0], [0, 0, 24]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.125]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.625]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si-elongated-variant-0

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-elongated-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[8, 0, 0], [0, 8, 0], [0, 0, 24]]
[[system.atoms]]
species = "Si"
fractional = [0.9975255828034126, 0.9990397931302788, 0.9994659893370663]
[[system.atoms]]
species = "Si"
fractional = [0.2523548849308699, 0.2545664059075983, 0.13281813341733814]
[[system.atoms]]
species = "Si"
fractional = [0.00025473072909386156, 0.003799132844612179, 0.49207826948123207]
[[system.atoms]]
species = "Si"
fractional = [0.253269405536378, 0.2554272321123899, 0.6301576477291659]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si-k222-variant-0

Expected failure matched: False

```text
error: NonInsulating { gap_hartree: -0.01776592622951892, required_hartree: 1e-7 }
error: NonInsulating { gap_hartree: -0.01776592622951892, required_hartree: 1e-7 }
```

```toml
[system]
name = "si-k222-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[0, 5.130606428358966, 5.130606428358966], [5.130606428358966, 0, 5.130606428358966], [5.130606428358966, 5.130606428358966, 0]]
[[system.atoms]]
species = "Si"
fractional = [0.9901791708512576, 0.9909263712779467, 0.0002771014629703479]
[[system.atoms]]
species = "Si"
fractional = [0.258425442427287, 0.25650160843613207, 0.2556894057370224]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kgrid = [2, 2, 2]
kshift = [0, 0, 0]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si-skew-variant-0

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-skew-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[0.7, 5.130606428358966, 5.130606428358966], [5.130606428358966, 0, 5.130606428358966], [5.130606428358966, 5.130606428358966, 0.3]]
[[system.atoms]]
species = "Si"
fractional = [0.006406135754818533, 0.9965865765337489, 0.00582652923618567]
[[system.atoms]]
species = "Si"
fractional = [0.24085008964721813, 0.258193477722265, 0.2554338479773732]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si-slab

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-slab"
lattice_unit = "bohr"
lattice_vectors = [[8, 0, 0], [0, 8, 0], [0, 0, 30.67671349550924]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0.065196032172507]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0.195588096517521]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0.260784128690028]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si-slab-variant-0

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-slab-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[8, 0, 0], [0, 8, 0], [0, 0, 30.67671349550924]]
[[system.atoms]]
species = "Si"
fractional = [0.9903558171159005, 0.00481854395299695, 0.9935601642756073]
[[system.atoms]]
species = "Si"
fractional = [0.49239479052909824, 0.5046855669476572, 0.06625351773791652]
[[system.atoms]]
species = "Si"
fractional = [0.9910689754523436, 0.9956268286636651, 0.19401598590446395]
[[system.atoms]]
species = "Si"
fractional = [0.498453158533235, 0.5038734662221545, 0.25517467541270844]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si16-gamma

Expected failure matched: False

```text
error: Numerical("k index 0: not converged after 400 iterations; max residual 1.828694229425743e-2")
error: Numerical("k index 0: not converged after 400 iterations; max residual 1.828694229425743e-2")
```

```toml
[system]
name = "si16-gamma"
lattice_unit = "bohr"
lattice_vectors = [[10.261212856717933, 0, 0], [0, 10.261212856717933, 0], [0, 0, 20.522425713435865]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0.0]
[[system.atoms]]
species = "Si"
fractional = [0, 0.5, 0.25]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0, 0.25]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0.0]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.125]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.75, 0.375]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.25, 0.375]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.75, 0.125]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0, 0.5, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.625]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.75, 0.875]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.25, 0.875]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.75, 0.625]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si16-gamma-variant-0

Expected failure matched: False

```text
error: Numerical("k index 0: not converged after 400 iterations; max residual 1.2279455730020664e-2")
error: Numerical("k index 0: not converged after 400 iterations; max residual 1.2279455730020664e-2")
```

```toml
[system]
name = "si16-gamma-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[10.261212856717933, 0, 0], [0, 10.261212856717933, 0], [0, 0, 20.522425713435865]]
[[system.atoms]]
species = "Si"
fractional = [0.996891782293039, 0.9922205129524883, 0.9955197848416522]
[[system.atoms]]
species = "Si"
fractional = [0.9910069849059968, 0.5046099819984217, 0.25019658660589056]
[[system.atoms]]
species = "Si"
fractional = [0.4906593039631164, 0.9981744418164215, 0.24914199367422213]
[[system.atoms]]
species = "Si"
fractional = [0.5074565697459381, 0.5094036413744051, 0.007190239001362985]
[[system.atoms]]
species = "Si"
fractional = [0.25979918161596555, 0.24008828561534662, 0.12586820972588308]
[[system.atoms]]
species = "Si"
fractional = [0.24606129093697268, 0.7564788556060342, 0.3678321908151959]
[[system.atoms]]
species = "Si"
fractional = [0.7431190054413599, 0.2505653205046344, 0.3762617103339205]
[[system.atoms]]
species = "Si"
fractional = [0.7445407274279814, 0.7553515984402055, 0.13060040798169184]
[[system.atoms]]
species = "Si"
fractional = [0.0007410761818156461, 0.9955397392412683, 0.5025721787732014]
[[system.atoms]]
species = "Si"
fractional = [0.005643671717503568, 0.4952915421445436, 0.7562702758250317]
[[system.atoms]]
species = "Si"
fractional = [0.4984144961391635, 0.9937950577355975, 0.7577229155098875]
[[system.atoms]]
species = "Si"
fractional = [0.49417011233161956, 0.5025010045202687, 0.5046121924202047]
[[system.atoms]]
species = "Si"
fractional = [0.25140403716671705, 0.25424973320077066, 0.6279743190935048]
[[system.atoms]]
species = "Si"
fractional = [0.2561430372783232, 0.7402943494416117, 0.8682473186966103]
[[system.atoms]]
species = "Si"
fractional = [0.7442981500767792, 0.2515186853703962, 0.8832171790831813]
[[system.atoms]]
species = "Si"
fractional = [0.7519768324513582, 0.7564331585090319, 0.6181936201198243]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si8-gamma

Expected failure matched: True

```text
error: Numerical("k index 0: not converged after 400 iterations; max residual 4.238824601895987e-6")
error: Numerical("k index 0: not converged after 400 iterations; max residual 4.238824601895987e-6")
```

```toml
[system]
name = "si8-gamma"
lattice_unit = "bohr"
lattice_vectors = [[10.261212856717933, 0, 0], [0, 10.261212856717933, 0], [0, 0, 10.261212856717933]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0]
[[system.atoms]]
species = "Si"
fractional = [0, 0.5, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.25]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.75, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.25, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.75, 0.25]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### si8-k222

Expected failure matched: False

```text
error: Numerical("k index 0: not converged after 400 iterations; max residual 8.064599029131822e-6")
error: Numerical("k index 0: not converged after 400 iterations; max residual 8.064599029131822e-6")
```

```toml
[system]
name = "si8-k222"
lattice_unit = "bohr"
lattice_vectors = [[10.261212856717933, 0, 0], [0, 10.261212856717933, 0], [0, 0, 10.261212856717933]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0]
[[system.atoms]]
species = "Si"
fractional = [0, 0.5, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.5, 0.5, 0]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.25]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.75, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.25, 0.75]
[[system.atoms]]
species = "Si"
fractional = [0.75, 0.75, 0.25]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 30
kgrid = [2, 2, 2]
kshift = [0, 0, 0]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### syn-100-bands

Expected failure matched: False

```text
error: Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")")
error: Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")")
```

```toml
[system]
name = "syn-100-bands"
lattice_unit = "bohr"
lattice_vectors = [[16, 0, 0], [0, 16, 0], [0, 0, 16]]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.1, 0.075]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.1, 0.325]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.1, 0.575]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.1, 0.825]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.6, 0.075]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.6, 0.325]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.6, 0.575]
[[system.atoms]]
species = "Xx"
fractional = [0.05, 0.6, 0.825]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.1, 0.075]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.1, 0.325]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.1, 0.575]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.1, 0.825]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.6, 0.075]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.6, 0.325]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.6, 0.575]
[[system.atoms]]
species = "Xx"
fractional = [0.55, 0.6, 0.825]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
nbands = 100
```

### syn-100-bands-variant-0

Expected failure matched: False

```text
error: Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")")
error: Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")")
```

```toml
[system]
name = "syn-100-bands-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[16, 0, 0], [0, 16, 0], [0, 0, 16]]
[[system.atoms]]
species = "Xx"
fractional = [0.05969137832563314, 0.10838302283431736, 0.06592249980704892]
[[system.atoms]]
species = "Xx"
fractional = [0.05444731401498902, 0.09294838517419952, 0.33152017325138217]
[[system.atoms]]
species = "Xx"
fractional = [0.059781256086911656, 0.10653208216322728, 0.5822675286486751]
[[system.atoms]]
species = "Xx"
fractional = [0.0471515168029032, 0.10657572163914163, 0.8305819784576626]
[[system.atoms]]
species = "Xx"
fractional = [0.05086739506838213, 0.6018486115588041, 0.07996367925363437]
[[system.atoms]]
species = "Xx"
fractional = [0.04537623456658771, 0.5901191022326344, 0.3164011041398362]
[[system.atoms]]
species = "Xx"
fractional = [0.05548594562061217, 0.5934561079842774, 0.5802047882368947]
[[system.atoms]]
species = "Xx"
fractional = [0.049758629112830174, 0.5991880545065725, 0.8190112505002038]
[[system.atoms]]
species = "Xx"
fractional = [0.5585926190266389, 0.09368425617057158, 0.07017627525976962]
[[system.atoms]]
species = "Xx"
fractional = [0.5487993287715339, 0.09692755088513202, 0.33085720467528895]
[[system.atoms]]
species = "Xx"
fractional = [0.5492188759877578, 0.09246706968317847, 0.5751793519055348]
[[system.atoms]]
species = "Xx"
fractional = [0.5551944405757998, 0.10924694790066015, 0.8262287080880256]
[[system.atoms]]
species = "Xx"
fractional = [0.5548979847646993, 0.6023427065689841, 0.07760494995531664]
[[system.atoms]]
species = "Xx"
fractional = [0.5471298766597343, 0.602495635169312, 0.31828269810587706]
[[system.atoms]]
species = "Xx"
fractional = [0.5446626947091955, 0.5993037698997039, 0.5825052413475285]
[[system.atoms]]
species = "Xx"
fractional = [0.5549061894320149, 0.591892317877482, 0.8196319610571026]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
nbands = 100
```

### syn-s3

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-s3"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Xx"
fractional = [0.13, 0.21, 0.32]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 3
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### syn-s3-variant-0

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-s3-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Xx"
fractional = [0.13875838037215502, 0.21831426196286577, 0.3194596330380067]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 3
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### syn-z14

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-z14"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Xx"
fractional = [0.13, 0.21, 0.32]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 12
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
nbands = 12
```

### syn-z14-variant-0

Expected failure matched: False

```text
error: NOT CONVERGED: iteration budget exhausted
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-z14-variant-0"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Xx"
fractional = [0.12383112564778931, 0.20577802353128208, 0.31715227475418123]
[pseudopotentials]
Xx = "synthetic.UPF"
[basis]
ecutwfc_ry = 12
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
nbands = 12
```

## Notes added 2026-11-04 (dev08 run by the remote measurement session)
- Running the default tier without `--robustness` aborts with exit 3 at the first case that does not converge (`si-elongated`, `si-slab` and variants are known to hit the 100 iteration limit with the default linear/Pulay 0.5 settings, on both backends). Use `--robustness` to record them as findings and continue. These are genuine SCF stagnation cases and are inputs for ticket P11M.
- Synthetic cases failed the input preflight (ADR 0022): species label `Xx` is not an element symbol. The synthetic UPF is now written with `element="Si"` and the species label is `Si` (non-physical potential, label only satisfies the identity check). Input hashes of synthetic cases changed; the stored `suite-baseline.json` entries for them must be regenerated before using `--baseline` on those cases.
