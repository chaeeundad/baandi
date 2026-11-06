# Generality benchmark inputs

Run from the repository root with Python 3.9 or later (standard library only):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
python3 -I scripts/bench_suite.py \
  --upf /absolute/path/to/Si.pz-vbc.UPF \
  --robustness --output /tmp/suite.json
```

The default selection contains cases measured to finish in less than five minutes
**per case** on the shared Mac. A case may finish with a numerical failure.
The timeout is 240 seconds per default run and 900 seconds per long run.
`--include-long` enables all cases; `--cases si8-gamma,syn-spd` selects exactly
those IDs, regardless of tier. The release build is separate from timed runs.
No potential is stored here. `--upf` is read-only and required only when Si cases
are selected. Synthetic files are generated in the run directory and explicitly
labelled non-physical. They exercise the NC parser and operator interface; they
are not transferable atomic potentials or validated norm-conserving solutions.

Each run creates a JSON record and adjacent Markdown report. Inputs, generated
synthetic potentials, stdout/stderr and original CLI output remain in a fresh
temporary directory printed at completion. Use `--work-dir /tmp/new-run-dir`
to choose its location. Existing output/work paths are rejected. Remove your
own temporary run directory when no longer needed; real UPFs are never copied.

## Case format

`cases/*.toml` deliberately uses a flat TOML subset compatible with Python 3.9:
one `key = value` per line, with JSON-compatible strings, numbers, booleans and
arrays, and whole-line `#` comments. No tables, trailing commas, inline comments,
single-quoted strings, multiline values, dates or special floating values.
Unknown/duplicate keys are rejected. These files are templates consumed by the
suite, **not** CLI inputs. Generated `input.toml` files are full CLI inputs.

Required fields: `id`, `purpose`, `system_class`, `tier` (`default` or `long`),
`potential` (`si` or `synthetic`), `lattice` (three user lattice vectors in Bohr),
`atoms` (fractional three-vectors), `cutoff_ry`, `timeout_s`.
Optional fields: `kgrid`, `kshift`, `fft_grid`, `nbands`, `max_iterations`,
`deterministic` (default true), `expected_error`; synthetic controls `valence`,
`channels` (angular momentum per projector), `radius` (Bohr), `shape`.
SCF settings match the existing benchmark: Pulay beta 0.5, 100 iterations,
energy tolerance 1e-8 Ha and density residual tolerance 1e-7 Bohr^(-3/2).
No eigensolver tolerance, extra-band count or engine default is changed.

## Reproduction and comparison

```sh
# Deterministic repeated measurements and baseline comparison:
python3 -I scripts/bench_suite.py --upf /path/Si.pz-vbc.UPF \
  --robustness --repeats 3 --threads 1 --seed 1 \
  --baseline docs/performance/suite-baseline.json --require-identical \
  --output /tmp/compare.json

# Independent random displacement streams for each selected base case:
python3 -I scripts/bench_suite.py --cases si-gamma,si8-gamma,syn-spd \
  --upf /path/Si.pz-vbc.UPF --robustness --seed 7 --variants 3 \
  --displacement 0.01 --output /tmp/variants.json

# Regenerate a synthetic potential independently (write outside repository):
python3 -I scripts/gen_synthetic_upf.py /tmp/synthetic.UPF \
  --valence 14 --channels 0,1,2 --radius 0.5 --shape gaussian
```

Variants perturb each fractional coordinate uniformly within the stated
amplitude and wrap into [0,1). The stream depends on case ID, seed and variant
index, so changing selection/order does not change a variant. This is a
robustness search, not a thermodynamic ensemble, geometry optimizer or collision
filter. Variant failures carry their actual coordinates, not only a seed.

Failure reproducers in JSON and Markdown contain the complete single-case CLI
input. Replace `SI_UPF_PATH` with the original validation UPF path (verify the
stored SHA-256); regenerate `synthetic.UPF` with the embedded case parameters
and generator version/hash. Runtime `input.toml` files already contain absolute
paths. No automatic input minimization is claimed: each reproducer is one
failing calculation, without the rest of the suite.

For limitations, comparison semantics, first measurements and findings see
[the benchmark report](../docs/performance/benchmark-suite.md).

## Catalog

| ID | Tier | Purpose |
|---|---|---|
| si-elongated | default | Four atoms, aspect ratio 1:1:3; robustness probe |
| si-gamma | default | Existing primitive Si Gamma, 8 default bands; insulating reference |
| si-high | long | High cutoff and larger FFT/basis |
| si-k222 | default | Existing unshifted 2x2x2 mesh |
| si-k222s | default | Existing half-shifted 2x2x2 mesh |
| si-k444s | long | Existing shifted 4x4x4 mesh (32 retained k points) |
| si-low | default | Low cutoff sensitivity, not a converged physical cutoff |
| si-skew | default | Sheared primitive lattice, broken degeneracies |
| si-slab | long | Two bilayers; 12 Angstrom vacuum between outer layers across boundary |
| si16-gamma | long | 1x1x2 conventional supercell, 16 atoms |
| si8-gamma | default | P7PC top-band degeneracy/stagnation regression; default bands |
| si8-k222 | long | Conventional cell with mesh and default bands |
| syn-100-bands | long | 16 Z=12 ions, 96 occupied plus 4 buffer bands; large nonlocal workload |
| syn-cubic-degenerate | default | Cubic single ion; symmetry produces degenerate excited multiplets |
| syn-d3 | default | Three coupled projectors in d channel |
| syn-fft-primes | default | Explicit 7x11x13 grid: chooser otherwise rounds to 2/3/5 smooth sizes |
| syn-gaussian-fcc | default | Two ions, Gaussian local potential and two s/one p projectors; fixture-like gap |
| syn-odd-shift | default | Odd 3x1x1 mesh with half shift; sampled-gap rejection probe |
| syn-p3 | default | Three coupled projectors in p channel |
| syn-s1 | default | Small soft-Coulomb fixture with one s projector |
| syn-s2 | default | Two coupled s projectors |
| syn-s3 | default | Three coupled projectors in s channel |
| syn-spd | default | One projector each in s,p,d channels |
| syn-z1 | default | Two light-looking Z=1 ions, even total electron count |
| syn-z14 | default | Heavy-looking Z=14, narrow local core, 12 bands |

## P25A large and second-phase candidates

New opt-in `long` cases: `si32-gamma` (32 atoms, 2x2x1 conventional cells)
and `si16-k222` (16 atoms, 1x1x2 conventional cells, unshifted 2x2x2 mesh).
Each has `-forces`, `-stress`, and `-relax1` companions. The relaxation
companions displace the first atom by 0.001 along the first fractional axis.
All timings and convergence qualifications await dev08; these are inputs only.
A 2x2x2 conventional diamond cell would contain 64 atoms, so it is not used.

Optional `phase` accepts `scf` (omitted default), `forces`, `stress`, or `relax1`.
Old files, generated inputs, fingerprints and report layout stay unchanged when
phase is omitted. Explicit phases appear in JSON and the Markdown phase list.
`scf` uses the existing command; derivatives enable `[output] forces=true` or
`stress=true`. Their `runs[].phase_seconds` is the CLI's derivative timer;
`wall_*_s` still includes SCF, process startup and output. Stress retains the
CLI's Pulay assessment cost and must not be compared to a differently scoped
QE timer without accounting for that work.

`relax1` uses `dft-cli relax` with BFGS and `max_steps=1`: initial certified SCF
plus one ionic move and its certified SCF/forces. The CLI does not expose a
separate move timer, so `phase_seconds` is null; wall time is the full workload.
Exactly two frames, one move and no failures/retries are required for a passing
workload. Expected CLI exit 3 (ionic budget exhausted) is accepted in this case;
`converged` remains false. A pass does not claim geometry convergence. Relaxation
repeat identity uses the trajectory hash; density hashes remain null because
this command does not write density.bin. `--require-identical` compares the
trajectory for relaxation. SCF iteration totals sum the two frames.

Example (external read-only Si UPF, no QE execution):

```sh
python3 -I scripts/bench_suite.py --upf /path/Si.pz-vbc.UPF \
  --cases si32-gamma,si32-gamma-forces,si32-gamma-stress,si32-gamma-relax1 \
  --threads 8 --repeats 3 --robustness --output /tmp/large-phase.json
```
