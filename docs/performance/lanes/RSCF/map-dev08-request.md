# R-SCF-MAP: same-input output-density diagnostic

Issue #25 patch 2B. Source base: `5c3c9db49fbdb754e7f9bfcc76ddb2cebb8e4b0d`;
changes are uncommitted for review. This is a correctness experiment, with **no
speed comparison**. Defaults and the production SCF dispatcher are unchanged.

## Interface and interpretation

```sh
dft-cli refscf-map INPUT.toml --qe TRAJECTORY --steps 6 \
  --ethr-from-stdout --solver legacy
dft-cli refscf-map INPUT.toml --qe TRAJECTORY --steps 6 \
  --ethr-from-stdout --solver reference
```

`legacy` is the existing real-Gamma Davidson map, with its standard Ritz-change
stopping policy and no residual certification. `reference` calls the separately
available `ref_davidson`, using the existing Hamiltonian, atomic-start, density,
XC and energy kernels. The latter supplies the full plane-wave Hamiltonian
preconditioner diagonal (kinetic + local mean + nonlocal), shared by both ORC
coordinates without another Gamma weight. Its algorithm and stopping flags
remain those of the standalone solver; see [its contract](../../../design/reference-scf.md#inner-solver).
Neither arm performs mixing, verification or an SCF convergence decision.

At **every** k the input is QE's own density, matched by Miller triple with exact
membership/cell checks: `(mix1 - beta*out1)/(1-beta)` at k1, preceding QE `mix`
thereafter. Beta comes from the input reference profile (default 0.7); beta=1
cannot recover k1 and is rejected. NC, unpolarized, integer occupations, Gamma
and symmetry off are required. Density cutoff is 4 times the wavefunction
cutoff; other density spheres reject by membership. The input's geometry, UPF,
XC, grid and band count must match QE. In particular the Mac Si8 fixture has
**20 bands**; its old mixing-only TOML has no band count and must be augmented.

The orbital start is BAANDI atomic at k1, then the preceding BAANDI orbitals.
`--cold-start` uses fresh BAANDI atomic orbitals at every k. QE wavefunctions
are not available in this bundle. A finite inner threshold defines an
approximate map that depends on the starting subspace and solver; a density
mismatch alone does **not** establish an error in the potential. In particular,
the legacy and reference buffer-band policies need not be identical.

`--ethr-from-stdout` (also the default when `--ethr` is absent) reads the last
printed ethr in iteration k of **k{k}-out** stdout. First-step redo thresholds,
if printed, are all replayed at the same density with the preceding returned
orbitals. QE prints only three significant digits: e.g. 1.37e-3 Ry has a
rounding half-width of 5e-6 Ry. This command uses those printed values and
reports their precision; it does not claim the unprinted internal threshold
was recovered. `--ethr X` overrides all steps with an exact value **in Ry**;
conversion to Ha is exactly X/2. Missing, nonfinite or nonpositive thresholds
fail closed. No threshold is guessed from the step number.

JSONL: header/provenance, one complete row per k, final summary. File/solver
errors stop immediately; numerical disagreements retain all requested rows
and exit **4**, with the earliest failing k/term on stderr and in the summary.
Exit 0 means the density comparison passed its stated tolerance, not that an
SCF was certified. Check ordering: input charge, import roundtrip, output
charge, analytic Hartree consistency, occupations, band count, inner flags,
then output-density maximum. Charge allowance is 1e-10 electrons, analytic
Hartree check 1e-12 Ha. Energy/band differences are diagnostic; no energy
threshold is inferred from a density floor.

Each row includes:

- Reciprocal maximum/RMS, G=0, charge, ten lowest nonzero shells and real-space
  maximum/RMS. Reciprocal RMS counts the stored half-sphere entries once.
- `density_distance_ry` is **the mixer's squared norm**,
  `8*pi*Omega*sum_{half,G!=0}|rho_B-rho_QE|^2/G^2`, not its square root;
  the Ha value is half this. The two separate input/output dr2 values also appear.
- All printed final energy components, bands in eV, occupations (QE's printed
  fractions multiplied by 2), solver flags, residuals and solve work counts.
  H counts exclude atomic-start construction. Forced-stop `k-out` total energy
  uses the band-form functional; the direct orbital functional is separate.
- Input Hartree potential against the analytic 4*pi*rho(G)/G^2 formula. The
  QE XC/local potential grids are unavailable and reported null. Hartree/XC
  energies recomputed on **QE's output density** provide a separate density-only
  check against printed QE terms. The Hartree/XC/local *output response* split
  applies BAANDI kernels to the two output densities: it is explicitly **not**
  a comparison to independently exported QE potentials. The local response is
  zero because the ionic potential is density independent.

The diagnostic materializes full densities and final orbitals on every rank.
The reference solve itself uses distributed plane-wave rows/H for MPI; legacy
uses the existing dispatcher. This interface makes no scaling or memory-use
acceptance claim. Production defaults and hot loops are untouched.

## Mac fixture evidence (2026-11-23)

QE 7.6 trajectories were produced on dev08, not on this Mac. See the unchanged
[fixture provenance](../../../../tests/reference/qe/i25-si8/MANIFEST.txt),
[inputs/commands](../../../../tests/reference/qe/i25-si8/README.md), and
[stored complete map rows](../../../../crates/dft-cli/tests/fixtures/refscf-map-si8.json).
MANIFEST SHA-256: `8594d9e038dcab9c9a2f1b47340f93be2fdb4327d9c9a23e12064bdf9c7ef814`.
External Si.pz-vbc.UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

Si8, 20 bands (16 occupied), 30/120 Ry, 36^3 grid, 11,936 stored density G rows,
Omega=1080.4286448273454 bohr^3, seed 801, one worker. The comparison tolerance
is **2.72e-16 e/bohr^3**, four times the upper end of the **user/dev08-relayed**
QE np8-vs-np1 output floor (1.7e-17..6.8e-17). The np8 fixture is not on the Mac;
the floor was not independently remeasured here. This is a roundoff equality
diagnostic, not a justified cross-code loose-solver error allowance; no
physical tolerance was enlarged to turn the following discrepancies into passes.

| solver | k | printed ethr (Ry) | max rho difference (e/bohr^3) | RMS (e/bohr^3) | mixer distance (Ry) | inner iterations |
|---|---:|---:|---:|---:|---:|---:|
| legacy | 1 | 1.00e-2 | 1.412541591e-4 | 3.767465655e-6 | 1.387939016e-3 | 3 |
| legacy | 2 | 1.37e-3 | 4.016489166e-5 | 1.099621700e-6 | 1.068642280e-4 | 2 |
| legacy | 3 | 6.11e-5 | 1.180480939e-5 | 2.770282465e-7 | 2.627257499e-5 | 4 |
| reference | 1 | 1.00e-2 | 1.411562359e-5 | 3.582885826e-7 | 1.966403362e-5 | 3 |
| reference | 2 | 1.37e-3 | 1.026520345e-5 | 2.762159766e-7 | 1.282576976e-5 | 1 |
| reference | 3 | 6.11e-5 | 4.814860351e-6 | 1.552336363e-7 | 5.596930412e-6 | 4 |

Both commands exit 4, first failing term **k=1 rho_out all G**; all steps and
shells remain available. Input roundtrip maximum <=4.423e-18 e/bohr^3;
charge deviation <=3.553e-14 electrons; input Hartree analytic difference
<=8.442e-17 Ha. Occupations and all inner convergence flags pass.
BAANDI Hartree and XC energies evaluated on QE's output agree with each printed
QE value within its 5e-9 Ry decimal half-width. Ewald differs by 1.203e-8 Ry;
that independently reported term is not hidden in a density tolerance.

Internal tight comparison (`--ethr 1e-12`, all three common inputs): largest
occupied-band difference between reference and legacy is 1.016e-11 eV.
Direct orbital functional differences are 6.127e-8, 6.550e-10, 9.309e-9 Ha.
At these non-self-consistent inputs the direct functional is not stationary
in the input-H eigenvectors; eigenvalue agreement does not imply the same
second-order bound for that functional. The test limits (set before measuring)
are 1e-8 eV on occupied eigenvalues and 1e-6 Ha on the direct functional.
This is an internal regression check, not the strict SCF accuracy gate.

Reference k1 on simulated ranks 2/3 differs from its serial output by at most
2.914e-16 / 2.498e-16 e/bohr^3 in **real space**, with identical iteration/H
counts. This is a separate numerical quantity from the relayed reciprocal QE
floor. Native MPI/MKL and Si128 remain unmeasured here.

Default-path A/B against a pre-edit same-HEAD binary: Gamma, k222s, Si8 all have
**zero energy difference, identical density SHA-256, unchanged 4/8/4 SCF
iterations** and fixed-point verification enabled. Single-run timings under
compilation load are deliberately excluded from the evidence record.

## Local reproduction and gates

Scratch stays inside this worktree. Never copy the external UPF into git.

```sh
mkdir -p target/scratch
export DFT_UPF=/absolute/path/to/Si.pz-vbc.UPF
python3 - <<'PY'
import os
from pathlib import Path
s = Path('docs/performance/lanes/RSCF/si8-qe-replay.toml').read_text()
s = s.replace('/path/to/Si.pz-vbc.UPF', str(Path(os.environ['DFT_UPF']).resolve()))
Path('target/scratch/si8-map.toml').write_text(s + '\nnbands = 20\n')
PY
cargo build --release -p dft-cli
target/release/dft-cli refscf-map target/scratch/si8-map.toml \
  --qe tests/reference/qe/i25-si8/si8-np1 --steps 3 \
  --ethr-from-stdout --solver legacy > target/scratch/map-legacy.jsonl
# Expected diagnostic exit 4; run the reference arm separately.
target/release/dft-cli refscf-map target/scratch/si8-map.toml \
  --qe tests/reference/qe/i25-si8/si8-np1 --steps 3 \
  --ethr-from-stdout --solver reference > target/scratch/map-reference.jsonl
TMPDIR="$PWD/target/scratch" cargo test --release -p dft-cli \
  --test refscf_map -- --include-ignored --nocapture
cargo test -p dft-cli --bin dft-cli refscf_map::tests -- --nocapture
CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh full
```

Release diagnostic tests: exit 0, 4 passed (material/simulated-rank tests are
ignored by default to preserve the debug gate budget); parser tests: 2 passed.
The measured reciprocal max/RMS discrepancies are regression-checked against
the stored evidence within 1e-12 e/bohr^3, separately from the unchanged CLI
equality tolerance of 2.72e-16 e/bohr^3. Full gate: exit 1;
build/fmt/clippy pass, 827 tests pass, 188 ignored, **3 existing pseudopotential
cache tests fail** because TMPDIR is inside the repository. The sandbox also
forbids the `ps` load probe; the gate continues. No scratch was moved outside
the worktree to bypass these restrictions. README freshness test: exit 1,
8 pass / 1 existing assertion requiring the obsolete phrase
`standard profile is unsupported`. README, its generator, dashboard/registry
artifacts and delegation script are outside this ticket's owned files and
were not changed.

## dev08 commands (master session only)

Run after review/merge, from that checkout. Save the commit, build command,
binary SHA-256, inputs, UPF hashes and complete stdout/stderr/exit status.
Use the existing MPI/MKL installation; no installs or QE source edits.
Input `[scf].threads` must be 1 for these one-worker-per-rank probes, with
identical atom coordinates, cutoff, FFT grid, UPF and bands to the QE run.

```sh
mkdir -p target/scratch/rscf-map
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli \
  --features mpi,mkl --target-dir target/rscf-map
MAP_BIN="$PWD/target/rscf-map/release/dft-cli"
sha256sum "$MAP_BIN"
git rev-parse HEAD

# Both solver arms, each k uses the ethr printed in that k-out stdout.
# Si128 np1 and np8 intentionally consume the SAME QE np8 input trajectory:
# this isolates BAANDI rank dependence; there is no cross-model speed ratio.
for solver in legacy reference; do
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -np 1 --bind-to core \
    "$MAP_BIN" refscf-map /data/sig-dft/i25/si128/input.toml \
    --qe /data/sig-dft/i25/si128 --steps 6 --ethr-from-stdout --solver "$solver" \
    > "target/scratch/rscf-map/si128-np1-$solver.jsonl" \
    2> "target/scratch/rscf-map/si128-np1-$solver.stderr"
  printf '%s\n' "$?" > "target/scratch/rscf-map/si128-np1-$solver.exit"
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -np 8 --bind-to core \
    "$MAP_BIN" refscf-map /data/sig-dft/i25/si128/input.toml \
    --qe /data/sig-dft/i25/si128 --steps 6 --ethr-from-stdout --solver "$solver" \
    > "target/scratch/rscf-map/si128-np8-$solver.jsonl" \
    2> "target/scratch/rscf-map/si128-np8-$solver.stderr"
  printf '%s\n' "$?" > "target/scratch/rscf-map/si128-np8-$solver.exit"
  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 mpirun -np 1 --bind-to core \
    "$MAP_BIN" refscf-map /data/sig-dft/i25/si8-np1/input.toml \
    --qe /data/sig-dft/i25/si8-np1 --steps 6 --ethr-from-stdout --solver "$solver" \
    > "target/scratch/rscf-map/si8-np1-$solver.jsonl" \
    2> "target/scratch/rscf-map/si8-np1-$solver.stderr"
  printf '%s\n' "$?" > "target/scratch/rscf-map/si8-np1-$solver.exit"
done
```

Do not run that loop under `set -e`: exit 4 is an expected numerical finding,
not missing output. If an actual QE **Si128 np1 trajectory** is already
available, repeat the np1 probe with its own matching input and directory,
label it separately, and retain both datasets; do not relabel the np8 source.

The exact Si128 and k4..k6 stdout thresholds were not delivered to this Mac.
No numbers are invented for them: `--ethr-from-stdout` selects each supplied
value, and each JSON row contains `qe_printed_ethr_ry`, rounding half-width,
`solve_thresholds_ry` and the Ha conversion. The delivered Si8 k1..k3 values
are **1.00e-2, 1.37e-3, 6.11e-5 Ry**. Return this per-step table for every
remote series, including any first-step redo. If stdout is absent, provide it
as `k{k}-out/stdout.txt` (or root `k{k}-out.out`/`.stdout`); do not silently
borrow a mixing-run threshold.

Return all JSONL rows and the first failing term, all shells, charge,
occupations, bands, energy terms and potential checks. Also return the original
QE commands/version/stdout/input/UPF SHA-256, and SHA-256 of every density file,
so the floor and threshold provenance can be audited. No UPF payload in git.

If the same-input mismatch remains: repeat both arms with `--ethr 1e-12`, then
with `--cold-start`. Compare only correctness here. For attribution beyond this
experiment, the remote session must supply a tight **QE map at each same
input**, or compatible start orbitals plus exact unrounded thresholds; tight
BAANDI versus the existing loose QE outputs alone cannot distinguish potential
from finite inner-solve error. Independently exported QE Hartree/XC/local
potential grids would complete the missing direct potential split.
