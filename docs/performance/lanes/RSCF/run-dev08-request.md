# R-SCF-2: full reference SCF trajectories on dev08

Execute only through the master/dft5090 session after review. No server action
was performed by the implementer. This request collects correctness and work
records; **no timing comparison, speed claim or default promotion** is authorized
by the Mac evidence. The Si16 simulated-rank tail exception in
[patch2-validation.md](patch2-validation.md) remains open. Do not hide or relax it.

## Jobs and identity

One case per job; every job name starts `final-`, wall limit **7200 s**. Do not
bundle multiple cells in one job. Use the existing dev08 launch environment,
MPI/MKL installations and socket/core binding. No installs or server changes.

| Job | Cell | BAANDI ranks | Matching QE ranks | Threads per rank |
|---|---|---:|---:|---:|
| final-rscf2-si128-np8 | Si128 | 8 | 8 | 1 |
| final-rscf2-si64-np8 | Si64 | 8 | 8 | 1 |
| final-rscf2-si216-np8 | Si216 | 8 | 8 | 1 |
| final-rscf2-si32-np1 | Si32 | 1 | 1 | 1 |
| final-rscf2-si64-np1 | Si64 | 1 | 1 | 1 |

Locate the matching **full** QE 7.6 runs already on dev08 through their manifests;
record actual paths and job IDs. The truncated k1..k6 mixing snapshots under
`/data/sig-dft/i25/` are not full-run substitutes. Do not invent a full-run path.
If a full reference is absent or its input differs, return that fact and run a
matching QE case as its own `final-rscf2-qe-<cell>-np<ranks>` job (7200 s).

Use exactly the QE lattice/coordinates, cutoff, FFT grid, NC UPF bytes, band
count, Gamma point, integer occupations, start, beta/history, threshold and
binding. Preserve the original geometry displacement and its seed; do not
substitute a nominal cell. Reference implementation currently supports LDA,
unpolarized NC without NLCC, no symmetry, final-cutoff atomic density and
atomic orbitals. QE must have the matching atomic start, no symmetry, plain
Broyden, beta 0.7 and history 8. Record QE's atomic-orbital count and any random
completion. A start mismatch is a diagnostic result, not matching identity.

Reference `[scf]` keys (merge into the existing table):

```toml
driver = "reference"
mixing = "broyden"
mixing_beta = 0.7
mixing_history = 8
standard_threshold_ha = 5e-9 # conv_thr = 1e-8 Ry
verify = "fixed_point"
threads = 1
max_iterations = 200
```

Set `nbands` explicitly to the QE value. Set `basis.fft_grid` explicitly. Do not
set a cutoff ladder or legacy tuning overrides. `BAANDI_REFERENCE_SCF=on` is an
alternative selector, not an extra solver configuration. All other legacy
SCF switch values must be recorded but have no control over this driver.

## Build and launch template

Build the reviewed patch with the existing toolchain (`cargo build --release
-p dft-cli --features mpi,mkl`); record commit plus dirty-diff SHA-256, binary
SHA-256, toolchain, MPI and MKL versions. Record QE executable/version/hash,
full input/hash, UPF hash, launch command and exit code. Si UPF identity used on
the Mac was `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Never copy UPFs into git.

For each separately scheduled case, bind `ranks` cores and run the equivalent
of the following with the matched absolute input/output paths and core list:

```sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR_990="$PWD/target/scratch"
mkdir -p "$TMPDIR_990/rscf2"
taskset -c "$cores" mpirun -np "$ranks" --map-by core --bind-to core \
  target/release/dft-cli run "$input" --output-dir "$output"
```

Use the existing rank wrapper to measure elapsed process lifetime and **each
rank's** maximum RSS, including verification, output serialization/write and
MPI finalization. Do not mistake the launcher's RSS for the ranks' RSS. Keep
rank-local results/stdout/stderr and all exit codes, also on failure. JSON
records Linux VmHWM through output preparation; external wrappers close the
remaining final-write/lifetime boundary. Internal array-payload bounds are
not RSS. Do not compare timings from different parallel models anywhere.

## Returned tables and acceptance evidence

Return raw `results.json`, `scf.log`, density hashes, rank files, manifests and
stdout, plus machine-readable and Markdown side-by-side tables. For every
SCF iteration pair, include:

- QE `estimated scf accuracy` and BAANDI dr2 in Ry and Ha; preserve QE printed
  precision, and mark unavailable full-precision data explicitly.
- Every requested/used ethr in Ry, first-step redo, each Davidson call and
  `avg # of iterations`; active-band widths/sum, H columns, reduced orders,
  restarts, unconverged flags and Ritz-change arrays on a rank mismatch.
- Mixing history occupancy, condition, beta and first failed term.
- Total variational energy in Ry/Ha, output functional energy, band sum,
  Hartree, XC, double counting and descf; do not compare unlike energy forms.
- Phase wall seconds: start, solve, density, potential, audit, mix, metric,
  transform, energy correction, verification and output. Keep inclusive map
  clocks separate. QE has no independent BAANDI certificate; its matching
  table entry must say unavailable, never zero fabricated work.
- First estimator crossing, all fixed-point attempts, final pass/fail,
  Coulomb error, residual bound, requested tolerance, work and seconds.
- Per-rank peak payload bytes and actual peak RSS, and final charge/density
  agreement. Verify all rank-local results, not just rank zero.

Check final energies against the existing strict legacy reference at the same
input (standard class 5e-9 Ha). Any missing strict oracle is reported as pending.
On a mismatch report the **first** divergent phase/iteration, with preceding
rows; retain failed evidence. The Mac Si8 comparison already differs in the
initial density and hence dr2 at step 1, while ethr/iteration counts are 0.01 Ry
and 3 in both codes. Request initial rho(G), potential diagonal, initial Ritz
values, first-step Ritz changes and output rho(G) to isolate startup/radial
transform effects from the new solver. QE source is never sent to Codex.

Patch-1 replay is passed by the user's measured dev08 floor: same-rank
`1e-17 + 1e-12 * amplitude` e/Bohr^3 for max/RMS; cross-rank diagnostic
`1e-16 + 1e-12 * amplitude`. G0 and charge checks remain exact/1e-12 e as
documented. These are density-replay limits, **not** permission to alter the
SCF threshold, Davidson `< ethr` rule, energy class or rank-trajectory gate.

A complete report states convergence/certificate and all gate outcomes,
including open failures. Keep the legacy default until acceptance is decided.
