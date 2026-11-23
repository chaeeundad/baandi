# R-SCF-1 / issue #25 patch 1: dev08 density replay

Execute through the master/dft5090 session. No server change or QE source access
was performed by this implementer. This is a numerical comparison request,
not a speed acceptance request. Do not compare MPI with threads in any column.

## Prerequisites and identity

Use the reviewed R-SCF-1 patch (record git SHA and working-tree diff hash).
Use the existing QE 7.6 trajectories under
`/data/sig-dft/i25/{si128,si8}/k{1..6}-{mix,out}/`. Each snapshot must contain
exactly one `charge-density.dat`, directly or under a `.save` subdirectory.
Place the corresponding BAANDI input at `/data/sig-dft/i25/si128/input.toml`
and `/data/sig-dft/i25/si8/input.toml` before the commands below. These must
match the **actual QE trajectory**, including displaced coordinates, lattice,
cutoff, FFT grid, NC UPF bytes, occupations, nbands, beta, history length and
symmetry; do not use a nominal benchmark geometry as a replacement. The replay
only needs geometry/cutoff/grid/mixing, but the SCF smoke needs the real UPF.
Keep both inputs, QE version/command/stdout, UPF SHA-256 and all density hashes
in the returned manifest. Copy/reference each k-step stdout as `k{k}-mix/stdout`
if its current name is not one of the names in the design note.

Each QE `k{k}-mix` must be the unconverged mixed x_next; `k{k}-out` must be the
corresponding y at the **same** x and inner-solve convention, not a different
trajectory's kth output. Document whether one-step extraction tightened ethr:
a changed y then measures extraction differences, not just mixer differences.
Supply `k1-in` when available and a beta=1 first-output directory named
`k1-beta1` when available. These optional files are validated, never invented.
The already relayed Si128 probe is ngm=191712, bytes=5368052, rho0=512/Omega;
independently check it in this delivery.

## Build and exact replay commands

Run from the reviewed checkout with the existing dev08 MPI/MKL environment.
The following only builds and writes task-local evidence; no bulk install.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export TMPDIR_990="$PWD/target/scratch"
mkdir -p "$TMPDIR_990/rscf"
cargo build --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
mpirun --version

taskset -c 0 mpirun -np 1 --bind-to core --map-by core \
  target/release/dft-cli refscf-check /data/sig-dft/i25/si128/input.toml \
  --qe /data/sig-dft/i25/si128 --steps 6 --atol 1e-12 --rtol 1e-9 \
  > "$TMPDIR_990/rscf/si128-np1.jsonl" 2> "$TMPDIR_990/rscf/si128-np1.stderr"
echo $? > "$TMPDIR_990/rscf/si128-np1.exit"
taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core \
  target/release/dft-cli refscf-check /data/sig-dft/i25/si128/input.toml \
  --qe /data/sig-dft/i25/si128 --steps 6 --atol 1e-12 --rtol 1e-9 \
  > "$TMPDIR_990/rscf/si128-np8.jsonl" 2> "$TMPDIR_990/rscf/si128-np8.stderr"
echo $? > "$TMPDIR_990/rscf/si128-np8.exit"
taskset -c 0 mpirun -np 1 --bind-to core --map-by core \
  target/release/dft-cli refscf-check /data/sig-dft/i25/si8/input.toml \
  --qe /data/sig-dft/i25/si8 --steps 6 --atol 1e-12 --rtol 1e-9 \
  > "$TMPDIR_990/rscf/si8-np1.jsonl" 2> "$TMPDIR_990/rscf/si8-np1.stderr"
echo $? > "$TMPDIR_990/rscf/si8-np1.exit"
taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core \
  target/release/dft-cli refscf-check /data/sig-dft/i25/si8/input.toml \
  --qe /data/sig-dft/i25/si8 --steps 6 --atol 1e-12 --rtol 1e-9 \
  > "$TMPDIR_990/rscf/si8-np8.jsonl" 2> "$TMPDIR_990/rscf/si8-np8.stderr"
echo $? > "$TMPDIR_990/rscf/si8-np8.exit"
```

These provisional limits are diagnostic. Do not widen them to obtain exit 0.
On a failure return the **first failing term** and all preceding rows before
changing anything. Assess whether a mismatch comes from units, normalization,
G0, membership, source trajectory/redo, conditioning or arithmetic. Preserve
all raw failures. Native MKL factorization/reduction differences are unmeasured
on the Mac; compare np1 and np8 errors on identical external x/F data.

## Independent QE floor and compact fixture

Return QE-against-QE repeated same-rank and 1-vs-8-rank density differences,
with the commands, printed precision, binary/version, inputs and UPF digest.
The note requires this floor before a measured tolerance can be named. Return
per-k dr2, next-density max/RMS, G0, total charge, lowest ten shells, gamma
norms and A condition estimates for the BAANDI replay. Include the k1 input
recovery and beta=1 cross-check results. A null printed estimate must remain
explicitly unavailable.

Deliver Si8 on `data/i25-si8-qe` under `tests/reference/qe/i25-si8/` with
`input.toml`, `k1..k6-{out,mix}` snapshots, stdout and a provenance manifest.
No UPF is committed. The test uses the actual CLI and emits an explicit skip
while absent. Input paths for the replay fixture need not resolve a UPF (the
comparison parses density data without running the Hamiltonian).

## New-driver smoke (not verification or timing acceptance)

For an explicitly unverified trajectory, put `driver="reference"` and
`verify="off"` in the matched input's existing `[scf]` table. Leave beta/history
at their QE values, start density atomic, and use the final cutoff. Then run:

```sh
BAANDI_REFERENCE_SCF=on taskset -c 0 mpirun -np 1 --bind-to core --map-by core \
  target/release/dft-cli run /data/sig-dft/i25/si8/input.toml \
  --output-dir "$TMPDIR_990/rscf/si8-driver-np1"
BAANDI_REFERENCE_SCF=on taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core \
  target/release/dft-cli run /data/sig-dft/i25/si8/input.toml \
  --output-dir "$TMPDIR_990/rscf/si8-driver-np8"
```

Also run the same input with verification omitted (the default): it must stop
with `converged=false`, exit 3, and the explicit patch-2 verification TODO.
Return complete results.json/scf.log, reference_scf records, density hashes,
per-rank agreement, real-Gamma evidence and output/error exit codes. The MPI
skeleton uses owned LDA kernels; unsupported NLCC/PBE/symmetry MPI cases must
produce input errors, not fallback. Initial construction remains shared engine
code; only the new density/mixer and iterative map inputs are distributed.
No performance conclusion is allowed before patch 2 supplies verification
inside the total measured wall and the generality/accuracy gates pass.

## Delivered-data follow-up

Job `final-i25-100-qe-ref` is complete. Use the actual directories
`/data/sig-dft/i25/{si128,si8,si8-np1}/k{1..6}-{mix,out}/` and optional
`k1-mix-beta{0.3,1.0}/`; the CLI recognizes these beta directory names.
The committed Si8 one-rank k1..3 regression passes at `--atol 1e-19 --rtol
1e-13`. For the full delivered trajectories, execute after the existing
`cargo build --release -p dft-cli --features mpi,mkl` command:

```sh
export TMPDIR_990="$PWD/target/scratch"
mkdir -p "$TMPDIR_990/rscf"
for system in si128 si8 si8-np1; do
  for np in 1 8; do
    mpirun -np "$np" --bind-to core --map-by core \
      target/release/dft-cli refscf-check "/data/sig-dft/i25/$system/input.toml" \
      --qe "/data/sig-dft/i25/$system" --steps 6 --atol 1e-19 --rtol 1e-13 \
      > "$TMPDIR_990/rscf/$system-np$np.jsonl" \
      2> "$TMPDIR_990/rscf/$system-np$np.stderr"
    echo "$system np=$np exit=$?"
  done
 done
```

Prepare the matched `input.toml` as described above; those TOMLs are BAANDI
inputs, not delivered QE artifacts. For Si8, the committed
`docs/performance/lanes/RSCF/si8-qe-replay.toml` supplies the exact cell,
positions and 36^3 grid; replace its UPF path. Return every exit code,
JSONL/stderr, input and UPF hashes, binary build/commit, MPI/MKL versions,
thread/binding settings, and any first failing term. Do not loosen the
measured tolerance on failure: first diagnose ordering/units/G0 and report
rank-count rounding or conditioning. These are correctness replays, not
cross-model speed comparisons.
