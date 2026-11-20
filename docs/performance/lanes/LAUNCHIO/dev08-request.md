# LANE-LAUNCHIO Phase 2: executable digest cache A/B request

Status: **DEV08_AB_REQUEST**. Candidate remains **default off**. No Mac timing,
remote speedup or adoption is claimed. Diagnosis and alternative candidates:
[diagnosis.md](diagnosis.md). Measurement jobs always retain binary provenance.

## Candidate contract

`BAANDI_BINARY_HASH_CACHE=on` (or `[performance] binary_hash_cache = true`)
enables a Unix-only cache next to the executable, at
`<current_exe>.baandi-sha256.json`. Off/unset uses the unchanged fresh streaming
SHA-256 implementation (once per process), creates/reads no sidecar, and retains
all existing output behavior. Non-Unix builds retain fresh hashing.

The key contains device, inode, byte length, nanosecond mtime **and ctime**, owner
UID and mode. It is read from an open executable file descriptor, then compared
with that descriptor and the current path before accepting a hit or publishing
a fresh digest. Content replacement, same-length rewriting with restored mtime,
chmod and timestamp changes invalidate the key. At most 4096 sidecar bytes are
read. A hit must have the recognized schema, exact metadata key, lowercase
64-character hexadecimal SHA-256 and a checksum over its key/digest. The sidecar
must be a regular non-symlink file owned by the executable owner and not group/
other writable. The checksum detects accidental record edits, not malicious ones.

A missing, malformed, stale, oversized or unsuitable sidecar causes a **fresh
hash**, preserving the provenance value. A fresh record is best-effort published
via a unique create-new, mode-0600 temporary and atomic rename. Concurrent MPI
misses may each compute the correct hash, but no partial record is visible.
No synchronization/collective, root-only output or extra runnable thread is
introduced. Failed publication leaves the freshly computed SHA-256 in results.
Only the created temporary is removed on error; a symlink/directory/foreign or
writable sidecar is not replaced. No fsync: loss after a crash is a cache miss.

A changed executable detected during hashing yields unavailable identity rather
than publishing a stale record; the measurement gate rejects missing provenance.
The ordinary unchanged-binary cache path never omits `run.binary.sha256`.

### Trust and scope limits

Metadata is a filesystem coherence key, **not a cryptographic proof of content**.
Use only a trusted, local build directory and a filesystem with coherent inode,
nanosecond timestamp and ctime semantics. Deliberate sidecar forgery by the same
UID/root, adversarial directory replacement, metadata-preserving storage damage,
or weak/incoherent network metadata are outside this opt-in cache contract.
A record checksum is not authentication. Use OFF in those environments. The
measurement master independently hashes the executable before and after each
job and requires **every rank's recorded digest to equal it**. This is required
for acceptance, not replaced by trusting the cache's own checksum.

The prior process OnceLock remains only on the off path. The on path rechecks
metadata for every identity request, so scoped library runs cannot reuse a cache
hit from another configuration without validation. Binary hash timing still
includes lookup/fresh hash/publication; it is not moved outside the process timer.

Output addition, on only: `run.binary.hash_cache` contains:

- `status`: `hit`, `miss` or `unavailable`.
- `hashed_bytes`: **0 on hit**, actual streamed executable bytes on miss; null
  when provenance cannot be computed.
- `published`: whether a fresh sidecar was successfully published.

The provenance path/digest and all density/scientific files are unchanged. Only
switch manifests/agreement hashes, nondeterministic timestamps/timers/run IDs,
and explicit cache diagnostics may differ in A/B output. A/B uses one binary.

## Pre-registered gates

1. **Mechanism:** every retained warm B rank must report `status=hit` and
   `hashed_bytes=0`; A uses fresh hashing (sidecar ignored). Cold B is separately
   retained with full hashing on cache misses. Diagnostic rank-max
   `launch.seconds.binary_sha256` median is expected to decrease by at least **75%** versus
   A. This is a diagnostic forecast, **not an adoption percentage gate**: a
   smaller timer drop still qualifies if whole-run gain is resolved. Investigate
   metadata/cache overhead when the forecast misses.
   A warm cache that misses or returns no digest is not accepted as this mechanism.
2. **Correctness/provenance:** all recorded digests equal the independently
   computed binary SHA-256, before and after the job. Input/UPF hashes unchanged;
   converged, verified results; exact density.bin/density.json and full normalized
   results JSON per rank. Preserve full scientific histories, work counts and
   verification. No tolerance is relaxed. Reject non-distributed MPI fallbacks.
3. **A/A first:** six retained off/off pairs after warmups. Return paired SD, SE,
   mean drift and every sample. If absolute mean/SE >=3 for A/A, stop before A/B
   and repeat a complete job in a new quiet window; do not subtract drift away.
4. **A/B verdict:** 12 crossed pairs (AB, BA alternating), same binary/model/
   ranks/cores/binding/input and output filesystem. Let savings `d_i=A_i-B_i`.
   Report mean, sample SD, `SE=SD/sqrt(n)`, medians and raw pairs. Require positive
   mean, lower B median and **abs(mean)/SE >=3** on external launch-to-exit wall.
   Zero/undefined SE is unresolved rather than infinite evidence. Compare the
   savings to the observed A/A floor; report any unexpected floor/variance change.
   No minimum percentage gain. An unresolved result stays off; no cherry-picked
   extension or favorable subset. A resolved regression or correctness/provenance
   failure withdraws this implementation pending a new premise.
5. All four requested cells must pass correctness. Adoption requires resolved
   improvement without regressions beyond run-to-run noise across the measured
   set. Whole-process wall includes verification, serialization and teardown;
   improvement in a nested hash span alone cannot establish adoption.

Expected warm saving **0.05–0.10 s**, an estimate from the relayed .094–.105 s
fresh hash cost. Cold execution pays the hash plus publication and may regress.
This request qualifies **repeated runs of the same unchanged executable on a
warm sidecar**. It cannot establish a benefit for one-use or continually rebuilt
executables. Retain the cold measurement instead of hiding initialization cost.

## Exact master commands

Build the reviewed source with the **same** MKL/MPI toolchain and release flags
as e947a35, outside the exclusive measurement window. No remote action is
performed by this implementer.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cargo build --offline --release -p dft-cli --features mpi,mkl
sha256sum target/release/dft-cli
rustc --version
git rev-parse HEAD
git diff
lscpu -e
```

Use the **same input files and UPFs as the relayed profile**, one thread per rank,
with unchanged grid/cutoff/bands/start/geometry/threshold/verification. Stage
absolute-UPF-path inputs at the paths below. These are staging contracts, not
claims of existing remote files. Preserve inherited production switches in
manifest; do not introduce another optimization in either arm. Confirm CPU IDs
0–7 are distinct physical cores on socket 0 and inspect binding reports. No
cross-model comparison or new QE measurement is requested.

One command per exclusive job; output paths must be fresh:

```sh
python3 -I docs/performance/lanes/LAUNCHIO/phase2_job.py --case si16 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si16.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/phase2-si16-np8
python3 -I docs/performance/lanes/LAUNCHIO/phase2_job.py --case si32 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si32.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/phase2-si32-np8
python3 -I docs/performance/lanes/LAUNCHIO/phase2_job.py --case si64 --ranks 8 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si64.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/phase2-si64-np8
python3 -I docs/performance/lanes/LAUNCHIO/phase2_job.py --case si64 --ranks 1 --input /data/sig-dft/lanes/LAUNCHIO/inputs/si64.toml --binary target/release/dft-cli --output /data/sig-dft/lanes/LAUNCHIO/phase2-si64-onecore
```

The script is Linux/dev08-only and requires Python 3.11+. It copies the reviewed
executable into the **same private path for A and B**, independently hashes it,
and records the original binary/input/UPF hashes and complete argv/environment.
No sidecar is shared with unrelated builds. It runs:

- A warmup, one cold B (misses and any concurrent hits allowed), one warm B.
- Six A/A pairs, then twelve crossed A/B pairs with
  `BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off` in both arms. A is
  `BAANDI_BINARY_HASH_CACHE=off`, B is `on`. These are the acceptance samples.
- Three separate crossed **diagnostic** pairs with phase/stage timing on in both
  arms; these profiles are never pooled with acceptance wall samples.

MPI uses `mpirun -np 8 --cpu-set 0-7 --map-by core --bind-to core --report-bindings`;
one core uses `taskset -c 0` directly with the same MPI-capable binary. Load must
be <0.30 before each launch (bounded ten-minute wait). No overlapping jobs or
builds. Failed gates preserve outputs and stop. The master's normal exclusive
window checks remain necessary; affinity cannot be inferred from low load alone.

## Return package and unresolved attribution

Return all original result/density files for **every rank**, stdout/stderr with
binding, `manifest.json` (wall samples and summary), the sidecar, independent
checksums, build log, source commit/diff and hardware mapping. Never add nested
span maxima. Diagnostic comparison: binary_sha256/run_block, density_bin_write,
materialize, results_write, mpi_finalize and external wall separately.

Also return `launchio.sizes.density_bytes`, grid dimensions, executable byte size,
output filesystem/mount type and any existing resource-rank/page-fault data from
Phase 1. They were omitted from the textual relay and are needed to put numeric
bandwidth units on the exact workload formulas in diagnosis.md. There is no
unmeasured gather inside density_bin_write to optimize. If write attribution
remains necessary, ask for a separate filesystem/syscall diagnostic after this
single-lever A/B; do not add tracing to the acceptance samples.

Stop after this request. Default stays off until dev08 approval.
