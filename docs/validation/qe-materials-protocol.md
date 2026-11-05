# Reproducible material comparisons (P16A)

This protocol implements input generation and numerical reporting for fixed,
unpolarized NC calculations, including multi-species PBE and NLCC. It does not
qualify a dataset or choose acceptance caps. MF05/MF06 in
[stage-oracles.md](stage-oracles.md) own that decision. Spin and smearing requests
fail explicitly until a separate ticket extends the protocol.

Follow [reference-environment.md](../reference-environment.md): only the remote
`dft5090` session runs QE on dev08. Use the original QE 7.6 MPI/DFTI build,
`~/fq/builds/orig-mpi-dfti/bin/pw.x`, after `source ~/fq/env.sh`. Check node load
before running. No QE source or UPF payload belongs in this repository.

## Evidence to send before running

Send the generator and its `qe_common.py` sibling, the exact dataset manifest,
the per-run BAANDI `results.json`, generated QE input and sidecar metadata,
expected input/UPF digests, and exact run command. Send `density.json` and
`density.bin` too if density comparison will be performed on the remote host.
Archive the BAANDI input TOML and executable commit/version as well.

`tests/reference/materials/results.json` is a six-case smoke **summary**, not
CLI output: it lacks structure, FFT dimensions and spectra. Both tools reject
it with a remedy. Obtain each full result from a retained P13D run directory
(e.g. the `--work-dir` supplied to `scripts/run_materials_smoke.py`). If that
record is unavailable, rerun BAANDI with the corresponding material TOML and
preserve its complete output directory; do not reconstruct a result by guessing
the grid or filling spectra from the summary. The existing runner normally uses
a temporary directory unless `--work-dir` is provided.

The runner substitutes absolute UPF paths into a working TOML, so the full
result's `input_sha256` identifies that **executed** TOML; the summary identifies
the repository template. They need not match. Preserve both without rewriting.

## Pseudopotentials on dev08

The user authorized the remote session to download the PseudoDojo NC-SR v0.4
PBE standard UPF2 archive to `/scr/mino/rustdft/pseudo`. Its source URL is in
[the manifest](../data/nc-datasets.json). This task performs no download. The
manifest currently records per-file hashes, not a trusted archive digest: record
the downloaded archive digest as observed and compare it to an independently
provided expected digest if available; never call an invented digest verified.
Place the extracted files directly at:

| File | Expected SHA-256 |
| --- | --- |
| `/scr/mino/rustdft/pseudo/Si.upf` | `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282` |
| `/scr/mino/rustdft/pseudo/Mg.upf` | `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14` |
| `/scr/mino/rustdft/pseudo/O.upf` | `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11` |

Use `sha256sum` remotely to verify these bytes before a run. Return the output.
Si and O carry NLCC; the manifest marks Mg as **not** NLCC. QE gets that behavior
from each UPF, with no invented input switch. The generator never opens a UPF:
it checks the recorded result hash and basename against the manifest and copies
only metadata. The remote verification is essential to verify the files actually
used. Keep PseudoDojo's CC-BY-4.0 attribution/license and citations with the archive.

## Generate, then run without editing the input

From the repository root, with a retained BAANDI output directory at
`/scr/mino/rustdft/baandi/si-20`, run:

```sh
mkdir -p /scr/mino/rustdft/qe/si-20
python3 -I scripts/make_qe_input.py \
  --results /scr/mino/rustdft/baandi/si-20/results.json \
  --pseudo-dir /scr/mino/rustdft/pseudo \
  --out /scr/mino/rustdft/qe/si-20/qe.in \
  --prefix si20 --outdir /scr/mino/rustdft/qe/si-20/work \
  --disk-io nowf --conv-thr 1e-12
sha256sum /scr/mino/rustdft/qe/si-20/qe.in
```

Repeat using `si-30`, `si-40`, `mgo-40`, `mgo-50`, `mgo-60` and distinct prefixes
and output directories. The `--conv-thr` value is QE's SCF threshold in Ry; it is
recorded, not an acceptance cap. Do not type or edit a QE input by hand. Change
an option and regenerate it. Comments record the BAANDI commit, executed input
hash, manifest and UPF hashes, and a shell-quoted invocation (argv normalized to
`python3 -I`). The sidecar records generated input and source-result digests.

For each case, run this command with the absolute paths for that case and retain
its literal invocation in `command.txt`:

```sh
source ~/fq/env.sh
export PATH="$HOME/fq/builds/orig-mpi-dfti/bin:$PATH"
command -v pw.x
mpirun --version
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/time -v -o /scr/mino/rustdft/qe/si-20/time.txt \
  mpirun -np 1 --map-by core --bind-to core pw.x \
  -in /scr/mino/rustdft/qe/si-20/qe.in \
  > /scr/mino/rustdft/qe/si-20/qe.out 2>&1
```

Record the exit code immediately, UTC date, host, resolved executable, QE version
from the actual output, MPI version, environment/thread settings, build/library
versions and load observation. A version reported in an older environment note
is not a substitute for the current observation. Return `/usr/bin/time`'s wall
time and maximum RSS verbatim, including the timing file.

Generated inputs use `ecutwfc_Ry = 2 * ecutwfc_Ha`, `ecutrho = 4 * ecutwfc`,
BAANDI `nbnd` and all three FFT dimensions, `occupations='fixed'`, `nosym=.true.`,
and `diago_full_acc=.true.`. Materials reject cutoff and automatic FFT overrides.
PBE records request `input_dft='pbe'`; LDA retains the existing functional rule.
Atomic masses are 1.0 (irrelevant to this static SCF; no forces, stress or dynamics).
Cell matrix columns become QE rows; positions are fractional crystal coordinates.

The existing k-point rule is unchanged: automatic `(m+s)/n` grids use shifts
0 or 0.5 and disable inversion unless time reversal was requested; explicit
lists always set `noinv=.true.` because any reduction already happened upstream.
An automatic 1x1x1 grid uses complex Gamma, which the existing density reader
requires. This protocol is for accuracy comparison, not QE's best-performance
real Gamma benchmark.

**`disk_io='none'` is forbidden:** it suppresses the saved charge-density file
needed for density comparison. Manifest mode defaults to `nowf`, which retains
density without wavefunctions. `low`, `medium` and `high` are explicit options.
For byte compatibility, legacy calls with `--upf-name` retain their former `low`
default and comment layout. Omit `--upf-name` for materials: identities must come
from the manifest and results, never a hand-supplied filename.

## Files to return unchanged

Return the generated `qe.in` and `.meta.json`, `qe.out`, `command.txt`, exit code,
versions/environment/host/date records, observed input and UPF SHA-256 output,
archive digest and download provenance, `time.txt`, and
`work/<prefix>.save/charge-density.dat`. Retain the save directory XML alongside
the density if available. Preserve the BAANDI full results, executed TOML,
density metadata and binary in the same evidence bundle. Nothing is estimated,
rewritten, rounded further, or relabeled as a reference without this provenance.
If density is absent, report that fact and reason; do not invent it.

## Compare returned evidence

```sh
python3 -I scripts/compare_materials_qe.py \
  --ours /scr/mino/rustdft/baandi/si-20/results.json \
  --qe-out /scr/mino/rustdft/qe/si-20/qe.out \
  --qe-density /scr/mino/rustdft/qe/si-20/work/si20.save/charge-density.dat \
  --json /scr/mino/rustdft/qe/si-20/comparison.json
```

Omit `--qe-density` only when unavailable. When supplied, the BAANDI density
files must be siblings of its results. An unreadable/inconsistent density fails;
it is never silently skipped. Reports refuse to overwrite existing evidence.
The density path reuses `qe_rhog.py`'s independent discrete transform and binary
reader and `compare_rhog.py`'s metrics: coefficient differences, rho(0), Parseval,
and relative L2 on identical G sets. Its supported binary layout is little-endian
Fortran records, complex Gamma (`gamma_only=0`), one spin channel; other variants
fail explicitly. The transform is a correctness tool and can be slow for large grids.

The JSON reports signed BAANDI minus QE energy differences in Ha and Ry,
occupied/first-empty eigenvalue differences for each matched k point, direct and
indirect sampled band-gap differences, configured electron counts, both iteration
counts, and optional density metrics/integrals. Spectra are compared without an
energy-zero fit. QE stdout prints rounded eV eigenvalues, not full precision.
No empty band yields an unavailable gap, not a fabricated value.

The parser rejects incomplete or concatenated SCFs, unrecognized output layouts,
wrong grids/cutoffs/XC/band counts, and inconsistent k lists. Matching is modulo
reciprocal lattice integers, with time reversal only under the generator's rule.
The fixed coordinate/weight bounds describe stdout precision, not scientific
acceptance caps. The configured QE electron count is not an independent integral;
only the density path supplies one. Stdout alone cannot verify the executed input
or UPF bytes: manually audit the returned provenance against the sidecar before
classifying a comparison as a reference.

Exit 0 means **reported**, not scientifically accepted. Exit 3 means inconsistent
evidence, exit 1 means unsupported/malformed input or I/O/usage failure. There
are no adjustable scientific tolerances or PASS labels in this tool; apply the
separately frozen MF05/MF06 caps to its stored numbers.

## Local regression evidence

Run `python3 -I scripts/tests/test_materials_qe.py` and the existing QE, k-point,
PBE and density script tests. The new test stores 12 pre-change legacy input
byte goldens plus Si-20/MgO-40 inputs under
`tests/reference/materials/qe_inputs/`. Their source records are literal field
projections of retained P13D output, not new QE results. See
`scripts/tests/fixtures/README.md` for provenance. Synthetic stdout and binary
fixtures test known perturbations and failures; they are never QE references.
