# R3 real Gamma Hamiltonian: component evidence and Mac observation

2026-11-06, P26C-REAL-GAMMA-R3. Library-only addition, uncommitted on
`task/P26C-REAL-GAMMA-R3`, base/current-main commit
`e9febbd7ed54a61cfd9c42d02126c5d8ce57a9de` (baseline and main are the same commit).
Derived from `real-gamma-design.md` sections 4–7, itself derived from the
architect's Gamma study note. No QE source was read and no QE run was made.

## Mechanism and scope

A snapshot constructor validates exact basis identity and Gamma symmetry, then
builds P from the public complex projector representatives. The origin check
rejects imaginary magnitude above 1e-14 times the column maximum; exact zero
columns are accepted. P^T X and P W use the existing deterministic real GEMM;
row-major D is applied separately for each atom. Kinetic values are copied into
both ORC slots. R2 applies fixed adjacent local pairs with a zero odd partner.
The source Hamiltonian and NonlocalOperator kernels are unchanged; three
read-only Hamiltonian accessors are additive. There is no engine caller.

This is not a retry of A43's rejected complex row-blocked traversal: the new
operator uses real arithmetic, a band-block reuse axis and a numerical-equivalence
contract. A28's many-band FFT grouping is not used: two real bands share one
complex transform. No tolerance, SCF default, solver, or existing kernel changed.
Local pairs are serial (R7 owns pair scheduling); only GEMM uses the supplied
persistent worker pool. For m bands the local part performs ceil(m/2) transform
pairs, versus m for the complex path. With n rows and p projectors the two GEMMs
have nominal work 4*n*p*m real FLOPs, excluding the small per-atom D products.

All numerical scratch is constructed once. Applications use slices, matrix
views, fixed loops and the existing R1/R2 storage contracts. No new allocation,
logging, I/O, timer or unsafe appears in production application code. Output
storage retention is tested. This is a source/storage-contract audit, not a
measurement by a global allocator; hidden backend allocations were not measured.
Scratch for R1 Gram tiles scales as ceil(n/256)*max(p,capacity)^2 doubles.
Potential/projectors are snapshots; updating the source does not update them.

## Component evidence

`cargo test -p dft-physics --test real_gamma_hamiltonian -- --nocapture`
passed (exit 0). The ten cases have respectively (atoms, rows, bands):
(1,19,3), (2,21,4), (4,23,3), (8,71,4), (12,19,3), (16,25,4),
(20,23,3), (3,77,4), (6,27,3), (10,31,4). Geometries cover cubic,
skewed, elongated and slab-like cells, each with a different metric. Synthetic
species use repeated s/p/d/f radial projectors and a second species with s/p
channels; same-l off-diagonal D entries are nonzero. No UPF file is stored.

- Full H maximum normwise relative difference: 6.316026193846391e-16.
- Nonlocal maximum normwise relative difference: 5.613885122435422e-16.
- Maximum |x^T H y - y^T H x| / max(|x^T H y|, |y^T H x|):
  1.2849945776886336e-15.
- All output bits match at T=1,2,4,8 on every small case.

The fixed 1e-13 gates allow FFT pairing and long-sum GEMM rounding, as prescribed
in design section 11. They are not SCF energy tolerances. Comparisons unpack the
real image back to the full complex sphere; they do not discard the imaginary
part of the complex reference. Empty nonlocal, zero-width blocks, bad dimensions,
nonfinite input and the G=0 rejection helper are also tested.

`cargo test --release -p dft-physics --test real_gamma_hamiltonian large_parallel -- --ignored --nocapture`
is the reproducible large test (also executed together with the timing test,
exit 0). It uses 3,035 rows, 360 projector columns, 17 bands and 20 atoms:
18,574,200 multiply-adds per GEMM, crossing the actual four-panel / 1,048,576-MAC
parallel threshold. At T=1,2,4,8 output bits agree, and full-H relative difference
is 9.765725874762848e-16 at every budget. This case establishes the parallel
branch, whereas the ten small cases execute GEMM inline.

## Timing observation (no gain claim)

Host: macOS 27.0.1 (26A434), arm64, Rust 1.99.1, release pure-Rust backend.
CPU model query was denied by the sandbox. Shared-host load averages before and
after: 4.95, 6.00, 4.24. This task ran no concurrent build or other benchmark
during the retained timing run; other host activity was not controlled. An
initial timing invocation overlapped other task work and is not retained here.
No baseline/main/branch speed ratio, scaling claim or dev08 acceptance is made.

Command (exit 0):

```sh
cargo test --release -p dft-physics --test real_gamma_hamiltonian timing_blocks -- --ignored --nocapture --test-threads=1
```

One warmup per path; 30 alternating real/complex pairs; T=1; construction,
projector transforms and allocation excluded. Both paths apply the same H to
the same normalized real-structured bands. Timers surround whole block calls.
These synthetic shapes use the actual counts below, not the design's illustrative
257/2969/16385-row shapes. They are not converged material or SCF benchmarks.

| Atoms | Rows | Bands | FFT grid | Real median, s/block | Complex median, s/block |
|---|---|---|---|---|---|
| 4 | 81 | 9 | [12, 12, 12] | 6.21875e-05 | 0.0001078545 |
| 8 | 519 | 20 | [24, 24, 24] | 0.000876771 | 0.0024944585 |
| 20 | 3035 | 41 | [40, 40, 40] | 0.015950021 | 0.059384958 |

Raw paired samples in seconds (index identifies the alternating pair):

Atoms 4, rows 81, bands 9:

```text
real_s = [0.000143958, 6.0375e-05, 6.0583e-05, 6.0334e-05, 6.0375e-05, 5.975e-05, 6.0167e-05, 6.0209e-05, 6.0041e-05, 6.2292e-05, 6.225e-05, 6.175e-05, 6.25e-05, 6.25e-05, 6.2209e-05, 6.1459e-05, 6.1583e-05, 6.2333e-05, 6.2542e-05, 6.225e-05, 6.2542e-05, 6.1292e-05, 6.2166e-05, 6.125e-05, 6.2417e-05, 6.2625e-05, 6.25e-05, 6.2417e-05, 6.2417e-05, 6.2042e-05]
complex_s = [0.000103791, 0.000104375, 0.000104458, 0.000104625, 0.000105042, 0.000104791, 0.000105375, 0.000104459, 0.000108125, 0.000108959, 0.000108125, 0.0001085, 0.00010775, 0.000107959, 0.000108792, 0.00010825, 0.000107042, 0.000108292, 0.000107708, 0.00010725, 0.000109, 0.000106833, 0.000106625, 0.000106875, 0.000108916, 0.00010825, 0.000110542, 0.000108458, 0.000108583, 0.000110583]
```

Atoms 8, rows 519, bands 20:

```text
real_s = [0.000873, 0.000894875, 0.000901333, 0.000882, 0.000881958, 0.000881292, 0.000881625, 0.00087525, 0.000871667, 0.000876, 0.000885417, 0.000874833, 0.000865208, 0.000883666, 0.000891542, 0.000888125, 0.000871625, 0.000873958, 0.000885041, 0.00088275, 0.00083075, 0.000877625, 0.000877083, 0.000852708, 0.000880167, 0.000876459, 0.000828167, 0.000866333, 0.000866209, 0.000828625]
complex_s = [0.002489125, 0.0025395, 0.0024495, 0.002460083, 0.002497209, 0.002499458, 0.002490875, 0.002515791, 0.002519125, 0.002486542, 0.002535125, 0.002515834, 0.002498958, 0.002504375, 0.002468791, 0.002488084, 0.002499375, 0.002494292, 0.002527334, 0.002494, 0.002488166, 0.002495, 0.002535083, 0.0024115, 0.002507625, 0.002482791, 0.002422333, 0.002478333, 0.002494625, 0.002374292]
```

Atoms 20, rows 3035, bands 41:

```text
real_s = [0.0159505, 0.015984459, 0.015854166, 0.01616425, 0.015943875, 0.015948083, 0.015915833, 0.01595775, 0.0159, 0.015921875, 0.015860209, 0.016108042, 0.015867083, 0.016015167, 0.015900708, 0.016072084, 0.015927958, 0.016034041, 0.015957375, 0.015970958, 0.015865042, 0.016010333, 0.015917125, 0.015953209, 0.015879709, 0.016014625, 0.015886708, 0.015959209, 0.015949542, 0.015958083]
complex_s = [0.059290584, 0.059252, 0.059125958, 0.059323417, 0.059398959, 0.059291792, 0.059436, 0.059087958, 0.059160625, 0.059283584, 0.059900041, 0.060093291, 0.059375458, 0.059422334, 0.059552291, 0.059377292, 0.059382667, 0.059426709, 0.059384375, 0.059318667, 0.059385541, 0.05936125, 0.059489917, 0.059390291, 0.059416958, 0.059430666, 0.059503292, 0.059417625, 0.059476208, 0.059377375]
```

## Verification and handoff

`scripts/check.sh dft-physics` and `scripts/check.sh` both passed build, test,
fmt and clippy, exit 0. The first crate check exited 1 for an unused import;
the import was removed, and both final gates passed without relaxed warnings.
`python3 -I scripts/tests/test_readme_figures.py` passed nine tests (exit 0);
existing resource warnings came from the unchanged comparison reader.
`git diff --check` passed (exit 0).

The explicit ownership list excludes README generation, reference artifacts,
ledger/backlog, dashboard/session registry and optimization-log edits. The
orchestrator can link this evidence and add an R3 ledger row after review;
those files were not changed. No new stored external-reference dataset or
README numerical claim is introduced. The session registry has no R3 entry in
this worktree; no delegate report or recurring orchestration timer was present.

Dev08 request: run the two ignored release tests and the full checks on pure
and MKL builds, record commit/binary hashes and host load, and repeat component
identity at T=1,2,4,8 before any integration. For performance acceptance use
quiet-node alternating comparisons, at least 30 samples per named site, and
include the larger intended shapes; Mac observations cannot decide a default.
The external R0 profile, MKL runtime allocation audit, real Davidson and engine
SCF equivalence remain outside this library-only ticket.

## Existing-engine identity evidence

The pre-edit executable and final executable were built from the same HEAD,
with only this isolated library addition between them:

- Baseline binary SHA-256:
  `ed4a4a6b8ef80208ec4d90544fa985113382e6634f51c62d1dbc66ef4fddc4cf`.
- Final binary SHA-256:
  `6a1d8c8e2801d039f44c6f1293685cfe006af8102752545dfc8fb3c375ec2ae3`.
- External Si UPF SHA-256:
  `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

Both `cargo build --release -p dft-cli` commands exited 0. For each T in
1,2,4,8,12, the following baseline and comparison commands exited 0 (ten commands,
25 matched case/budget comparisons). The baseline was retained in `/tmp` before
building the changed source. All energy differences are exactly 0 Ha; density
SHA-256 and iteration counts agree with the corresponding baseline. Timing
columns emitted by this regression tool are not used as performance evidence.

```sh
python3 -I scripts/bench_scf.py --binary /tmp/p26c-baseline-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 1 --threads "$T" \
  --record "/tmp/p26c-threads-$T-baseline.json"
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 1 --threads "$T" \
  --baseline "/tmp/p26c-threads-$T-baseline.json" --require-identical \
  --record "/tmp/p26c-threads-$T-current.json"
```

Retained final energy and density evidence (same values across all five budgets):

| Case | Energy, Ha | SCF iterations | Density SHA-256 |
|---|---|---|---|
| gamma | -7.2955304682107602 | 10 | `6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26` |
| k222s | -7.9194234124690865 | 10 | `9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104` |
| k222 | -7.8312300571211351 | 10 | `21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe` |
| k444s | -7.9259761266480915 | 10 | `767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660` |
| si8 | -31.32180360103127 | 15 | `826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e` |
