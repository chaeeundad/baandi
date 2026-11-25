# G2-PBE: legacy/reference PBE at MPI8

Status: implementation and Mac correctness validation; remote measurements
pending. Master session delegates these jobs to baandi-dev08. Codex performs
no remote login, install, QE source access or UPF check-in.

## Full-convergence screen: one run per driver per cell

Use a reviewed commit containing G2 for both arms. Build once with the pinned
release MPI/MKL configuration (`cargo build --release -p dft-cli --features
mpi,mkl`); retain the exact command, source/diff SHA-256, commit, compiler/MPI/
MKL versions and binary SHA-256. Do not time concurrent builds.

1. Si Gamma PBE NC **without NLCC**, from a locally available licensed PBE
   dataset with atomic charge and atomic wavefunctions. Validate its header
   and record full UPF SHA-256 before execution. Do not substitute Si.pz-vbc
   (LDA), remove a core table, or silently override the functional. If no
   qualifying real file exists, mark this real-Si screen blocked on the file;
   the synthetic diamond Mac test is separate internal evidence.
2. After G1 merges, a PseudoDojo PBE NLCC cell (prefer MgO for non-Si coverage;
   Si is an additional alternative). Use the existing manifest and archived
   input under `docs/data/nc-datasets.json`; record actual input and UPF hashes.
   First merge G1's combined-density XC handling with G2's combined-density
   gradient, valence-only Hartree and valence-only double counting. Confirm
   the merged start/map path and any density-floor restrictions locally before
   running this arm. G2 alone still rejects core charge.

For each eligible cell execute legacy then reference, exactly once each,
MPI8 with one thread per rank and identical socket-0 physical-core binding.
Use the archived MPI launcher/binding and retain its actual binding report.
Both inputs must specify identical cell/positions, cutoff, FFT dimensions,
UPFs, band count, atomic start, PBE, beta/history, convergence class and
fixed-point verification. Set explicit final cutoff (no ladder); use
`scf.driver="legacy"` versus `scf.driver="reference"`; explicitly set
`BAANDI_REFERENCE_SCF=off` and `BAANDI_REFERENCE_STRICT=off` in both arms so
that launch environment cannot override driver selection. The reference mixer/inner solver differ by design;
record both resolved configurations. Keep the optional RMM, verification
reuse and warm-verification variants off, and record LEAN_ENERGY. Use the
standard 5e-9 Ha class (reference threshold 1e-8 Ry), no tolerance relaxation.

```sh
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -p \
  mpirun -np 8 <archived-socket-0-binding> <candidate> run <legacy.toml> --output-dir <legacy-out>
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /usr/bin/time -p \
  mpirun -np 8 <archived-socket-0-binding> <candidate> run <reference.toml> --output-dir <reference-out>
<candidate> refscf-check --stop-to-verified <reference-out>/results.json
```

Each full-convergence job is fast class, timeout 600 s. No warmup SCF,
repetition, extra QE run or new seed campaign is requested. Preserve failures
and timeouts. Report unavailable matching strict-oracle evidence as missing;
do not label legacy/reference agreement alone as independent true accuracy.

Return external wall including startup, verification failures/recovery,
output and MPI exit; SCF/Davidson/H counts; verification attempts and verdicts;
energy terms, total E and signed E_reference-E_legacy in Ha; matching archived
strict error if available; density/eigen residuals, grid, charge, input/binary/
UPF hashes, load/binding, exit codes and raw stdout/stderr/results. Standard
agreement gate is 5e-9 Ha; a separate tight comparison, when matching archived
strict data are available, is 1e-10 Ha. Report rank-local field/FFT counters and
peak RSS if already available; do not add rank maxima into a fictitious wall.
One run each is a functionality/cost screen, not a median or speed adoption.
There is no QE speed ratio and no comparison across parallel models.

## Task-local attempt and integration record

This is capability generalization, not a promoted performance optimization.
The library formulas, legacy path, LDA arithmetic, LEAN_ENERGY decisions and
verification criteria remain unchanged. See the [PBE design section](../../../design/reference-scf.md#pbe)
for derivative normalization, Nyquist handling, exact transform counts and
energy conventions. No rejected optimization is retried.

Shared README generator/evidence, dashboard/session JSON, attempts ledger and
optimization log are outside this ticket's OWNED_FILES. The master should
reconcile this record after review; no new generated README result or remote
acceptance is claimed by this worktree.

## Mac validation

Commands, exit codes and measured numerical results are recorded below. Full SCF matrix tests are ignored in routine debug gates:

```sh
cargo test --release -p dft-engine --test refscf pbe_legacy -- --ignored --nocapture --test-threads=1
cargo test -p dft-engine --test refscf pbe_distributed -- --nocapture
```

### Completed SCF comparison and limitations

On the Mac pure-Rust backend, the final release matrix passed 24 reference
SCFs (two synthetic systems x ranks 1/3/8 x standard/tight x LEAN_ENERGY
off/on), plus two independently started legacy strict oracles. Both reference
profiles passed fresh fixed-point verification. LEAN_ENERGY preserved every
delivered density value and total-energy bit; full band/Harris–Foulkes audits
passed. No production parameter or tolerance changed.

| Synthetic case | Standard max absolute delta vs legacy (Ha) | Tight max absolute delta (Ha) | Maps standard/tight |
| --- | ---: | ---: | ---: |
| Two atoms in Si diamond geometry, Gaussian local PP, z=2, rc=1.8 Bohr | 3.927607128e-10 | 1.687538997e-14 | 8 / 13 |
| One atom, skew 7-Bohr cell, same local PP | 1.379742987e-10 | 1.598721155e-14 | 7 / 10 |

Thresholds are 1e-8 and 1e-12 Ry; comparison caps are respectively 5e-9 and
1e-10 Ha. Maximum same-profile rank-count energy difference is
1.332267630e-15 Ha. Legacy uses its strict deterministic independent orbital
start, reference uses its existing atomic start. Synthetic species are Xx,
not a real Si PBE dataset. These figures are not QE evidence or timing claims.

An additional two-atom skew fixture with strong synthetic nonlocal projectors
fails closed when Broyden produces negative PBE density (both beta=0.7 and an
explicit trial beta=0.3). The collective rejection is now a regression test;
no absolute-value, clipping, floor, damping default or tolerance workaround
was introduced. Nonlocal PBE *start fields* still agree with the serial oracle
on all grid/rank fixtures. Successful nonlocal full-SCF generality is **not**
established by the two successful local-PP examples.

A trial widening the Gaussian PP/atomic density to rc=2.5 Bohr did not solve
that nonlocal rejection. Its local diamond case converged to an energy
0.004378060656018645 Ha away from the independently started legacy result;
that trial is not accepted as a unique-solution oracle. Its discrepancy is
unresolved (different fixed points are a possibility, not established here).
The final tests retain rc=1.8, and the failed trial is recorded, not counted
as agreement. All four exploratory failing commands exited 101; raw logs
remain under `target/g2-scf*.log` in this worktree.

For the remote request use modified PW92 and no `compat.mode="reference"`
override requiring original PW92; original PW92 remains explicitly rejected
by this reference start. Default Broyden can fail on signed mixed PBE density;
archive that failure instead of treating it as permission to alter physics.

### Field oracles and local command exits

Three skew-cell grids (8x10x12, 9x11x7, 7x8x9), each on simulated ranks 1/3/8,
include Nyquist modes, uneven/empty ownership and synthetic nonlocal start.
A communicator wrapper prohibits gather/broadcast and nonscalar reductions
during GGA evaluation. Against analytic derivatives and the unchanged serial
GGA evaluator, measured maximum errors are:

- gradient: 7.307522643e-17 Bohr^-4 (cap 3e-15);
- potential: 1.665334537e-16 Ha (cap 2e-12);
- XC energy integral: 6.661338148e-16 Ha (cap 2e-12);
- central energy derivative versus distributed Vxc inner product:
  6.585958445e-10 Ha (cap 1e-8, step 1e-6 Bohr^-3; the tested directional
  derivative is explicitly nonzero).

Start potentials, including Hartree plus the full GGA divergence, agree with
the serial oracle to 2e-12 Ha. A negative input on only one rank fails on all
ranks before FFTs. No real-space full-grid payload is gathered by production
GGA; full arrays are constructed only by test oracles.

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace/tests compile |
| `cargo test -p dft-engine --test refscf -- --nocapture` | 0 | 62 passed, 9 ignored |
| `cargo test -p dft-physics` | 0 | 97 passed, 18 ignored across binaries; doc tests pass |
| `cargo test --release -p dft-engine --test refscf pbe_legacy -- --ignored --nocapture --test-threads=1` | 0 | Both fixtures, 24 reference SCFs, legacy oracles |
| `cargo fmt --check && cargo clippy -p dft-engine -p dft-physics --all-targets -- -D warnings` | 0 | Format and changed-crate all-target lints |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh` | 1 | Build/fmt/clippy PASS; 944 tests passed, 4 failed, 200 ignored; failures described below |

Raw local logs are `target/g2-check-final.log`, `g2-refscf-final2.log`,
`g2-physics.log`, `g2-scf-final.log`, and `g2-clippy.log` (the latter four
also under `target/`). They are uncommitted scratch, not published evidence.

The full workspace script additionally exposes an obsolete assertion at
`crates/dft-engine/tests/refstart.rs:519`: it explicitly enables functional
override, selects PBE and demands rejection. G2 intentionally makes that
start succeed. This file is outside OWNED_FILES; its assertion requires
master scope extension/reconciliation. Three unrelated CLI pseudo-cache
tests reject the sandbox-writable in-checkout scratch directory, because
those tests require a cache outside any repository. The script's `ps` load
probe is also sandbox-denied (nonfatal). Do not interpret these failures as
a passing whole-workspace gate.

Base revision: `f9cfd1a93798cf125d69b752921be7e34c5bedc7`. All changes remain
uncommitted. No shared registry/dashboard/README artifacts were edited.
