# P3-RESET: guarded ideal Davidson contraction

Status: default-off experiment, Mac correctness/work evidence only. Master
reviews and commits; `baandi-dev08` executes this request. No remote run,
install or QE source access occurred here. Derived from the own-words
[Davidson iteration note](../../../research/qe-study/davidson-iteration-structure-gamma.md),
implemented using this engine's existing rotations and factor adapters.

## Premise and implementation

Ticket-relayed si128 MPI8: 26 restarts, restart stage about 1.7 s, plus
measured refresh work attributed to border/factor. These are motivation,
not fresh measurements. The historical unguarded contraction/whitening
experiment drifted 1.0154099783221682e-11 Ha on the degenerate fixture. This
re-evaluation is explicitly authorized by P3-RESET under strict 1e-10 Ha /
standard 5e-9 Ha, convergence and verification gates, rather than bit identity.
No existing tolerance or default is relaxed.

`BAANDI_REF_DAV_IDEAL_RESET=off|on`, or `[performance] ref_dav_ideal_reset`,
default **off**. Only a successful Ritz contraction may reset S/T/C/Y to I
and reduced H/A to diag(theta). Repair/whitening still measures S/H, factors
S and transforms H. No H call or EVD is added by this experiment.

The on arm checks global column squared norms at each candidate contraction
(O(local rows * bands), one bands-wide sum). Nonfinite norms or |norm²-1|
above 1e-12 trigger measured refresh and disable ideal resets for the rest
of that solve. Every fourth contraction instead performs the complete
measured refresh and checks max|S-I|, including off-diagonal entries, with
the same threshold; an excessive defect also disables the experiment for
that solve. Thus no more than three consecutive contractions lack a full
measurement. The interval is a conservative experiment safeguard, not a
material-tuned performance choice or an error bound. The norm guard alone
cannot detect unit-norm correlated columns. Independent verification remains
necessary even if all guards pass. A fresh solve starts a fresh guard history.

Per-solve and per-iteration `Record.resets` record attempts, ideal resets,
periodic refreshes, guard trips, subsequent disabled refreshes, and the
maximum periodic measured defect. The existing run manifest entry for
`ref_dav_ideal_reset` includes cumulative successful-solve counts for the
calling thread/rank and run ID in `note.path`, even with diagnostics off;
`note.fallback` states a guard trip. The bounded thread-local summary resets
when a different run ID publishes; interleaved runs must use their complete
per-solve `Record.resets` instead. Failed solves do not publish partial counts. Root and replicated reduced placement follow the same guard decisions.
For a periodic guard trip, both periodic and guard counters increment:
measured restarts = attempts - ideal, not the sum of those two counters.

For an ideal reset, nominal saved Gram/transform work per rank is
4*r*n² + 4*n³ FLOPs plus factor work. Its norm check costs about 2*r*n+n
scalar FLOPs and 8*n reduction-input bytes; two 8*n² refresh sums disappear.
Compared with adopted LEAN_RESET, five additional n² matrices are written
(40*n² bytes). Periodic/disabled refreshes retain measured arithmetic. Actual
savings must use observed reset counts, not all restart counts.

## Exactly one profiled run per arm and case

| Case | Model | A: off | B: on |
| --- | --- | --- | --- |
| si64 | MPI8, one thread/rank | 1 | 1 |
| si128 | MPI8, one thread/rank | 1 | 1 |
| elongated | MPI8, one thread/rank | 1 | 1 |
| rskew | MPI8, one thread/rank | 1 | 1 |

**Eight executions total**, full reference-SCF convergence with independent
verification charged. No extra warmup, repetitions, A/A, medians, crossed
pairs, or replacement of failures in this request. One case/arm per 600 s
fast job, `p3-reset-<case>-<arm>-mpi8`; retain timeout/failure and exit code.
Do not truncate convergence or silently extend a job. A larger precision or
promotion campaign is separate work. A single profiled execution is a
mechanism screen, not uncertainty-qualified speed or default-adoption evidence.

Use the exact frozen P3l case/input/UPF/start and tight-reference manifests
(`/data/sig-dft/p3l/REPORT.md` and its referenced P3i/P3g/vc records). Resolve
actual paths and SHA-256 on dev08; case labels alone are insufficient. Keep
cell/coordinates, cutoff, FFT grid, bands, occupations, potential bytes,
start, threshold, verification and output policy identical. Pin all other
switches to the current accepted P3l manifest, including adopted Davidson
and streaming switches; do not combine a P3-ADDV candidate in this screen.
Keep REFSCF_VERIFY_ENERGY off. Use full diagnostics in both arms.

Build one reviewed candidate in the established native MPI/MKL environment,
without installs:

```sh
CARGO_TARGET_DIR=target/p3-reset cargo build --release -p dft-cli --features mpi,mkl
```

Record source commit/patch hash, compiler/MPI/MKL versions, native/build flags,
build exit, binary SHA-256, frozen inputs/UPF/start/reference hashes, complete
effective switch manifest, exact commands and node load. Serialize builds,
gates and measurements. Bind to socket-0 physical cores 0-7, one rank/core,
with actual affinity evidence; generic mapping flags alone are not proof.
Invocation skeleton (master resolves frozen inputs and unique output paths):

```sh
BAANDI_REF_DAV_IDEAL_RESET=off BAANDI_DAVIDSON_STAGES=on \
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-reference-input-with-unique-output-path>
```

Arm B changes only IDEAL_RESET to on and the output destination. Capture the
existing caller-resolved CPU profile during that same execution, including
rotation, border Gram, factor/transform, guard scans and collectives. Record
profiler/rate/sample counts; fewer than 30 relevant samples do not justify
quantitative attribution. A profiler failure does not authorize a rerun.
No new QE execution or QE speed ratio is requested.

## Required output

Retain stdout/stderr, exit status, external delivered wall (including failed
verification attempts/recovery), results.json and density/eigenvector hashes.
For each arm report SCF steps, verification attempts/pass/fail, total H
calls/columns, Davidson iterations/solves/restarts/repairs/dropped columns,
active-band and reduced-order traces, ideal/periodic/guard/disabled counts,
and the first changed branch/threshold operand. Changed counts are reported,
not silently rejected as bit-identity failures or hidden by rounded energy.

| case/arm | verified | E (Ha) | E - E_strict (Ha) | E_B - E_A (Ha) | H columns | iterations/restarts | ideal/periodic/guard/disabled | delivered s | verification s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| pending dev08 | | | | | | | | | |

Use an independently established strict reference for exactly the same
physical conditions, with command/version/input/UPF hashes. If missing,
report the missing oracle; do not invent E_strict or add an unauthorized run.
Strict class is 1e-10 Ha, standard class 5e-9 Ha. Neither Ritz stopping nor
orthogonality is a fixed-point certificate. Preserve convergence failures.

For `convergence.work.davidson_stages.work`, retain rank0 and per-rank values
and componentwise maxima separately:

| case/arm/rank/stage | calls | exclusive s | GEMM/scalar FLOPs | copy/zero/stream bytes | collective calls/bytes/s |
| --- | --- | --- | --- | --- | --- |
| restart | | | | | |
| border_s | | | | | |
| border_h | | | | | |
| factor | | | | | |
| transform | | | | | |

Also report coarse border/factor/restart seconds and the sum of affected
exclusive stages per rank. Coarse refresh transforms are charged to factor;
fine `transform` is separate. Do not sum componentwise rank maxima as a
critical path or double-count collective time. Full work traces show the
removed full n-by-n borders/refresh transforms. Use actual counts and
per-rank local rows for the work prediction above.

## Local attempt record

Base: `941a67fd68f9a5619a8476ea9206b747dfbaf887`, plus the uncommitted patch.
Local evidence below is fixed-H numerical correctness, not SCF total-energy
or dev08 speed acceptance. The central ledger/log, dashboard and generated
README evidence are outside this ticket's OWNED_FILES; master may transfer
this record after review. No new external-reference claim is made.

Local targeted command (debug, under 1 s test execution after compilation):

```sh
cargo test -p dft-solvers --lib ideal_reset -- --nocapture
```

Two tests pass: 24 fixture A/B pairs plus injected norm/correlation failures
on ranks 1/3/8. Four fixtures x three rank counts x two placements; all 48
solves converge. The on/off active-band/order/repair/drop/restart decisions,
H calls and H columns match for every pair. Tests use the original DCT
fixture construction from `tests/ref_davidson.rs`, not a silicon-only matrix.
Max root difference on/off: **2.886579864025407e-15 Ha**. Max doubly occupied
band-sum difference (six occupied roots for the 32-row fixtures, three for
the seven-row fixture): **9.325873406851315e-15 Ha**. Max on/dense-EVD root
error: **7.838174553853605e-14 Ha**. These pass 1e-10 Ha without relaxing any
existing test. Fresh H applications agree with returned H X within 2e-11.

The test-only Gram audit prints `orth_off` / `orth_on` for **every restart**,
including before a measured fallback, with full precision. It is absent
from production; audit collectives must not be used as production work
counts. Scratch transcript: `target/tmp/reset-test-final.log`; rerun the
command above to regenerate it after cleanup. Root/replicated results and
counts match; the table shows each rank-count result once.

| Spectrum, dimension/bands/capacity | ranks | H columns off/on | iterations off/on | restarts off/on | ideal/periodic/guard/disabled | max restart defect on |
| --- | --- | --- | --- | --- | --- | --- |
| spaced, 32/8/16 | 1 | 233/233 | 41/41 | 29/29 | 0/0/1/28 | 3.507e-12 |
| spaced, 32/8/16 | 3 | 233/233 | 41/41 | 29/29 | 0/0/1/28 | 1.692e-12 |
| spaced, 32/8/16 | 8 | 233/233 | 41/41 | 29/29 | 22/7/0/0 | 6.227e-13 |
| degenerate, 32/8/16 | 1 | 237/237 | 48/48 | 30/30 | 0/0/1/29 | 1.388e-11 |
| degenerate, 32/8/16 | 3 | 237/237 | 48/48 | 30/30 | 0/0/1/29 | 1.846e-12 |
| degenerate, 32/8/16 | 8 | 237/237 | 48/48 | 30/30 | 0/0/1/29 | 3.648e-12 |
| degenerate, 32/8/24 | 1 | 115/115 | 28/28 | 6/6 | 5/1/0/0 | 1.554e-15 |
| degenerate, 32/8/24 | 3 | 115/115 | 28/28 | 6/6 | 5/1/0/0 | 5.107e-15 |
| degenerate, 32/8/24 | 8 | 115/115 | 28/28 | 6/6 | 5/1/0/0 | 3.442e-15 |
| small, 7/4/7 | 1 | 7/7 | 2/2 | 1/1 | 1/0/0/0 | 6.661e-16 |
| small, 7/4/7 | 3 | 7/7 | 2/2 | 1/1 | 1/0/0/0 | 4.441e-16 |
| small, 7/4/7 | 8 | 7/7 | 2/2 | 1/1 | 1/0/0/0 | 4.441e-16 |

Global maximum restart defect in either arm: **1.3880785409980945e-11**,
dimensionless, before fallback. The 2n degenerate fixture's first norm guard
trips at all tested rank counts, so this candidate does not re-establish
unguarded reset accuracy there. The 3n fixture exercises actual ideal resets
on the same degenerate spectrum. The conservative guard may prevent savings
on other problems; report its dev08 rate rather than changing it to obtain
a positive timing. Band energies alone do not establish SCF verification.


### Synthetic full-SCF verification (additional Mac probe)

The standard `refscf` integration tests build their own switch configuration;
setting the environment cannot exercise the new on arm in those tests. A
temporary probe therefore uses an explicit RunConfig for each arm, calling
`driver::run_with` on the existing synthetic atomic-start fixture. It changes
no crate dependencies or engine tests. Rustc 1.99.1, arm64, pure Rust backend.
Three two-atom cells, five bands, 2 Ha cutoff, broad Gaussian synthetic NC
potential with a nonlocal projector, default reference profile and independent
fixed-point verification (`refscf_verify_energy=off`). All 18 runs (3 shapes x
3 simulated rank counts x 2 arms) converge and verify. Per-rank-count counts
match between arms. No timing conclusion is drawn.

| cell | ranks | E_on - E_off (Ha) | maps off/on | Davidson iterations off/on | H columns off/on | ideal/periodic/guard/disabled |
| --- | --- | --- | --- | --- | --- | --- |
| cubic | 1 | 1.7319479184152442e-14 | 21/21 | 96/96 | 356/356 | 31/6/0/0 |
| cubic | 3 | -7.5495165674510645e-15 | 21/21 | 96/96 | 356/356 | 31/6/0/0 |
| cubic | 8 | -2.7533531010703882e-14 | 21/21 | 96/96 | 356/356 | 31/6/0/0 |
| skew | 1 | -2.2204460492503131e-16 | 13/13 | 47/47 | 237/237 | 21/4/0/0 |
| skew | 3 | 0 | 13/13 | 47/47 | 237/237 | 21/4/0/0 |
| skew | 8 | -6.6613381477509392e-16 | 13/13 | 47/47 | 237/237 | 21/4/0/0 |
| elongated | 1 | 8.8817841970012523e-16 | 9/9 | 43/43 | 169/169 | 19/2/0/0 |
| elongated | 3 | -4.4408920985006262e-16 | 9/9 | 43/43 | 169/169 | 19/2/0/0 |
| elongated | 8 | -4.4408920985006262e-16 | 9/9 | 43/43 | 169/169 | 19/2/0/0 |

Maximum total-energy change: **2.7533531010703882e-14 Ha**, below 1e-10 Ha.
This is on/off agreement plus independent fixed-point verification on synthetic
systems, not a new QE/material reference or an E_strict measurement. The
cumulative manifest counts are actual probe output; no guard trips occurred.

### Reproduce the SCF probe without editing an engine test

After `cargo test -p dft-engine --test refscf`, save this Rust block as
`target/tmp/reset_scf_probe.rs` and the following Python block as
`target/tmp/run_reset_scf_probe.py`. The Python runner finds the exact numerics
artifact linked by the engine using Cargo fingerprints (separate test targets
can produce several numerics artifacts), compiles with the pinned rustc and
runs the probe. Run from the worktree root. Logs: `target/tmp/reset-scf.log`.

```rust
#[allow(dead_code)]
#[path = "../../crates/dft-engine/tests/common/mod.rs"] mod common;
use dft_engine::{ScfOptions, GammaReal, SpeciesPotential, refscf::{driver, profile::Profile}, refstart::StartInput};
use dft_geometry::{Cell, Atom};
use dft_numerics::switches::{RunConfig, Id};
use dft_mpi::Collective;
fn main() {
    let pp = common::potential_with_atomic_start(2.0, 1.8, true, 161, "SLA PZ");
    let atoms = [Atom::new("Xx", [0.13,0.21,0.32]), Atom::new("Xx", [0.57,0.68,0.71])];
    let potentials = [SpeciesPotential { label: "Xx", potential: &pp }];
    let options = ScfOptions { cutoff_hartree: 2.0, bands: Some(5), max_iterations: 100, gamma_real: GammaReal::On, ..Default::default() };
    for (label, lattice) in [
        ("cubic", [[7.,0.,0.],[0.,7.,0.],[0.,0.,7.]]),
        ("skew", [[7.,0.,0.],[0.2,7.,0.],[0.,0.1,7.]]),
        ("elongated", [[5.,0.,0.],[0.,5.,0.],[0.,0.,10.]])
    ] {
        let cell=Cell::from_user_rows_bohr(lattice).unwrap();
        for ranks in [1,3,8] {
            dft_mpi::simulate(ranks, |c| {
                let mut baseline = None;
                for arm in ["off", "on"] {
                    let config=RunConfig::from_pairs(&[("ref_dav_ideal_reset",arm),("refscf_verify_energy","off")]).unwrap();
                    let _scope=config.enter();
                    let run=driver::run_with(StartInput {cell:&cell,atoms:&atoms,potentials:&potentials,options:&options},Profile::default(),c).unwrap();
                    assert!(run.outcome.converged);
                    assert!(run.outcome.result.verification.as_ref().unwrap().passed());
                    let e=run.outcome.result.energy.total_hartree();
                    let base=*baseline.get_or_insert(e);
                    assert!((e-base).abs()<1e-10, "{label} ranks={ranks} delta={:e} Ha",e-base);
                    let calls:Vec<_>=run.records.iter().flat_map(|r|&r.calls).chain(&run.verification_calls).collect();
                    if c.rank()==0 {
                        println!("RESET_SCF case={label} ranks={ranks} arm={arm} E_Ha={e:.17e} delta_Ha={:.17e} maps={} Dav={} H={} verified=true resets={:?}",e-base,run.records.len(),calls.iter().map(|c|c.iterations).sum::<usize>(),calls.iter().map(|c|c.h_columns).sum::<usize>(),config.path_note(Id::REF_DAV_IDEAL_RESET));
                    }
                }
            });
        }
    }
}
```

```python
from pathlib import Path
import json
import subprocess

root = Path.cwd()
deps = root / "target/debug/deps"
engine = max(deps.glob("libdft_engine-*.rlib"), key=lambda p: p.stat().st_mtime)
fingerprints = root / "target/debug/.fingerprint"
info = json.loads((fingerprints / ("dft-engine-" + engine.stem.split("-")[1]) /
                   "lib-dft_engine.json").read_text())
num_hash = next(d[3] for d in info["deps"] if d[1] == "dft_numerics")
num = next(p.parent.name.split("dft-numerics-")[1]
           for p in fingerprints.glob("dft-numerics-*/lib-dft_numerics")
           if int.from_bytes(bytes.fromhex(p.read_text()), "little") == num_hash)
cmd = ["rustc", "--edition=2024", "-L", "dependency=" + str(deps),
       str(root / "target/tmp/reset_scf_probe.rs"), "-o",
       str(root / "target/tmp/reset-scf-probe")]
for name in ["engine", "numerics", "geometry", "physics", "upf", "mpi"]:
    lib = (deps / ("libdft_numerics-" + num + ".rlib") if name == "numerics"
           else engine if name == "engine"
           else max(deps.glob("libdft_" + name + "-*.rlib"), key=lambda p: p.stat().st_mtime))
    cmd += ["--extern", "dft_" + name + "=" + str(lib)]
subprocess.run(cmd, check=True)
with (root / "target/tmp/reset-scf.log").open("w") as log:
    subprocess.run([str(root / "target/tmp/reset-scf-probe")],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
```


### Verification commands and limits

All commands below were executed on this Mac; no command is inferred passing.

| Command | Exit | Result |
| --- | --- | --- |
| `cargo check --workspace --tests` | 0 | Final workspace and test targets compile. |
| `cargo test -p dft-solvers --lib ideal_reset -- --nocapture` | 0 | 2 tests; 24 fixed-H A/B pairs, every-restart audit and injected guards. |
| `cargo test -p dft-solvers` | 0 | 180 passed, 25 ignored; existing tolerances unchanged. |
| `cargo test -p dft-engine --test refscf` | 0 | 28 passed, 1 ignored; baseline regression, not an on-arm test. |
| `cargo test -p dft-numerics switches` | 0 | 14 passed; registry contains 159 switches. |
| `cargo fmt --check` | 0 | Formatting passes. |
| `cargo clippy -p dft-solvers -p dft-numerics --all-targets -- -D warnings` | 0 | No warnings. |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh fast dft-solvers dft-numerics` | 0 | Build/test/fmt/clippy PASS; 269 passed, 47 ignored across 23 binaries. |
| `python3 target/tmp/run_reset_scf_probe.py` | 0 | Exact documented probe compiles and all 18 full SCFs converge/verify. |
| `git diff --check` | 0 | No whitespace errors. |

The gate's optional load-report `ps` invocation was denied by the sandbox;
its four actual gates completed successfully. No full workspace execution,
release/MKL/real-MPI acceptance or dev08 profiling was done. Ignored tests
remain ignored. The task-owned fast gate plus the requested workspace check
and targeted commands are the local validation scope. README/datasets,
central ledger/log and session/dashboard updates are left to the master
because they are outside OWNED_FILES. Changes remain uncommitted.
