# P3-A2A: padded transpose accounting and stick-cap experiment

Base: `ee6a0c40d7ce8c197e570681247bc6744ff45c67`; changes uncommitted.
Default-off candidate, no measured dev08 speed claim. Master executes through
`baandi-dev08`; implementer makes no remote changes. Registry/dashboard,
README generation, shared ledger/log updates and commit are master integration
work outside this ticket's OWNED_FILES. The local registry has no P3-A2A run
or `.delegate/P3-A2A/report.json`; this report does not invent registry state.
No recurring orchestrator timer or dashboard upload is active in this implementer.

## Finding and mechanism

[Exact accounting](../../../design/reference-scf.md#p3-a2a-payload-minimum-and-capped-stick-ownership)
shows that sphere-only sticks and whole-H-block batching already exist.
Two real bands packed into complex coefficients need both mirror sticks.
The largest removable component in the existing equal-block interface is
rank padding, not half of the physical payload.

`BAANDI_FFT_A2A_BALANCE=on` minimizes the maximum sticks/rank while retaining
mirror units. It keeps the original owner assignment unless the equal block
strictly shrinks. Pair units are assigned before singletons, decreasing
coefficient weight first, choosing the least coefficient-loaded rank with
capacity. The closed-form capacity bound is geometry independent; there is
no material/threshold tuning. Collective calls and arithmetic within each
FFT do not change. Metadata is construction-only and existing profile/counter
paths measure actual reduced send buffers.

This changes row owners in active cases: nonlocal and solver reductions are
rounding class. Possible tradeoff: more coefficient-load imbalance despite
less FFT-work imbalance. Report both, plus trajectory/count changes. If the
block is already optimal, the original row owners and arithmetic remain.

**Negative result for the main motivation:** catalog Si64 and Si128 MPI8 are
already optimal under this interface/mirror constraint. Off/on payload per
pair/rank/direction is respectively 145152/145152 and 290304/290304 bytes.
All-padding-removal ceilings are only 0.2976% and 0.5456% respectively.
The switch is expected to be inactive there if remote inputs match these
catalog spheres. Do not advertise this as a Si128 speedup or tune a partition
to force an effect. Elongated/skew support may benefit; the exact remote
supports were not available locally. The eight requested executions remain
a screen including the two expected no-change controls, not a timing claim.

No async wrapper is added. The ticket aggregate cannot separate alltoall
startup, bandwidth, arrival wait and projector reduction time. Si128 alone
does not yet establish the prior asynchronous-overlap reopening premise.
The accounting section lists the missing wrapper/lifetime/progress support.

## Exactly eight profiled executions

The explicit ticket protocol overrides historical repeated-run recipes:
**ONE profiled run per arm**, no warmup executions, repetitions, A/A series,
medians or additional QE jobs. Use one reviewed binary for both arms.

| Case | Model | A | B |
| --- | --- | --- | --- |
| si64 | MPI8, one thread/rank | FFT_A2A_BALANCE=off, once | FFT_A2A_BALANCE=on, once |
| si128 | MPI8, one thread/rank | off, once | on, once |
| elongated | MPI8, one thread/rank | off, once | on, once |
| rskew | MPI8, one thread/rank | off, once | on, once |

Use the current P3 converging reference-driver input manifests from
`/data/sig-dft/vc/` and the existing P3j campaign, not newly generated nominal
cells. Freeze cell/coordinates, cutoff, FFT grid, bands, UPF bytes, XC,
start, convergence thresholds, verification and output policy identically.
Record actual absolute paths and input/UPF/start SHA-256. Preserve the
adopted defaults identically; force M2_STICKS=mirror and M2_BATCH=on in both
arms. Do not reuse row-local snapshots between changed layouts: initialize
both arms from the same physical start. This ownership/rounding change calls
for full convergence, including verification, rather than a fixed-row replay.

Build after master review with the established baseline native flags:

```sh
CARGO_TARGET_DIR=target/p3-a2a cargo build --release -p dft-cli --features mpi,mkl
```

Record source commit/patch hash, binary SHA-256, compiler/MPI/MKL versions,
features and build exit. One case/arm per 600 s fast job; any failed or timed
out run is retained, without automatic replacement or a fabricated convergence
claim. Use the established socket-0 physical cores 0-7 placement, no overlapping
builds/jobs, and record actual rank binding and load. `--map-by core` alone
is not proof of placement. Command skeleton (resolve existing paths before
queuing, separate output destinations):

```sh
BAANDI_FFT_A2A_BALANCE=off BAANDI_M2_STICKS=mirror BAANDI_M2_BATCH=on \
BAANDI_DAVIDSON_STAGES=on BAANDI_PHASE_PROFILE=on BAANDI_FFT_PROFILE=on \
BAANDI_STAGE_TIMERS=on OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
mpirun -np 8 --map-by core --bind-to core --report-bindings \
  <reviewed-binary> run <frozen-reference-driver-input-for-arm-A>
```

Arm B changes only FFT_A2A_BALANCE and output destination. Retain the same
reference-driver selection as the established campaign. If the CLI input
selects switches, explicitly check the resolved manifest and agreement hash.
Detailed profiling is on in both arms and these numbers do not belong in a
QE speed table. No mixed MPI/threads comparison is requested.

## Required evidence, including no-change controls

- Complete stdout/stderr, results JSON, exit code, resolved switch manifest,
  run completion marker, source/binary/input/UPF hashes, binding and load.
- Per rank and H call: band/pair width, grid, first axis, S, owned stick and
  coefficient/ORC row counts, plane counts, maximum counts and effective
  smaller-block/fallback decision. Counts may be derived outside timing from
  the exact ordered basis with `position_owners`; never log inside hot loops.
  For a nonempty plan, `s_max = bytes/(16*q*R*p_max)` independently verifies
  the padded block from the existing stage counters. Report coefficient
  imbalance as well as the reduction in maximum stick count.
- H-scoped alltoall **calls, outgoing bytes and seconds**, per rank, inverse
  and forward separately. Use FFT stages 3 and 13 for transpose time; the
  reference work collective time includes tree sums and is not an alltoall
  timer. Keep projector sums and any mirror exchanges separate. Check two
  transposes/H call, zero mirror exchanges and one projector sum on the
  mirror/nonlocal path. Density transposes are outside H accounting.
- FFT split, per rank: first/middle/last 1D kernels, coefficient scatter,
  transpose pack/unpack, middle pack/unpack, slab gather/scatter, multiply,
  final gather/scale; bytes and calls alongside seconds. Do not add unrelated
  rank maxima as wall time. Existing arrival sites 104/114 can diagnose waits
  if COLLECTIVE_TRACE is enabled identically in both arms; no extra run is
  authorized for it. Do not infer pure startup latency from residence time.
- Full delivered wall, H wall, SCF maps, Davidson iterations, H columns,
  verification attempts/status, energies in Ha, density differences/hash,
  eigenvalue differences and all failures. On fallback cases require identical
  science/counts. On changed layouts retain the existing strict 1e-10 Ha /
  standard 5e-9 Ha accuracy classes and the 1e-12 Ha unique-solution regression
  check; do not relax any tolerance. Report any count/trajectory change.

A single instrumented result is directional, not a precision measurement.
No default promotion is part of this request. A regression or neutral result
is retained. A byte saving with worse solve time is not adopted silently.

## Local validation

New tests pass on pure Rust/Mac (no native MPI or MKL speed measurement):

- `a2a_balance_preserves_fft_bits_and_minimizes_padding`: seven supports,
  ranks 1/2/3/8, greedy/mirror, widths 5/1/0/3/5, odd/prime grids, empty and
  G=0-only supports, inverse planes, full local application and kept planes.
  Exact bit equality, stable warmed send/receive pointers and exact profile
  bytes. Eleven layouts reduce payload; examples below. Optimal-capacity and
  ownership/fallback/mirror invariants are independently asserted.
- `a2a_catalog_payload_accounting`: 23871/47851 PW, 1005/2005 sticks, exact
  Si64/128 catalog payloads and unchanged on/off owners as tabulated above.
- `a2a_balance_h_equivalence`: cubic/skew/elongated with grids 12^3,
  11x13x17 and 7x11x23, ranks 1/2/3/8, both stick policies, bands
  16/1/2/5/0, batched/unbatched, no projectors or mixed 1/3/5 channels with
  off-diagonal couplings, direct packing/sparse clearing enabled. Six layouts
  change ownership; local+kinetic terms are bit-identical, full H maximum
  absolute difference **1.1102230246251565e-16 Ha** (gate 1e-12).

| Synthetic support, ranks, policy | Outgoing bytes/pair/rank/direction off | on | Reduction |
| --- | ---: | ---: | ---: |
| 16^3 sphere, 8, mirror | 2048 | 1792 | 12.5% |
| 15x15x45 sphere, 3, mirror | 12480 | 12000 | 3.846% |
| 12x14x16 skew support, 3, mirror | 5040 | 4800 | 4.762% |
| 5x7x11 sphere, 3, mirror | 576 | 480 | 16.667% |
| 7x11x13 sphere, 8, greedy | 640 | 512 | 20% |

These are exact byte counts, not speed measurements. No largest-band or
silicon-only threshold is introduced.

The standalone reproduction below compares complete verified reference SCF
on/off at ranks **1/2/3/8**, using the existing synthetic atomic-start fixture
with nonlocal projectors, physical start and all other settings identical.
All ranks agree within each run; all 12 on/off comparisons retain counts and
pass verification. Maximum energy difference **1.56097357262297e-13 Ha**
(gate 1e-12 Ha). These are synthetic regression fixtures, not QE references.

| Fixture | SCF maps, both arms at every rank count | Davidson iterations | H columns | Verification attempts |
| --- | ---: | ---: | ---: | ---: |
| cubic | 21 | 96 | 356 | 1 |
| skew | 13 | 47 | 237 | 1 |
| elongated | 9 | 43 | 169 | 1 |

Requested commands: workspace check, physics, MPI, reference SCF (24 passed,
1 ignored), fmt and final clippy pass (exit 0). `cargo test -p dft-solvers`
returns 101: an existing `ref_davidson::diagnostic_tests::
davrem_one_correction_has_analytic_row_work_and_error_cleanup` assertion
expects 36 hypot calls, but the already-adopted REF_DAV_SQRT default gives 0.
The identical failure was reproduced from an unmodified `git archive HEAD`
under `target/tmp/a2a-baseline` (exit 101). This patch does not edit that
out-of-scope test or revert the adopted default. The unskipped check.sh result
is exit 1: build/fmt/clippy PASS; test FAIL (278 passed, 1 pre-existing
failure, 43 ignored, 23 binaries). No passing full-gate claim is made.

Development failures retained: the first new FFT-counter test omitted
PHASE_PROFILE and read zero counters; enabling the required profile scope
fixed the test. Initial clippy flagged complex test result types; named type
aliases fixed them. Final focused tests and clippy pass. Neither change
relaxes a numerical tolerance.


### Reproduce the reference-SCF on/off count check

Keep the existing engine tests unchanged (outside ownership). After
`cargo test -p dft-engine --test refscf`, save this block as
`target/tmp/a2a_refscf.rs` and compile against the just-built workspace rlibs.
No new dependency, potential file, network access or production code is used.

```rust
#![allow(dead_code)]
#[path = "../../crates/dft-engine/tests/common/mod.rs"]
mod common;
use dft_engine::refscf::profile::Profile;
use dft_numerics::switches::RunConfig;
fn main() {
    for (label, rows) in [
        ("cubic", [[7., 0., 0.], [0., 7., 0.], [0., 0., 7.]]),
        ("skew", [[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]),
        ("elongated", [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]]),
    ] {
        let pp = common::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
        let cell = dft_geometry::Cell::from_user_rows_bohr(rows).unwrap();
        let atoms = [
            dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32]),
            dft_geometry::Atom::new("Xx", [0.57, 0.68, 0.71]),
        ];
        let potentials = [dft_engine::SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let options = dft_engine::ScfOptions {
            cutoff_hartree: 2.,
            bands: Some(5),
            max_iterations: 100,
            gamma_real: dft_engine::GammaReal::On,
            ..Default::default()
        };
        for ranks in [1, 2, 3, 8] {
            let run = |arm| {
                dft_mpi::simulate(ranks, |comm| {
                    let cfg = RunConfig::from_pairs(&[
                        ("fft_a2a_balance", arm),
                        ("m2_sticks", "mirror"),
                        ("refscf_verify_energy", "off"),
                    ])
                    .unwrap();
                    let _scope = cfg.enter();
                    let result = dft_engine::refscf::driver::run_with(
                        dft_engine::refstart::StartInput {
                            cell: &cell,
                            atoms: &atoms,
                            potentials: &potentials,
                            options: &options,
                        },
                        Profile::default(),
                        comm,
                    )
                    .unwrap();
                    assert!(result.outcome.converged);
                    let verify = result.outcome.result.verification.as_ref().unwrap();
                    assert!(verify.passed());
                    let calls: Vec<_> = result
                        .records
                        .iter()
                        .flat_map(|r| &r.calls)
                        .chain(&result.verification_calls)
                        .collect();
                    (
                        result.outcome.result.energy.total_hartree(),
                        result.records.len(),
                        calls.iter().map(|c| c.iterations).sum::<usize>(),
                        calls.iter().map(|c| c.h_columns).sum::<usize>(),
                        verify.attempts.len(),
                    )
                })
            };
            let off = run("off");
            let on = run("on");
            for arm in [&off, &on] {
                for rank in &arm[1..] {
                    assert_eq!(arm[0], *rank);
                }
            }
            let a = off[0];
            let b = on[0];
            let delta = (a.0 - b.0).abs();
            println!("A2A_REFSCF {label} P={ranks} off={a:?} on={b:?} delta_ha={delta:e}");
            assert!(delta < 1e-12);
            assert_eq!((a.1, a.2, a.3, a.4), (b.1, b.2, b.3, b.4));
        }
    }
}
```

From the workspace root, the following compiles and runs it (the most recent
rlib of each crate must be from this checkout, not a concurrent baseline build):

```sh
python3 - <<'PY'
from pathlib import Path
import subprocess
names = ['dft_engine', 'dft_geometry', 'dft_physics', 'dft_upf',
         'dft_mpi', 'dft_numerics']
args = ['rustc', '--edition=2024', '-O', 'target/tmp/a2a_refscf.rs',
        '-L', 'dependency=target/debug/deps', '-o', 'target/tmp/a2a_refscf']
for name in names:
    libs = list(Path('target/debug/deps').glob('lib' + name + '-*.rlib'))
    lib = max(libs, key=lambda p: p.stat().st_mtime)
    args += ['--extern', name + '=' + str(lib)]
subprocess.run(args, check=True)
subprocess.run(['target/tmp/a2a_refscf'], check=True)
PY
```


### Final command record

| Command | Exit | Result |
| --- | ---: | --- |
| `cargo check --workspace --tests` | 0 | Workspace/tests compile |
| `cargo test -p dft-solvers --lib a2a_ -- --nocapture` | 0 | 2 new FFT/accounting tests pass |
| `cargo test -p dft-physics --lib a2a_balance -- --nocapture` | 0 | New H equivalence test passes |
| `cargo test -p dft-solvers` | 101 | Existing hypot-count assertion; new FFT tests pass |
| `cargo test -p dft-physics` | 0 | 96 passed, 18 ignored across test targets |
| `cargo test -p dft-mpi` | 0 | 5 passed |
| `cargo test -p dft-engine --test refscf` | 0 | 24 passed, 1 ignored |
| `cargo fmt --check` | 0 | Pass |
| `cargo clippy -p dft-solvers -p dft-physics -p dft-mpi --all-targets -- -D warnings` | 0 | Final pass after test type aliases |
| `cargo test -p dft-numerics --test all switches::legacy_defaults_are_unchanged -- --exact` | 0 | Updated registry count passes |
| `cargo test -p dft-cli --test all run::expressions_preserve_energy_and_record_resolutions -- --exact` | 0 | Updated switch manifest count passes |
| `target/tmp/a2a_refscf` (compiled as above) | 0 | 12 complete verified on/off SCF comparisons |
| `CHECK_SCCACHE=0 CHECK_NICE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh fast dft-solvers dft-physics dft-mpi` | 1 | Build/fmt/clippy pass; only the pre-existing test fails |
| Baseline `CARGO_TARGET_DIR="$PWD/target" cargo test --manifest-path target/tmp/a2a-baseline/Cargo.toml -p dft-solvers --lib ref_davidson::diagnostic_tests::davrem_one_correction_has_analytic_row_work_and_error_cleanup -- --exact` | 101 | Same 0-versus-36 failure at base HEAD |

Scratch logs are `target/a2a-*.log` and `target/a2a-validation.json`; they are
not committed scientific evidence. The persistent tests and SCF reproduction
above supply the rerun procedure. No remote speed, MPI ABI, MKL equivalence or
production multi-system promotion is established by these Mac checks.
