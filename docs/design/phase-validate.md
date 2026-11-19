# Owned projector phase validation (G-PHASE-VALIDATE)

Status: implemented, **default off**, pending same-model dev08 A/B. This is a
numerical-range optimization of setup validation, not a physics approximation.
Parent: `d6be1d0cab87cb0654d6b6de2c45f8c9a1fe6cfb`; uncommitted candidate.
Task: `docs/performance/plan-2026-11-19-followup.md`, section 4.2.

## Contract and proof

`BAANDI_PHASE_VALIDATE=on` (or `[performance] phase_validate="on"`) changes only
the owned constructor's finite-dot prevalidation. The default is `off`.
At the original validation boundary, after the species' factor/norm checks:

1. Scan every Cartesian momentum component once per constructor, retaining
   whether all are finite and their maximum absolute value is at most 2^500.
2. Scan each species' position components once, with the same sufficient test.
3. If both tests pass, skip that species' atom-by-momentum finite-dot loop.
   Otherwise execute the original loop, in atom/momentum order, including its
   nonfinite-position predicate and first failing pair.

2^500 is exactly representable in binary64. Each product has magnitude at most
2^1000, and the intermediate/final sums are below 2^1002, far below overflow.
The proof allows underflow and signed zero: it proves finiteness only. It never
supplies a rounded dot product or phase. The separate `phase()` computation,
complex multiplication, factor/norm checks and order, lazy recipe and table,
owned row selection, and MPI collective sequence are unchanged. No collective,
allocation, logging or I/O is added to validation. Failed sufficient tests do
**not** reject input: huge imbalanced or cancelling values can still pass the
original predicate. Empty arrays retain the original vacuous acceptance.
Species error order is retained; proof scans do not return errors. On a failing
species, only transient atom metadata construction is deferred until validation
has succeeded; no previously observable output is changed.

The engine call site remains restricted to eligible owned Gamma states. Ordinary
one-core/threads-only construction uses the eager constructor and never executes
this prevalidation. Enabling the switch does not change that dispatch.

## Instrumentation and units

`RunConfig` exposes these additive counters, also serialized by the existing CLI
`run.counters` writer. They are local to a run and, in real MPI, to a process/rank.
A simulated-rank caller must enter a separate run scope for each rank.

| Counter | Definition |
| --- | --- |
| `phase_validate.constructors` | Owned constructors reaching validation setup, including later failures |
| `phase_validate.calls` | Species validation calls after their factors have passed |
| `phase_validate.pair_checks` | Actual original-loop dot checks, including the first failing pair |
| `phase_validate.fast_accepts` | Species calls accepted by the sufficient range proof |
| `phase_validate.fallbacks` | Enabled-proof declines; off does not count as a fallback |
| `phase_validate.nanoseconds` | Sum of exclusive validation spans, only with `BAANDI_PHASE_PROFILE=on` |

The first species' timer includes the momentum scan; subsequent species reuse
its verdict. All position scans and fallback loops are included. Factor checks,
radial tables, phase evaluation, atom metadata allocation, materialization and
counter publication are excluded. The span includes its scalar bookkeeping and
elapsed-clock read. No per-pair clock or atomic operation is performed. Counters
are published once on constructor exit, including ordinary error returns.
Profile off reads no clock and records zero nanoseconds, while counts remain
available. This instrumentation is not a new global phase-profile bucket.
The ordinary CLI writes root's counters; a per-rank diagnostic must read each
rank's `RunConfig` explicitly. No extra reporting collective was inserted.

## Mac measurements, 2026-11-19

macOS 27.0.1 (26A434), arm64, pure Rust release, one worker per rank. CPU model
and process inventory were unavailable under the sandbox (`sysctl`/`ps` denied).
These are local observations, not dev08 adoption evidence or QE comparisons.
No MPI launcher/compiler was found on PATH or standard Homebrew MPI locations;
real two-rank MPI was **not run**, consistent with ADR 0041's Mac policy.

Si64: cubic side 20.522425713435865 Bohr, 64 atoms, full cutoff 15 Ha,
coarse cutoff 7.5 Ha (fraction 0.5), **23,871 + 8,409** momenta. Instrumentation
confirms **2 constructors, 2 species calls, 2,065,920 pair checks per owned rank**.
On: **0 pair checks, 2 fast accepts, 0 fallbacks**. Counts are identical in every
measured repetition. No pseudopotential is stored in the repository.

External `../pseudo/Si.pz-vbc.UPF` SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

### Isolated constructor A/A and crossed A/B

The committed ignored test `hamiltonian::phase_validation_tests::phase_validation_si64`
builds both owned operators with **all rows selected**, sequentially on one
thread. It excludes UPF parsing, basis construction and index generation.
The constructor denominator includes allocation, phases and destruction, but
**is not the engine's full initialization**. It verifies that no lazy full
table is materialized. Warm-up sample 0 is discarded; each of 12 further
samples runs `A0,A1,B0,B1,A2`, with A=off and B=on. A0/A1 provides the noise
floor; A1/B0 and B1/A2 provide 24 crossed comparisons. The noise definition
was fixed before inspection: sample SD of the 12 A1-minus-A0 differences;
resolved means median paired saving at least three times that SD.

| Quantity | Seconds |
| --- | ---: |
| Off validation median (24 arms) | 0.000780146 |
| On validation median (24 arms) | 0.0000717085 |
| Paired saving median | **0.0007075205** |
| Paired saving mean | 0.000712133625 |
| A/A noise SD | 0.0000703289222 |
| Three times noise | 0.000210986767 |
| Off constructor wall median | 0.020770188 |

Saving/noise = **10.0602**, so the isolated saving is resolved by the requested
criterion. Median off validation share of constructor wall = **3.77028%**;
this is an upper bound on what removing this scan could recover in that
constructor, not a whole-run speedup. No total-wall gain is claimed.

All retained exclusive spans (seconds), sufficient to recompute the result:

| Sample | A0 | A1 | B0 | B1 | A2 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.000754375 | 0.000787875 | 0.000070042 | 0.000069958 | 0.000787501 |
| 2 | 0.000791251 | 0.000758875 | 0.000078625 | 0.000076000 | 0.000818292 |
| 3 | 0.000821000 | 0.000885834 | 0.000072125 | 0.000074708 | 0.000827042 |
| 4 | 0.000740001 | 0.000736501 | 0.000068583 | 0.000068708 | 0.000738166 |
| 5 | 0.000761042 | 0.000783125 | 0.000073750 | 0.000070834 | 0.000765625 |
| 6 | 0.000763166 | 0.000828916 | 0.000076125 | 0.000078625 | 0.000824209 |
| 7 | 0.000826625 | 0.000813499 | 0.000076042 | 0.000076166 | 0.000818708 |
| 8 | 0.000926334 | 0.000741958 | 0.000068999 | 0.000070083 | 0.000739583 |
| 9 | 0.000767875 | 0.000758043 | 0.000072126 | 0.000071501 | 0.000777167 |
| 10 | 0.000853667 | 0.000757458 | 0.000069916 | 0.000071916 | 0.000820833 |
| 11 | 0.000796125 | 0.000738291 | 0.000068876 | 0.000070167 | 0.000753833 |
| 12 | 0.000765084 | 0.000766708 | 0.000070542 | 0.000074584 | 0.000802166 |

A baseline-only run was collected before implementing the proof; it already
confirmed 2,065,920 checks. The crossed experiment uses the final production
code with switch off/on in one binary, avoiding cross-build timing subtraction.

### Engine construction denominator and rank scope

An ephemeral harness in `target/scratch/phase-harness` reads the exact committed
`benchmarks/cases/si64-gamma.toml` geometry and creates `ScfProblem::new` at
15 and 7.5 Ha, `ScfOptions::default()` plus `threads=1`. Each arm runs 13
samples, discarding sample 0. The denominator is the elapsed time of both
problem constructors/destructors. It excludes CLI parsing, lazy atomic-start
Hamiltonian work, SCF, verification and output; repeated Ewald memoization is
retained. It is a **warm problem-init denominator**, not the entire profile
`init` bucket, which also includes deferred atomic-start work.

| Execution | Mode | Rank | Validation s (median) | Problem init s (median) | Median validation/init |
| --- | --- | ---: | ---: | ---: | ---: |
| Real one-core serial | off | 0 | 0 | 0.101990584 | **0%** |
| Real one-core serial | on | 0 | 0 | 0.100734104 | **0%** |
| Two simulated ranks | off | 0 | 0.0009547505 | 0.079239688 | 1.21008% |
| Two simulated ranks | off | 1 | 0.0009642085 | 0.079114188 | 1.21761% |
| Two simulated ranks | on | 0 | 0.0000815620 | 0.075084875 | 0.107290% |
| Two simulated ranks | on | 1 | 0.0000778335 | 0.075714771 | 0.099016% |

Serial records **zero owned constructors/calls/checks**, so its saving is
exactly zero for this loop. The simulated rows check activation and rank-local
instrumentation only: they are **not MPI wall or scaling evidence**, and their
separate off/on batches establish neither a speedup nor timer resolution.
No comparison between serial and simulated execution models is made.
The real-MPI critical-rank wall and fraction of the full init bucket remain
for dev08; do not substitute these simulated numbers.

## Correctness tests

Fast unit gates cover:

- An independent copy of the old finite-dot loop: NaN and both infinities in
  each component of G or position; signed zero; inclusive 2^500 boundary and
  its next representable value; 1e200/1e-200 imbalanced factors; finite large
  cancellation, intermediate product overflow despite mathematical
  cancellation, and overflowing sums. Exact success/error and checked-prefix
  counts match. Empty unused arrays retain their previous behavior.
- Synthetic s/p/d, repeated l and off-diagonal coupling, plus a local-only
  species: cubic, skewed, elongated and slab-like cells, Gamma and general k,
  two radial rules, empty and sparse owned row sets. Finite values (including
  huge positions), NaN, infinity and overflowing phases test the constructor.
  **13,360 full complex coefficients** are compared by real/imaginary bits;
  selected coefficients and coupling values are also compared. Full tables
  stay absent until requested, then build exactly once per run.
- Shape-before-volume and factor-before-phase error precedence, including
  overflowing radial factors with no selected rows; profiling-off nanoseconds
  remain zero. Switch-registry and CLI switch counts become 117; the new
  switch's default is explicitly tested as off.

No tolerance is relaxed. This is constructor identity evidence, not a new SCF
energy reference, QE result, or general performance adoption claim.

## Commands and evidence

```sh
export TMPDIR="$PWD/target/scratch" TMPDIR_990="$PWD/target/scratch"
cargo test --offline --release -p dft-physics --lib phase_validation -- --nocapture
DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" PHASE_VALIDATE_BENCH=cross \
  cargo test --offline --release -p dft-physics --lib phase_validation_si64 \
  -- --ignored --nocapture --test-threads=1
CARGO_TARGET_DIR="$PWD/target" cargo build --offline --release \
  --manifest-path target/scratch/phase-harness/Cargo.toml
# Each command was run sequentially; ranks=2 uses dft_mpi::simulate, not MPI.
target/release/phase-harness off 1 ../pseudo/Si.pz-vbc.UPF
target/release/phase-harness on 1 ../pseudo/Si.pz-vbc.UPF
target/release/phase-harness off 2 ../pseudo/Si.pz-vbc.UPF
target/release/phase-harness on 2 ../pseudo/Si.pz-vbc.UPF
CHECK_SCCACHE=0 scripts/check.sh fast dft-physics dft-engine dft-numerics dft-cli
```

Raw local evidence: `target/scratch/baseline-constructor.log`,
`cross-constructor.log`, `init-{serial,sim2}-{off,on}.log`,
`measurement-summary.json`, `unit-tests.log`, `check.log`, and the ephemeral
harness source/Cargo manifest. Preserve these before target cleanup.
Crossed log SHA-256:
`95b550e428349a02778a766d13c5046700d57694051a3957e9985093e54c6c3c`.
Engine-harness binary SHA-256:
`3a7e43d60cd7c16da00e216ecedce9f58acb66e0a89cb227cab50ab27de5f52c`.
The durable isolated timing data is reproduced above; the ignored test is the
committed reproducer. Initialization harness source is preserved below.

## dev08 handoff (not executed)

Leave default off. The master should run the same isolated A/A + crossed A/B
release test first, then matched production Si64 and additional cell/species
controls at the **same MPI rank count**, with separate diagnostic and
uninstrumented acceptance batches. Start with two ranks, then the established
eight-rank workload. In a diagnostic harness, enter one RunConfig per rank,
read counters around full and coarse construction, and record each rank's
validation nanoseconds and full profile-init denominator; report critical-rank
(maximum) wall without summing ranks. Keep ordinary collective ordering.
For production arms only vary `BAANDI_PHASE_VALIDATE=off|on`, using identical
input, binary, binding, PP, grids, bands, start, threshold and verification.
Require unchanged energy/density/orbital bits, work and certificate, and apply
the measured-noise rule to whole-run time. One-core serial is an inactive
control, never an MPI denominator. No default change is authorized by the Mac
isolated result. README/dashboard artifacts are outside this ticket's owned
files; this note and the ledger carry the local evidence for master review.

## Reproducing the initialization harness

Save the following template as `target/scratch/phase-harness/template.rs`.
It deliberately does no SCF work and leaves production source untouched.

```rust
use dft_engine::{ScfProblem, ScfOptions, SpeciesPotential, parallel};
use dft_geometry::{Cell, Atom};
use dft_numerics::switches::{RunConfig, RUN_COUNTER_LABELS};
use dft_mpi::Collective;
use std::time::Instant;
use std::sync::Arc;
fn main() {
 let args: Vec<_> = std::env::args().collect();
 let mode=&args[1]; let ranks:usize=args[2].parse().unwrap();
 let pp=dft_upf::parse_upf(&std::fs::read_to_string(&args[3]).unwrap()).unwrap();
 let cell=Cell::from_user_rows_bohr(LATTICE).unwrap();
 let atoms:Vec<_>=POSITIONS.into_iter().map(|p| Atom::new("Si",p)).collect();
 let work=|comm: Option<&dft_mpi::SimRank>| {
  let rank=comm.map_or(0, |c|c.rank());
  let run=RunConfig::from_pairs(&[("phase_validate",mode),("phase_profile","on")]).unwrap();
  let _scope=run.enter();
  parallel::set_m2_override(Some(comm.is_some()));
  parallel::set_communicator(comm.map(|c|Arc::new(c.clone()) as Arc<dyn Collective+Send+Sync>));
  for sample in 0..13 {
   let before=run.counters(); let start=Instant::now();
   for cutoff in [15.,7.5] {
    let p=ScfProblem::new(&cell,&atoms,&[SpeciesPotential{label:"Si",potential:&pp}],
      ScfOptions{cutoff_hartree:cutoff,threads:1,..ScfOptions::default()}).unwrap();
    std::hint::black_box(p);
   }
   let wall=start.elapsed().as_secs_f64(); let after=run.counters();
   let stats=RUN_COUNTER_LABELS.iter().zip(after.iter().zip(before)).filter(|(k,_)|k.starts_with("phase_validate")).map(|(k,(a,b))|format!("{k}={}",a-b)).collect::<Vec<_>>().join(" ");
   println!("INIT sample={sample} mode={mode} ranks={ranks} rank={rank} seconds={wall:.9} {stats}");
  }
  parallel::set_communicator(None); parallel::set_m2_override(None);
 };
 if ranks==1 {work(None);} else { dft_mpi::simulate(ranks,|c|work(Some(c))); }
}
```

Generate its manifest and exact fixture arrays from the repository root:

```python
from pathlib import Path
import ast
root = Path.cwd()
d = root / "target/scratch/phase-harness"
(d / "src").mkdir(parents=True, exist_ok=True)
crates = ["dft-engine", "dft-physics", "dft-numerics", "dft-geometry", "dft-upf", "dft-mpi"]
(d / "Cargo.toml").write_text(
    '[package]\nname="phase-harness"\nversion="0.0.0"\nedition="2024"\n'
    '[workspace]\n[dependencies]\n'
    + ''.join(f'{c}={{path="{root}/crates/{c}"}}\n' for c in crates))
source = (d / "template.rs").read_text()
for line in (root / "benchmarks/cases/si64-gamma.toml").read_text().splitlines():
    for key, token in [("lattice", "LATTICE"), ("atoms", "POSITIONS")]:
        if line.startswith(key + " ="):
            rows = ast.literal_eval(line.split("=", 1)[1].strip())
            source = source.replace(token, repr([[float(x) for x in row] for row in rows]))
(d / "src/main.rs").write_text(source)
```

## Validation outcomes

- Release phase-validation unit suite: exit **0**, 3 passed, benchmark ignored;
  13,360 complex coefficients checked bit for bit. An initial test compile
  error used `&Vec` for a slice-target `Deref`; corrected to `&[Complex64]`
  before the successful tests (initial command exit 101).
- Isolated release A/A + crossed A/B benchmark: exit **0**, 1 passed, all
  observed counters match the expected counts and lazy-table checks pass.
- Initialization harness release build and all four executions: exit **0**.
- Required `CHECK_SCCACHE=0 scripts/check.sh fast dft-physics dft-engine
  dft-numerics dft-cli`: exit **1**, build/fmt/clippy PASS; **559 passed,
  3 failed, 152 ignored**, 56 binaries, 177 s. The three failures are the
  existing `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`,
  `pseudo::cache_resolution_and_no_network_side_effects` and
  `pseudo::cached_fetch_checks_hash_and_provenance` cache tests, rejecting
  in-repository TMPDIR. The sandbox also denied the optional `ps` load
  inventory; gates continued. No workaround changed those tests or policy.

- Existing eager-coefficient oracle with the new path enabled:
  `BAANDI_PHASE_VALIDATE=on cargo test --offline --release -p dft-engine
  --test init_owned_states
  recipe_preserves_both_radial_rules_empty_channels_and_multiple_species
  -- --nocapture`: exit **0**, 1 passed. This compares owned coefficients
  directly with the existing eager constructor for both radial rules,
  empty/sparse selections, local-only channels and multiple species.
- `git diff --check`: exit **0**.

- Isolated rerun of the existing CLI bin cache failure:
  `TMPDIR="$PWD/target/scratch" cargo test --offline -p dft-cli --bin dft-cli
  pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites
  -- --nocapture`: exit **101**, same repository-cache policy error;
  `target/scratch/known-cache-failure.log` records the complete diagnostic.
