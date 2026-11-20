# LANE-FFT Phase 2 validation

Mac correctness only. No optimization timing or speed conclusion from this host.
Base 10c067e538d1d8775629df996616399593eedac2; changes uncommitted.
`BAANDI_FFT_KEPT_FUSE` defaults off; only M2 kept-plane loading changes when on.
Serial production FFT source, reduction order and physical settings are unchanged.

## Kernel and generality

`cargo test --offline --release -p dft-solvers --lib kept_fusion_bits_counters_and_rank_independence -- --nocapture`: exit 0, one test.
Five support/grid fixtures: cubic, elongated/slab, non-symmetric skew support,
first-axis-y with ranks lacking planes/sticks, and singleton/empty support.
Ranks 1/2/3/8, 1/3/8 complex pairs (up to 16 real bands), kept counts 0/1/all,
negative and signed-zero potential entries. OFF and ON agree bitwise with
full inverse/multiply/forward; assembled outputs also agree bitwise across
rank counts. Source kept planes are unchanged. Existing slab/send/recv pointers
remain stable after batch reservation. All counter rows except kept_copy are
identical; old copy bytes = 32 * kept points; new copy counter = zero.

`BAANDI_FFT_KEPT_FUSE=on BAANDI_STAGE_TIMERS=off cargo test --offline --release -p dft-physics --test distributed_real_gamma`: exit 0, 11 tests. Existing cubic/skew/slab Hamiltonian, multiple-projector, band-width and ranks 1–8 oracles pass.

`BAANDI_FFT_KEPT_FUSE=on BAANDI_STAGE_TIMERS=off cargo test --offline --release -p dft-engine --test distributed_scf m2c_stick_row_scf_matches_threads_only -- --nocapture`: exit 0, one test. Existing strict synthetic cubic/skew SCF rank-accuracy oracle passes; this is not a bit-identity claim for full SCF across rank counts.

## Real Si UPF, simulated MPI SCF

Build:
`cargo test --offline --release -p dft-engine --test m2_entry_cost --no-run --message-format=json`: exit 0.
Select its executable from the compiler-artifact JSON, then for each ARM=off,on:

```sh
DFT_SI_UPF=/Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
M2_PROFILE_CASE=si8 M2_PROFILE_RANKS=2,3,8 M2_PROFILE_SERIAL=off \
M2_PROFILE_MODE=production M2_PROFILE_LADDER=on \
BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off \
BAANDI_FFT_PROFILE=off BAANDI_PHASE_PROFILE=off BAANDI_FFT_KEPT_FUSE="$ARM" \
"$FFT_SCF_TEST_BINARY" --ignored --nocapture --exact m2_entry_cost_profile
```

Both processes exit 0. Their existing generic wall fields are not analyzed or
stored as performance evidence. Same binary/input in both arms:

| Simulated ranks | Energy Ha (OFF = ON bits) | SCF | Davidson | H columns |
|---:|---:|---:|---:|---:|
| 2 | -31.321784923608313 | 11 | 52 | 437 |
| 3 | -31.32178492360831 | 11 | 52 | 437 |
| 8 | -31.321784923656523 | 10 | 50 | 404 |

OFF/ON energy delta 0 Ha at every rank count; density FNV hashes and work
counts agree. Cross-rank energy spread 4.82138773350016e-11 Ha, within ADR
0041's strict/standard energy classes. The existing SCF rank-dependent
trajectory is preserved, not removed. M2SCREEN exports energy bits/density
FNV/counts, not full density files, trajectory bits or verification fields;
those remain mandatory in the real-MPI dev08 A/B. Do not present FNV as SHA-256.

## Actual CLI and fixed-H harness

`cargo build --offline --release -p dft-cli --bin dft-cli --example davidson_bench`: exit 0.
Using the retained Phase-0 Si8 input (`target/scratch/fft-si8.toml`), run each
ARM=off,on, profile and stage timers off:

```sh
BAANDI_FFT_KEPT_FUSE="$ARM" BAANDI_FFT_PROFILE=off BAANDI_PHASE_PROFILE=off \
BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off \
target/release/dft-cli run target/scratch/fft-si8.toml \
  --output-dir "target/scratch/fft-phase2-cli-$ARM"
```

Each exit 0. Serial control energy -31.321776563202498 Ha, verified in both;
energy/eigenvalues/occupations/scientific histories and density bytes identical.
This protects the unchanged serial default; it does not exercise the new branch.

The real-UPF fixed-H warm solve **does** exercise that branch on simulated MPI3:

```sh
BAANDI_FFT_KEPT_FUSE="$ARM" BAANDI_FFT_PROFILE=on BAANDI_PHASE_PROFILE=on \
BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off \
target/release/examples/davidson_bench --input target/scratch/fft-si8.toml \
  --kinds warm --repeats 1 --sim-ranks 3
```

Each exit 0. Identical eigenvalue FNV `0xfdff2330905ef17e`, eigenvalues,
residuals and work counters. Actual `phase_profile.fft_lane` JSON shows, on
every rank in fine: kept_copy 288 -> 0 calls, 11943936 -> 0 logical bytes.
All other stage counters (including multiply and collective payloads) agree.
Total traffic decreases by exactly 11943936 B/rank in this harness profile
(includes setup/re-arm/repeats; not just the measured solve). All lane clocks
are zero. Deterministic work reduction is evidence; no Mac speed claim.

[phase2-equivalence.json](phase2-equivalence.json) retains numerical evidence,
source/binary/UPF hashes. Recheck with:

```sh
python3 -I docs/performance/lanes/FFT/phase2_verify.py target/scratch \
  --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF
```

Exit 0. Raw local logs: `target/scratch/fft-phase2-{kernel,physics,synthetic-scf}.log`,
`fft-phase2-scf-{off,on}.log`, `fft-phase2-dav-{off,on}.log`, CLI output dirs.

## Required gate

`CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-numerics dft-solvers dft-cli`: **exit 1**. Build/fmt/clippy PASS. 358 tests pass,
87 ignored, three known Mac in-repository TMPDIR failures:

- pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites
- pseudo::cache_resolution_and_no_network_side_effects
- pseudo::cached_fetch_checks_hash_and_provenance

No failed test is reported as passing. No full/push gate or native MPI/MKL run
is claimed. Source/build/test logs are retained in target/scratch. Shell syntax
checks of the request/runner pass. Defaults remain OFF; nothing was committed
or run remotely. Stage-3 broader materials/cell/PP qualification and real-MPI
verification belong to adoption after a resolved dev08 speed result.
