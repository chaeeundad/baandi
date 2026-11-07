# P39B: caller-owned MKL in-place FFT

Baseline/current main: `010b4e73e54880a7b9d2872158881c729e62176b`.
Branch: `task/P39B-MKL-ZERO-COPY`; changes uncommitted. Local host: Apple
Silicon arm64, rustc 1.99.1 (58cf8f9c3 2026-10-13), no MKL. Remote runtime validation and speed acceptance remain pending.

## Mechanism and scope

The motivation is the **relayed** 2026-11-08 raw FFT study in
`tests/reference/perf/dev08_materials_qe_relayed.md`: wrapper/raw inverse medians
191.13/147.25 us at 36^3 and 486.30/323.29 us at 72x36x36; 25^3 showed no gain.
Those are prior measurements, not measurements of this patch. This follows B33
and PFFT-6; it does not revive A34 packed-lane threads or A35 scoped GEMM threads.
No threading, algorithm, default, normalization or operation order changes.
No QE source was consulted; no FFTW code or dependency was added.

`Plan3d::transform_in_place` validates the exact caller slice length before
calling the already committed in-place descriptor. Both directions retain its
existing scale (backward unscaled, forward 1/N). Exclusive mutable borrows cover
the descriptor and caller storage; no pointer is retained. DFTI errors propagate,
and on an error the caller's buffer contents are unspecified. The old `transform`
body is unchanged and remains the runtime bit-identity oracle.

The additive `Fft3d::transform_in_place(data, scratch, inverse)` accepts a second,
preallocated array so its default can call the existing out-of-place operation
and copy its result back without allocating or changing arithmetic. Both arrays
must have the plan's shape and C layout, including for the MKL override, which
leaves scratch unused and invokes the new primitive. The default adds a copy
for a custom/Rust `GridFusedPlan`; the production pure-Rust ordered and plane
plans and the pure-Rust pair plan are unchanged.

`GridFusedPlan` scatters into its real array, transforms that same storage,
multiplies by the potential in the same iteration order, transforms in place
again, and gathers from it. Single and pair density reuse the same inverse.
The reciprocal array remains fallback scratch. On MKL, each FFT eliminates two
N-element complex copies (32N bytes of copied payload, or 64N bytes counting
reads plus writes). At 25^3/36^3/72x36x36 this removes respectively
500000/1492992/2985984 bytes of copied payload per transform. Full-grid zero-fill
and scattering remain. Existing out-of-place callers with immutable input,
including the public pruned APIs, preserve their API and copying behavior;
there is no other eligible owned full-grid scratch transform in numerics.
No P39A-owned file is changed. No new unsafe code exists outside dft-mkl.

Only named owned files are edited. The dft-mkl lib addition is a private length
guard shared by the FFI entry and a no-MKL unit-test module under tests/contracts.
This tests the actual pre-FFI length predicate without pretending to execute MKL.
All numerical storage is reserved before calls. The old Plan3d buffer is retained
for compatibility; this change does not reduce plan allocation size.

## Tests and local validation

The no-MKL tests exercise exact/short/long/zero lengths, wrong shapes and
Fortran layout, default-method bit identity, repeated fused local application,
single density, pair density and stable array addresses. MKL-gated tests compare
both transform directions on deterministic xorshift random and sparse inputs,
including offset caller slices, at 3x4x5, 25^3, 36^3, 72x36x36 and 7x9x11.
The numerics MKL test separately checks a scattered pruned pattern through local
application and both density APIs against the old two-array wrapper.

MKL tests are compiled out without the feature, and automatically run with it.
The microbenchmark alone is ignored, requiring explicit release execution. It
prints all 32 samples and median/min/max us per call for each direction and
shape. Eight warmups precede 32 alternating wrapper/in-place pairs. Inputs are
reset outside the timed region to avoid overflow from repeated inverse calls;
allocation and plan creation are also outside timing. This is a warm-input
single-call benchmark; it differs from the earlier relay's repeated-transform
inner loop and must not be presented as an identical measurement protocol.

Local command exit codes and SCF evidence are recorded below after completion.
No artifacts are edited, per ticket instruction. README generation and freshness
checks run, but there is no new external performance claim to add to its data.

## dev08 commands (orchestrator only; not executed on the Mac)

Use the reviewed patch and retain the baseline executable from the exact parent
above before rebuilding into the same target tree. Record source/binary hashes,
rustc version, CPU, load, MKL version and full output. Do not run benchmark builds
concurrently. The user supplies MKL and the external Si UPF.

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
source "$MKLROOT/env/vars.sh"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
cargo test --offline --release -p dft-mkl --features mkl --test in_place in_place_matches_wrapper_bits -- --nocapture
cargo test --offline --release -p dft-numerics --features mkl --test in_place -- --nocapture
cargo test --offline --release -p dft-numerics --features mkl
# Pin the benchmark harness and its single test thread to one CPU.
taskset -c 0 cargo test --offline --release -p dft-mkl --features mkl --test in_place in_place_microbenchmark -- --ignored --nocapture --test-threads=1
```

Run the SCF command sequence below first on baseline MKL, then on candidate MKL,
using the same UPF, build flags and thread budget. `P39B_BASE_MKL` must name the
preserved baseline executable; `P39B_UPF` the external Si.pz-vbc.UPF. Repeat with
`--extra-toml 'gamma_real = "on"'` and distinct records for gamma and si8 to check
the MKL pair consumer as well as the default complex path. Compare energies by
`float.hex()` as well as density hashes; bench_scf's energy tolerance alone is
not a bit-identity assertion. Stop on any bit difference; do not relax the gate.

```sh
cargo build --offline --release -p dft-cli --features mkl
mkdir -p /tmp/p39b-mkl
for t in 1 2 4 8 12; do
  python3 -I scripts/bench_scf.py --binary "$P39B_BASE_MKL" --upf "$P39B_UPF" --cases gamma,k222s,k222,k444s,si8 --threads "$t" --repeats 1 --record "/tmp/p39b-mkl/before-$t.json" || exit $?
  python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf "$P39B_UPF" --cases gamma,k222s,k222,k444s,si8 --threads "$t" --repeats 1 --baseline "/tmp/p39b-mkl/before-$t.json" --require-identical --record "/tmp/p39b-mkl/after-$t.json" || exit $?
done
python3 - <<'PY'
import json
for t in [1,2,4,8,12]:
    before = json.load(open(f'/tmp/p39b-mkl/before-{t}.json'))['cases']
    after = json.load(open(f'/tmp/p39b-mkl/after-{t}.json'))['cases']
    for case, a in before.items():
        b = after[case]
        assert a['total_energy_ha'].hex() == b['total_energy_ha'].hex(), (t,case)
        for key in ['density_sha256','iterations','operator_applications','davidson_iterations','certificate_applications','certificate_iterations']:
            assert a[key] == b[key], (t,case,key)
print('25/25 MKL SCF comparisons bit-identical with unchanged work')
PY
```

For speed acceptance, rerun paired SCFs at least three times per condition on an
idle dev08, fixed CPU affinity, and report spread/work counts. The 32-sample FFT
timing establishes adapter cost only, not end-to-end gains or QE parity.

## Completed local evidence (2026-11-08)

All 25 comparisons pass: energy bits equal (maximum delta **0 Ha**), density
SHA-256 equal, and SCF/Davidson/operator/certificate work counts unchanged.
Gamma/k222s/k222/k444s take 10 SCF iterations; si8 takes 15, at every budget.
The rows below record values shared by baseline and candidate. This verifies
pure-Rust correctness, not performance: builds/tests overlapped some SCF runs,
and these single samples are deliberately not accepted as timing evidence.

Baseline release binary SHA-256:
`5559a33a2ac4b378c4d7fa68e5086ea37e4546df7c354f7aa6c4b42bca1d17d3`.
Candidate release binary SHA-256:
`ceecd6b37f3d8ff098ad7ce78ee850ae9ccc04aec848e5b77a550294829863cd`.
External Si.pz-vbc.UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Only the existing target directory was used. The unchanged-source baseline
was built first, and its executable copied to `/tmp/p39b-evidence/baseline-cli`
before source edits. Raw local command logs and before/after JSON records are
in `/tmp/p39b-evidence/` (ephemeral); the durable bit-identity table follows.

| T | Case | Energy (Ha) | Energy bits (hex float) | SCF iterations | Density SHA-256 |
|---|---|---|---|---|---|
| 1 | gamma | -7.2955304682107602 | -0x1.d2e9f89ffbf44p+2 | 10 | 6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26 |
| 1 | k222s | -7.9194234124690865 | -0x1.fad7d54beed00p+2 | 10 | 9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104 |
| 1 | k222 | -7.8312300571211351 | -0x1.f532df8db2663p+2 | 10 | 21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe |
| 1 | k444s | -7.9259761266480915 | -0x1.fb43315f3523cp+2 | 10 | 767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660 |
| 1 | si8 | -31.32180360103127 | -0x1.f5261b8862a12p+4 | 15 | 826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e |
| 2 | gamma | -7.2955304682107602 | -0x1.d2e9f89ffbf44p+2 | 10 | 6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26 |
| 2 | k222s | -7.9194234124690865 | -0x1.fad7d54beed00p+2 | 10 | 9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104 |
| 2 | k222 | -7.8312300571211351 | -0x1.f532df8db2663p+2 | 10 | 21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe |
| 2 | k444s | -7.9259761266480915 | -0x1.fb43315f3523cp+2 | 10 | 767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660 |
| 2 | si8 | -31.32180360103127 | -0x1.f5261b8862a12p+4 | 15 | 826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e |
| 4 | gamma | -7.2955304682107602 | -0x1.d2e9f89ffbf44p+2 | 10 | 6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26 |
| 4 | k222s | -7.9194234124690865 | -0x1.fad7d54beed00p+2 | 10 | 9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104 |
| 4 | k222 | -7.8312300571211351 | -0x1.f532df8db2663p+2 | 10 | 21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe |
| 4 | k444s | -7.9259761266480915 | -0x1.fb43315f3523cp+2 | 10 | 767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660 |
| 4 | si8 | -31.32180360103127 | -0x1.f5261b8862a12p+4 | 15 | 826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e |
| 8 | gamma | -7.2955304682107602 | -0x1.d2e9f89ffbf44p+2 | 10 | 6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26 |
| 8 | k222s | -7.9194234124690865 | -0x1.fad7d54beed00p+2 | 10 | 9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104 |
| 8 | k222 | -7.8312300571211351 | -0x1.f532df8db2663p+2 | 10 | 21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe |
| 8 | k444s | -7.9259761266480915 | -0x1.fb43315f3523cp+2 | 10 | 767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660 |
| 8 | si8 | -31.32180360103127 | -0x1.f5261b8862a12p+4 | 15 | 826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e |
| 12 | gamma | -7.2955304682107602 | -0x1.d2e9f89ffbf44p+2 | 10 | 6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26 |
| 12 | k222s | -7.9194234124690865 | -0x1.fad7d54beed00p+2 | 10 | 9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104 |
| 12 | k222 | -7.8312300571211351 | -0x1.f532df8db2663p+2 | 10 | 21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe |
| 12 | k444s | -7.9259761266480915 | -0x1.fb43315f3523cp+2 | 10 | 767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660 |
| 12 | si8 | -31.32180360103127 | -0x1.f5261b8862a12p+4 | 15 | 826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e |

Commands actually executed (PATH included `/opt/homebrew/opt/rustup/bin`):

| Command | Exit |
|---|---|
| `cargo build --offline --release -p dft-cli` before edits | 0 |
| `cargo build --offline --release -p dft-cli` after edits | 0 |
| `scripts/check.sh` (no argument means all workspace crates; build/test/fmt/clippy PASS) | 0 |
| `cargo check --offline -p dft-numerics --features mkl --all-targets` | 0 |
| `cargo clippy --offline -p dft-numerics -p dft-mkl --features dft-numerics/mkl --all-targets -- -D warnings` | 0 |
| `python3 -I scripts/make_readme_figures.py` (no generated changes) | 0 |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 |
| `git diff --check` | 0 |
| `cargo test --offline -p dft-mkl caller_buffer_length_contract` (1 contract test passed) | 0 |
| `cargo test --offline -p dft-numerics --test in_place` (2 tests passed) | 0 |
| `bench_scf.py`, baseline loop below, all five T values | 0 |
| `bench_scf.py`, candidate loop below, all five T values | 0 |
| Python post-check of energy `.hex()`, density SHA-256 and all five recorded work fields | 0 |

The MKL feature checks explicitly warn that MKLROOT is unset and executables
cannot link. No native link or MKL runtime pass is claimed. Both MKL crates'
new test targets are covered by the feature-enabled all-target clippy check.

Exact local SCF commands (two loops run before/after the candidate build):

```sh
for t in 1 2 4 8 12; do
  python3 -I scripts/bench_scf.py --binary /tmp/p39b-evidence/baseline-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,k444s,si8 --threads "$t" --repeats 1 --record "/tmp/p39b-evidence/before-$t.json" || exit $?
done
for t in 1 2 4 8 12; do
  python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,k444s,si8 --threads "$t" --repeats 1 --baseline "/tmp/p39b-evidence/before-$t.json" --require-identical --record "/tmp/p39b-evidence/after-$t.json" || exit $?
done
```

Scope limits: no new QE comparison, no measured MKL speedup, no proof yet that
MKL chooses bit-identical arithmetic for every caller alignment. The remote
bit tests explicitly exercise this risk; exact copies alone do not prove that
a vendor library selects the same code path for differently aligned storage.
Runtime MKL and its SCF identity are the remaining orchestrator gates.
