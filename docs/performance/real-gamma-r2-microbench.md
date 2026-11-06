# P26B / R2: paired full-sphere FFT library validation

2026-11-06. **Mac observation**, not a dev08 acceptance result or QE comparison.
Source base: `f176bc7d78e3cc989671ecf49cfcbb1d4712d396`, branch
`task/P26B-REAL-GAMMA-R2`, plus the uncommitted R2 changes. Host: macOS 27.0.1
(26A434), arm64; rustc 1.99.1 (58cf8f9c3 2026-10-13), ordinary release profile,
no native CPU flags. CPU model lookup was denied by the sandbox. No remote
machine or pseudopotential was changed. No QE source was read.

## Implementation and contracts

Derived from `real-gamma-design.md` sections 3, 4, 6 and 11, itself based on
architect-written Gamma/FFT study notes. This implements the new pair-transform
mechanism, not rejected band grouping A28 (which reused potential traversal
without reducing transform count). No rejected lane staging, threaded FFT,
scoped worker creation, nonlocal traversal, or solver tuning was repeated.

- `PairMap` accepts the arrays returned by `HalfSphere::plus/minus`, together
  with the ordered full-sphere FFT indices. It validates G=0, full coverage,
  unique indices, and inversion on the grid. Passing arrays avoids a dependency
  cycle: dft-basis already depends on dft-numerics.
- `PairPlan<P>` owns complex input/output scratch. ORC conversion uses sqrt(2),
  packs a+i*b, invokes the existing local application exactly once, then separates
  by parity. G=0 keeps only the real coordinate for each band. Odd final columns
  use the same routine with an implicit zero partner. Nonzero norms must be within
  a factor 1e3 in debug builds; an exactly zero partner is exempt.
- The additive `PairDensityOperator` implementations for `OrderedFusedPlan` and
  `GridFusedPlan<F>` reuse their existing inverse stages, without a forward
  transform. They add `(w_a*Re*Re + w_b*Im*Im)/volume` in e/Bohr^3. Only two
  exactly zero weights skip the transform, including signed zero. An isolated
  zero weight and a tiny nonzero weight do not skip the pair.
- `PairDensityWorkspace<P>` has four pairs per fixed piece, one plan and one
  density grid per piece. Each piece sums pairs in column order, followed by
  ascending-piece-order reduction. Recursive joins use the existing explicit
  WorkerPool, including nested calls. The grouping never depends on threads.
  Memory is O(ceil(bands/8) * (plan scratch + grid)); this simple library policy
  is not a tuned performance default or an engine integration.
- All new numerical storage is reserved at construction. Hot calls have no
  allocation sites, logging, I/O or timers, and no new unsafe code. A storage
  test checks pointers, lengths and capacities across repeated calls. Existing
  FFT allocation-free contracts are inherited; no global allocator instrumentation
  was added and external-backend hidden allocations were not measured.

The architect explicitly extended ownership to `fused.rs` for additive changes.
Only two new trait implementations were inserted there. Removing that insertion
reproduces the original file byte-for-byte, including all existing method bodies,
signatures, visibility and tests. `lib.rs` only adds a module declaration and
re-exports. No engine caller or existing numerical kernel was changed.

## Component evidence

`cargo test -p dft-numerics --test real_gamma_pair -- --nocapture`: exit 0,
5 passed and 1 ignored, 0.53 s in debug. The reusable-storage unit test also
passes in the crate/workspace gates. Both ordered and generic Rust fallback
plans produce the following relative L2 errors against two single complex
applications/densities:

| Synthetic reciprocal metric | Grid | Full coefficients | Local relative error | Density relative error |
| --- | --- | ---: | ---: | ---: |
| Cubic | 11x11x11 | 123 | 2.4514960468393703e-16 | 2.179748991757297e-16 |
| Skewed | 11x12x13 | 107 | 2.2341913080719364e-16 | 1.6554504145690005e-16 |
| Elongated | 11x12x19 | 107 | 2.3266948064601855e-16 | 1.767951968498147e-16 |
| Slab | 13x14x21 | 283 | 2.4372380784229483e-16 | 1.9449585680482307e-16 |

These are synthetic inversion-symmetric ellipsoidal supports, not material
calculations. Each uses 19 normalized ORC bands, fractional and zero weights,
a nonconstant real potential in Ha, and volume 120 Bohr^3. Relative errors are
dimensionless. The fixed gate is **1e-14**, allowing rounding of the combined
function in the FFT, sqrt(2) conversion, and the different density grouping.
The same gate is also checked against completed full complex output, including
imaginary coefficients. No tolerance was relaxed.

Density output bits agree across budgets **1, 2, 4, 8**, on repeated and nested
calls, for all four shapes and both backends. Other tests cover empty blocks,
G=0-only and odd blocks, G=0 imaginary leftovers, mirror/dimension/volume errors,
and the debug norm assertion. A counted generic FFT verifies two inverse and
two forward calls for three bands, and inverse-only density with exact-zero
skip behavior; invalid weights leave density and transform counts unchanged.

## Pair application microbenchmark (Mac observation)

```sh
cargo test --release -p dft-numerics --test real_gamma_pair \
  pair_microbenchmark -- --ignored --nocapture
```

Exit 0. The ignored test uses the unchanged ordered fused plan, normalized
synthetic real bands, a nonconstant potential, and spherical support of radius
one quarter of the grid dimension. Pair timing includes ORC packing and parity
unpacking. Single timing uses already completed complex coefficients and two
separate applications. Plan construction and allocation are outside timing.
Each sample has 20 pair applications or 40 single applications, divided by 40
for microseconds per band. Two alternating warmup rounds are discarded; 30
alternating measurement rounds are retained. Values are medians (average of
middle two). Our verification processes had completed before the final run;
the host was not reserved exclusively.

| Grid | Full coefficients | Pair us/band | Single us/band | Pair / singles per-band ratio |
| --- | ---: | ---: | ---: | ---: |
| 25^3 | 1021 | 64.850000 | 128.542712 | 0.504502 |
| 36^3 | 3071 | 176.410937 | 349.802088 | 0.504316 |
| 48^3 | 7153 | 419.278650 | 828.820837 | 0.505874 |

No gain gate is claimed. This measures only local application on synthetic
supports, not density throughput, SCF speed, MKL runtime, or general production
performance. dev08 decides performance acceptance. An earlier exploratory run
used 251 coefficients at all three grids while other checks ran; it was not
used for the final table.

## Complex-engine isolation

A release binary was built before adding the module and additive methods, then
another from the changed tree. Both were run with the same UPF and five cases.
The first command with a relative UPF path exited 3 because the harness resolves
it inside its temporary input directory; it was corrected to the absolute path.
No baseline numbers or tolerances were changed to obtain agreement.

```sh
cargo build --release -p dft-cli
# Before edits: copy target/release/dft-cli to /tmp/p26b-before-cli.
python3 -I scripts/bench_scf.py --binary /tmp/p26b-before-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 1 \
  --record /tmp/p26b-before.json
# After edits: rebuild, then compare against that baseline.
python3 -I scripts/bench_scf.py \
  --binary /Users/mino/Develop/rust-dft-wt/P26B-REAL-GAMMA-R2/target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 1 \
  --baseline /tmp/p26b-before.json --require-identical \
  --record /tmp/p26b-after.json
```

Both builds and both absolute-path harness runs exited 0. Every density SHA-256
matched and every energy difference was **0 Ha**:

| Case | SCF iterations | Total energy (Ha), before and after |
| --- | ---: | ---: |
| gamma | 10 | -7.295530468210760 |
| k222s | 10 | -7.9194234124690865 |
| k222 | 10 | -7.831230057121135 |
| k444s | 10 | -7.9259761266480915 |
| si8 | 15 | -31.32180360103127 |

SHA-256 provenance:

- UPF: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
- Before binary: `21ad13f204a44710849269946b05506d57cd00db9a576e535cbe65e24fbf3c9f`.
- After binary: `e58721a9d58cbdceba766f5a0ff275ff4d3ff9de5eb29ab1cbb9026f18d62288`.
- gamma density: `6077d62f5ddcf7d003e61b29f12e83d73e0dc332a208097da0836dfb63ec1b26`.
- k222s density: `9e311684370831efd1f4e5a3f3fd41479dcadfc4794b51e2db0183f05e7e9104`.
- k222 density: `21238c7026d96860233edfc427807389029229e7ce26162b5025d471e255bdbe`.
- k444s density: `767602ec4395d37bf7a9da12cf4a2a399c3fee89d7cdf4129a05b39dee127660`.
- si8 density: `826e1046923c4674f619f1ce2dbac27906bb9c6b66ae09b02aa1687d7680766e`.

Timing from these single SCF runs is not performance evidence. This is BAANDI
before/after isolation, not a new QE reference. Raw logs/JSON remain in `/tmp`
with prefix `p26b-`; the results above are the durable record within ownership.

## Gates and handoff

| Command | Exit | Outcome |
| --- | ---: | --- |
| `scripts/check.sh dft-numerics` | 0 | build, test, fmt, clippy pass |
| `scripts/check.sh` | 0 | workspace build, test, fmt, clippy pass |
| `cargo check -p dft-numerics --features mkl` | 0 | optional backend type-checks on this Mac; no MKL runtime/link claim |
| component test and release microbenchmark commands above | 0 | pass |
| `bench_scf.py ... --require-identical` above | 0 | all five cases identical |
| `python3 -I scripts/make_readme_figures.py` | 0 | generated artifacts unchanged |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 9 tests pass; existing reader ResourceWarnings only |

No new material/README evidence type was introduced: this ticket owns only the
R2 note and crate files, not the README evidence generator or datasets. The
existing generated validation block remains fresh. The orchestrator owns ledger,
backlog, dashboard/session registry updates and any later README integration.
Suggested ledger entry R2: library-only ORC pair application/density; local errors
<=2.46e-16 and density <=2.18e-16; T=1/2/4/8 identical; Mac pair/singles ratios
0.504502/0.504316/0.505874; five unchanged SCF cases; dev08 pending. No optimization
is adopted into the engine, and no performance default changes. Changes are left
uncommitted for architect review.
