# Backend spike results

Measured on 2026-11-02. These are dimensionless synthetic experiments, not physical energies or QE references. The reproducible test is `crates/dft-spike-backend/src/experiments.rs`; the complete final release output is [release.txt](../crates/dft-spike-backend/measurements/release.txt).

## Environment and procedure

Host probes returned macOS 27.0.1 (26A434), arm64, and rustc 1.99.1 (58cf8f9c3 2026-10-13). CPU model and logical CPU count probes were denied by the sandbox. `otool -L` confirmed that the release test executable links `/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate` (current version 4.0.0).

Exact versions read from Cargo.lock:

| Crate | Version |
| --- | --- |
| faer | 0.24.4 |
| lapack | 0.20.0 |
| blas-src | 0.14.0 |
| lapack-src | 0.13.0 |
| accelerate-src (transitive) | 0.3.2 |
| ndarray | 0.17.2 |
| num-complex | 0.4.6 |
| ndrustfft | 0.6.2 |
| rustfft | 6.4.1 |

No dependencies were added or upgraded. The source-selector Accelerate features still depend on `accelerate-src`; they do not eliminate that older transitive crate.

All reported timings are elapsed seconds, median of three calls after one warm-up, in a release test run with no other agent-launched build/test running. Inputs are fixed, and their construction is outside timing. Output allocation and adapter copies are inside timing. Eigensolve includes all eigenvectors, LAPACK workspace query/allocation and column-major conversion; quality checks are outside timing. GEMM includes result allocation. FFT handler construction is outside timing, while input clone and output allocation are inside. Timings describe these straightforward adapters, not isolated kernels. Three samples are not a statistical performance study.

Faer global parallelism and explicit matmul parallelism are `Seq` (one calling thread). ndrustfft uses serial `ndfft`/`ndifft`, not the parallel variants. `VECLIB_MAXIMUM_THREADS=1` is set before process startup; actual Accelerate worker count was not observable. No claim about default multithread performance is made.

The deterministic generator is SplitMix64, initialized once with seed `0x5eed`, taking the upper 53 bits and mapping to [-0.5, 0.5). Real and imaginary components are successive draws. Each random Hermitian matrix uses the lower triangle and conjugates it into the upper triangle; its diagonal is real. The generator stream continues across matrix and FFT cases. No rand API is used by the spike (faer has transitive random dependencies).

## Tolerances and definitions

All comparisons assert strict bounds without loosening the requested tolerances. Eigenvalues are in ascending order. For each eigensystem, `s=max |lambda|`; residual means `max_j ||A v_j-lambda_j v_j||_2`, and relative residual is residual/s. The eigenvalue cross-comparison uses the LAPACK spectral scale; known-value errors use each solver's scale. Orthonormality is `max_ij |(V^H V-I)_ij|`. Residual and Gram products use faer matmul, which is separately checked against the independent triple loop.

| Check | Tolerance (dimensionless) |
| --- | --- |
| Eigenvalue error / spectral scale | 1e-9 |
| Residual / spectral scale | 1e-9 |
| Orthonormality, absolute | 1e-10 |
| GEMM max entrywise error / max reference magnitude | 1e-12 |
| FFT round trip max entrywise error / max input magnitude | 1e-12 |
| Parseval: abs(sum(abs(F)^2)/M - sum(abs(x)^2))/sum(abs(x)^2) | 1e-12 |
| Naive DFT relative maximum entrywise error | 1e-12 |
| Plane-wave peak error/M and maximum off-peak magnitude/M | 1e-12 |
| Inverse unit coefficient gain: abs(M*gain-1) | 1e-12 |
| Handler reuse relative maximum entrywise difference | 1e-12 |

## Complex Hermitian eigensolve

The analytic case is [[2, 1-i], [1+i, 3]], with eigenvalues [1, 4]. The second case constructs `A=U D U^H`, with `U[j,k]=exp(2*pi*i*j*k/4)/2` and `D=[-3, 0.5, 2, 7]`. This complex Fourier unitary mixes every coordinate.

| Case/N | Max eigenvalue difference | Relative difference | Faer (s) | LAPACK zheev (s) |
| --- | --- | --- | --- | --- |
| analytic2/2 | 1.3322676295501878e-15 | 3.3306690738754706e-16 | 0.000000333 | 0.000000375 |
| unitary4/4 | 8.8817841970012523e-16 | 1.2688263138573217e-16 | 0.000001000 | 0.000001292 |
| seeded/100 | 4.2632564145606011e-14 | 5.4514498615573778e-15 | 0.000745250 | 0.001377333 |
| seeded/300 | 6.3948846218409017e-14 | 4.6428603043506600e-15 | 0.013085750 | 0.024361917 |
| seeded/600 | 1.1013412404281553e-13 | 5.5510471282892966e-15 | 0.085498250 | 0.152928917 |
| seeded/1000 | 2.1316282072803006e-13 | 8.2813047938168331e-15 | 0.368617583 | 0.628785708 |

| Case/N | Backend | Max residual | Relative residual | Max orthonormality error |
| --- | --- | --- | --- | --- |
| analytic2/2 | faer | 4.4408920985006262e-16 | 1.1102230246251565e-16 | 1.1146791852811679e-16 |
| analytic2/2 | lapack | 1.9389211565826713e-15 | 4.8473028914566801e-16 | 3.3315242806850840e-16 |
| unitary4/4 | faer | 1.6727693978859666e-15 | 2.3896705684085238e-16 | 5.5511151231257827e-16 |
| unitary4/4 | lapack | 2.0654203592503655e-15 | 2.9506005132148080e-16 | 4.4408920985006262e-16 |
| seeded/100 | faer | 5.4779127622368303e-14 | 7.0046377382621203e-15 | 5.5511151378530695e-15 |
| seeded/100 | lapack | 1.8962609629207558e-14 | 2.4247595168062277e-15 | 4.4408930823968676e-15 |
| seeded/300 | faer | 6.1122859398350465e-14 | 4.4376859688723364e-15 | 4.4408926996308334e-15 |
| seeded/300 | lapack | 6.7467422437969023e-14 | 4.8983185154626004e-15 | 5.3290705257733238e-15 |
| seeded/600 | faer | 8.0847293971347296e-14 | 4.0749145002064825e-15 | 3.5527136972672117e-15 |
| seeded/600 | lapack | 1.4156354450720492e-13 | 7.1351719009603040e-15 | 6.3282712455138778e-15 |
| seeded/1000 | faer | 1.5049548707125865e-13 | 5.8466996931005454e-15 | 6.2172490607966572e-15 |
| seeded/1000 | lapack | 3.2662625393940051e-13 | 1.2689321492822173e-14 | 1.0880185659696539e-14 |

| Case/N | Backend | Measured eigenvalues | Max known-value error |
| --- | --- | --- | --- |
| analytic2/2 | faer | [0.9999999999999998, 4.0] | 2.2204460492503131e-16 |
| analytic2/2 | lapack | [0.9999999999999993, 3.9999999999999987] | 1.3322676295501878e-15 |
| unitary4/4 | faer | [-3.0, 0.5000000000000004, 2.0000000000000004, 7.0] | 4.4408920985006262e-16 |
| unitary4/4 | lapack | [-3.000000000000001, 0.5000000000000002, 2.0000000000000004, 7.0] | 8.8817841970012523e-16 |

## Complex matrix multiply

Both ordinary `A B` and conjugate-transposed `A^H B` use independent dense complex inputs. The reference is a literal triple loop with complex products and sums. Both sizes are checked, including the timing case.

| N | A adjoint | Relative error | Faer (s) | Naive (s) |
| --- | --- | --- | --- | --- |
| 64 | false | 8.1340060812670084e-16 | 0.000042208 | 0.000195000 |
| 64 | true | 1.0607729545390141e-15 | 0.000042333 | 0.000195709 |
| 500 | false | 2.4771621510247132e-15 | 0.018948250 | 0.092120250 |
| 500 | true | 3.7153342879757139e-15 | 0.018918500 | 0.093906833 |

## Three-dimensional FFT

Arrays have standard ndarray C layout, with z contiguous and shape (nx, ny, nz). A transform calls the library along axes x, y, then z, swapping two array buffers. It does not implement its own per-lane rustfft loop. The reference on 4x5x6 explicitly sums over all spatial points for every output coefficient (quadratic in the total grid size).

For every grid, the input plane wave is `exp(2*pi*i*(x/nx-y/ny+2*z/nz))`. Its expected peak is at wrapped index `(1, ny-1, 2)`. The measured peak equals the grid size M, with only roundoff leakage. An inverse transform of a single unit coefficient yields 1/M at the origin. Combined with the measured identity round trip, these establish forward-unscaled and inverse-1/M normalization, rather than assuming rustfft conventions. The inverse applies scaling along each axis. No extra application scaling was added.

| Grid | Round trip relative | Parseval relative | Forward (s) | Inverse (s) | Reuse relative difference |
| --- | --- | --- | --- | --- | --- |
| 4x5x6 | 4.2051664774730458e-16 | 4.1319633560683846e-16 | 0.000001667 | 0.000001708 | 0.0000000000000000e0 |
| 15x20x24 | 7.9113646052180291e-16 | 1.9230058421511334e-15 | 0.000048500 | 0.000056042 | 0.0000000000000000e0 |
| 16x16x16 | 5.2875382000344749e-16 | 4.1218347524886851e-15 | 0.000021833 | 0.000025000 | 0.0000000000000000e0 |
| 27x27x27 | 1.6254867190301116e-15 | 6.9433982789568414e-15 | 0.000203333 | 0.000207084 | 0.0000000000000000e0 |
| 32x32x32 | 8.6496122222378266e-16 | 1.1648193178396638e-15 | 0.000231833 | 0.000251083 | 0.0000000000000000e0 |
| 48x48x48 | 1.3717616851043414e-15 | 1.9162370907909160e-14 | 0.001215375 | 0.001188917 | 0.0000000000000000e0 |

| Grid | Wrapped peak index | Peak real | Peak imaginary | Max leakage / M | Inverse unit gain |
| --- | --- | --- | --- | --- | --- |
| 4x5x6 | (1, 4, 2) | 1.2000000000000000e2 | -3.9497394250679243e-14 | 2.0899301900814605e-16 | 8.3333333333333332e-3 |
| 15x20x24 | (1, 19, 2) | 7.2000000000000000e3 | -9.3332699556560864e-13 | 2.3920403074011821e-16 | 1.3888888888888889e-4 |
| 16x16x16 | (1, 15, 2) | 4.0960000000000000e3 | -6.7252069806088614e-13 | 1.5793228746967149e-16 | 2.4414062500000000e-4 |
| 27x27x27 | (1, 26, 2) | 1.9682999999999996e4 | -4.2391437941537538e-12 | 2.2381112914806013e-16 | 5.0805263425290857e-5 |
| 32x32x32 | (1, 31, 2) | 3.2768000000000000e4 | -4.1087324599124956e-12 | 1.3326293213617709e-16 | 3.0517578125000000e-5 |
| 48x48x48 | (1, 47, 2) | 1.1059200000000000e5 | 9.9083249731378124e-12 | 2.9075522147056382e-16 | 9.0422453703703699e-6 |

| Grid | Naive DFT relative error | Naive DFT (s) |
| --- | --- | --- |
| 4x5x6 | 1.7631359185589452e-15 | 0.000062291 |

The same three handlers survive warm-up, all timed calls, round trip, plane-wave tests, inverse unit probe, and a final repeated forward transform. Reuse gives exactly zero difference on these inputs. This verifies plan reuse, not allocation-free execution. Vendored ndrustfft source handles strided axes with temporary storage; scratch/allocation counts were not instrumented.

## API and build observations

Faer owned matrices support complex values, adjoint views, a self-adjoint decomposition returning sorted diagonal values and eigenvector columns, and matmul with explicit parallelism. The adapter is safe Rust and uses column-oriented indexing. LAPACK requires an explicit contiguous column-major copy, 32-bit dimensions, triangle selection, workspace query, and INFO checking. All unsafe calls are confined to `experiments/native.rs`; initialized buffer lengths and query results are checked. SDK clapack.h declares the legacy LP64 zheev interface and marks it deprecated. The analytic genuinely complex case and the larger residual/orthogonality checks exercise the complex ABI; no modern Accelerate ILP64 interface was tested.

The pinned installed compiler was accessed through the temporary rustup-home alias procedure already documented in `docs/environment.md`, using `/private/tmp/p0c-rustup` in this run. No global toolchain changes or downloads occurred. Offline build, linking, formatting and clippy succeeded without dependency changes or source fixes. A convenience Python metadata probe failed because the system Python lacks `tomllib`; Cargo.lock versions were instead read directly. CPU sysctl probes were denied. Neither issue affected numerical execution.

## Recommendation and limits

Adopt **faer for complex Hermitian eigensolve** behind the planned narrow backend interface. It passes the same tolerances as zheev and is faster at every requested random size in this run, while avoiding the unsafe workspace/ABI adapter in production. At N=1000, medians were 0.368617583 s for faer and 0.628785708 s for zheev; the maximum eigenvalue difference was 2.1316282072803006e-13. Keep LAPACK as an independent test oracle. This comparison is specifically with zheev; it says nothing about zheevd, zheevr, modern Accelerate interfaces, or threaded defaults.

Adopt **faer for complex GEMM provisionally**: it passes the independent reference for both ordinary and adjoint products and is faster than the triple loop at N=500 (ordinary product: 0.018948250 s versus 0.092120250 s). An Accelerate BLAS zgemm benchmark was **not** performed: the detailed GEMM step requested faer versus the naive loop, and the supplied manifest has no BLAS function-wrapper crate. Thus these measurements do not establish a faer-versus-Accelerate GEMM performance winner. No hand-written BLAS FFI was added to bridge that gap.

Adopt **ndrustfft for the 3D FFT adapter**, using its existing axis APIs and reusable handlers. Every requested grid passes round trip, Parseval and signed-mode tests, including the non-power-of-two sizes, and the tiny independent DFT agrees. This removes the need to own lane traversal/strided-copy kernels. A custom rustfft axis loop was not benchmarked; the recommendation follows correctness and reduced maintenance, not a demonstrated speed advantage over that unmeasured alternative.

Not measured: multithread scaling or actual Accelerate worker count; NEON-on versus NEON-off effects; peak memory or allocation counts; other OS/CPU/compiler combinations; rectangular GEMM or arbitrary matrix strides; alternative LAPACK eigensolvers; nearly degenerate or badly scaled spectra; FFT planning time, parallel FFT APIs, custom normalizations or larger-grid naive DFT; real DFT/SCF observables; QE comparisons. Tiny timings should not drive backend selection. These are throwaway tests, not engine changes.

## Verification and reproduction

Use the verified installed toolchain as described above, then:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
export CARGO_NET_OFFLINE=true
export VECLIB_MAXIMUM_THREADS=1
cargo build --offline --workspace
cargo test --offline --workspace
cargo test --offline --release -p dft-spike-backend -- --nocapture
cargo fmt --check
cargo clippy --offline --workspace --all-targets -- -D warnings
scripts/check.sh
```

Exit codes from this run:

| Command | Exit code |
| --- | --- |
| cargo build --offline --workspace | 0 |
| cargo test --offline --workspace | 0 |
| cargo test --offline --release -p dft-spike-backend -- --nocapture | 0 |
| cargo fmt --check | 0 |
| cargo clippy --offline --workspace --all-targets -- -D warnings | 0 |
| scripts/check.sh | 0 |

The workspace build/test and script gate logs are preserved beside the release measurements. The exhaustive debug spike test took 707.29 s; the final release spike test took 6.81 s. The debug workspace test overlapped the script verification run, so its duration is a workflow observation, not an isolated performance comparison. No tolerance was changed. Existing root Cargo.toml/Cargo.lock changes, .cargo configuration and vendor files predated this task and were left intact. The absent ticket and crate AGENTS.md were replaced only as context by the explicit user task, not created. Artifact/session registry files are outside ownership and were not modified. All task changes are uncommitted for Claude review; no git metadata was written.
