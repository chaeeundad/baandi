TASK_ID: P10H-MKL-CALL-PATH
GOAL: Close the gap between the Rust MKL feature path and the C reference MKL kernels found on dev08 (2026-11-04, commit f10cc4e, 1 core): fft3d 1.3-1.7x slower than C, zgemm 750xN only 1.6x faster than pure Rust versus 4.1x for the C reference, zheev not routed through MKL at all.
OWNED_FILES: crates/dft-mkl/src, crates/dft-numerics/src (gemm and eigen adapters only), crates/dft-numerics/examples/lib_bench.rs, docs/performance/optimization-log.md
DEPENDENCIES: MKL link fix merged (e4660c3) and confirmed on dev08; the remote measurement is done by dft5090 and re-checked by Claude.
STEPS: profile the call path for zgemm (layout conversion, allocation, extra copies, transposes, threading layer) and fft3d (plan caching, in-place versus out-of-place, alignment); route the Hermitian eigenproblem (zheev or zheevd) through MKL behind the same feature with the same accuracy gate; keep pure Rust default and bit-identical behaviour per backend across thread counts.
GATES: kernel table (lib_bench) before and after on dev08; total energy difference pure Rust versus MKL under 1e-12 Ha on the protocol cases; density hash identical across threads within each backend.
NON_GOALS: new backends, GPU, changing defaults.
GIT: branch task/P10H-MKL-CALL-PATH, uncommitted
