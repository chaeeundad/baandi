# QE study note 2: FFT call structure in the wavefunction path (from the fastqe Phase 1 report)

Status: analysis note written by Claude Code in its own words from a prose report of the user's own earlier QE acceleration project (the fastqe session; QE 7.6, MKL DFTI backend, CPU), received on 2026-11-03. No code or excerpts. Policy: AGENTS.md, "QE source study policy". The statements are the report's measurements on other hardware (Xeon w5-3435X, 72^3 grids and larger); every hypothesis must be re-measured on our engine.

## 1. What the report found
1. The baseline problem: with the MKL DFTI backend QE transformed the whole 3D box for every band (the FFTW backend pruned zeros, DFTI did not), and copy, zero, scatter, multiply and gather were separate full-array passes.
2. The pruned and fused pipeline for one band, inverse direction:
   - scatter the plane-wave coefficients into sticks, the (x,y) columns that contain at least one point of the wavefunction sphere (about 20 percent of the columns at the usual density cutoff of four times the wavefunction cutoff), stored plane-major;
   - one batched 1D FFT call along z over all sticks (number of transforms = number of sticks, stride and distance chosen so that the transforms vectorize across sticks);
   - then, one z plane at a time while the plane stays in cache: scatter that plane's sticks into a plane buffer, 1D FFT in y only on the x columns that contain sticks (about 50 percent), 1D FFT in x on the whole plane out of place, multiply by the local potential (or accumulate the weighted squared modulus for the density), x transform back, y transform back on the active columns only, gather the sticks;
   - forward direction pruned the same way because only the sphere values are needed afterwards.
   The full 3D box is never formed; flops fall by about 40 percent and memory passes by much more. Because the x pass is out of place, the inactive columns of the plane buffer stay zero forever: no per-plane zeroing, only active columns are rewritten.
3. Gamma point: two real bands share one complex FFT; bands with exactly zero occupation weight are skipped in the density (bit-identical).
4. Pruned plane transform (x,y) only when at most about 60 percent of the x columns are active; otherwise the vendor's native 2D transform is faster (strided 1D y transforms are 1.3 to 1.4 times less efficient than the 2D kernel). The wavefunction grid has about 50 percent active columns, the density and potential grids about 100 percent.
5. Measured (their node, 72^3, 64 bands): per band 2.3 ms originally, 0.46 ms for the first pruned and fused version, 0.31 ms with plane-major sticks and no zeroing; the transform error against the full 3D one about 2e-16. Whole runs with one process and four threads: the local-potential step 3.0 to 5.0 times faster, whole run 1.11 to 1.87 times; energies equal to 8 decimals.
6. Vendor descriptor settings that worked: 1D, double complex; number of transforms equals the number of sticks or active columns; input and output distance and strides chosen per pass (z: stride 1 with distance equal to the padded stick count; y: stride equal to the row length, distance 1); out-of-place only for the x pass, in-place elsewhere; thread limit 1 (sequential descriptors, parallelism from outside); forward scale 1/N applied only where needed; descriptors cached per thread in a small cache keyed by (length, count, stride, distance, placement); creating and committing a descriptor costs about 21 microseconds, so sharing descriptors across k points by rounding sizes was not worth it.
7. Threads: the honest comparison is one thread per process with pinned cores; threaded runs needed pinning; an idle node is mandatory; vendor-internal threading on small strided batches was about 3 times slower.

## 2. What did not help (measured by them)
Grouping 2 to 8 bands per pass over the planes so that the potential is read once per group (+-1 percent; the multiply is compute-bound, not memory-bound); batching several planes into one vendor call (call overhead negligible); the code's own "many bands" batched 3D mode (18 to 55 percent slower than per-band pruned and fused); batching bands into one all-to-all (irrelevant for us); mixed-precision FFT (not done, at most about 15 percent of that step, changes the SCF path). Small grids gain least from plane fusion (the whole array fits in cache); for a 25^3 box (about 250 KB) they expect the pruning (about 40 percent fewer flops) to matter more than the fusion.

## 3. Differences from rust-dft today (verified in our code)
- We already skip zero lanes (O3), but through a full-grid array: every application zero-fills the whole grid (15,625 complex numbers), scatters, and the transform copies the input grid into the output grid before working in place; the second transform copies again. The density path (`dft-physics/src/density.rs`) uses the full transform per band, not the pruned one.
- One strided lane at a time per axis; no batched 1D calls (rustfft has none; MKL has them).
- Bands are applied one vector after another; the local-potential multiply is a separate pass over the full grid.
- The Gamma-only real-wavefunction pairing is not available: the engine uses complex wavefunctions at every k including Gamma by design (validation against QE complex Gamma); it would apply only to a future dedicated Gamma-real path.

## 4. Hypotheses for rust-dft (to measure; the vendor backend is ADR 0019)
S-A. Stick layout for the wavefunction transform in the vendor backend: one batched z transform over the sticks, then plane-by-plane y (active columns) and x passes with the potential multiply fused, as in item 2; descriptors owned by the plan (no per-thread global cache), sequential threading layer.
S-B. Remove the full-grid zero-fill, the scatter into a full grid and the in-transform copy in the pure Rust path too: the same plane-buffer technique works with rustfft lanes (a plane buffer with permanently zero inactive columns), bit-identical to the current transforms.
S-C. Density accumulation through the same pruned and fused path (inverse only, then add the weighted squared modulus), skipping bands with zero weight.
S-D. Keep per-band processing (no "many bands" 3D batching); keep sequential vendor descriptors with our own outer threads (consistent with P7P6, P7PB).
S-E. Pruned 2D only when the active column share is at most about 60 percent; use the vendor 2D or 3D transform otherwise (density, potential, Hartree grids).
Measurement protocol: pinned cores on the idle node, compare per-band transform time at 25^3 and at a larger grid (the 8-atom cell), then the whole SCF with the benchmark cases; gates as in docs/performance/README.md.
