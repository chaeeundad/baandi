# LANE-LAUNCHIO Phase 2 diagnosis (2026-11-20)

Source: user/master relay of dev08 commit **e947a35**, five cells, median of
**six on-* runs**, instrumented, span values rank maxima. Raw result files,
executable byte length, grid sizes and hardware bandwidth measurements were not
included in the relay. Do not manufacture them or treat this as reproduced data.
Current implementation base is `0413485` (Phase 0). Spans nest and overlap;
none of the columns below is summed into a wall-time closure.

| Span, seconds | Si64 MPI8 | Si32 MPI8 | Si16 MPI8 | Si64 one core | Si32 one core |
| --- | ---: | ---: | ---: | ---: | ---: |
| External minus max main | .1245 | .0824 | .0689 | .0899 | .0414 |
| SCF driver | 9.5205 | 1.7879 | .5901 | 44.5407 | 7.8312 |
| Materialize including outputs | .23333 | .11964 | .06132 | .05312 | .02864 |
| density.bin write | .20749 | .10612 | .05372 | .03186 | .01699 |
| Engine construction | .19118 | .09805 | .05572 | .23107 | .10852 |
| Binary SHA-256 | .10446 | .10474 | .10456 | .09387 | .09424 |
| Preflight including preparation/UPF | .05849 | .02682 | .01284 | .04131 | .02139 |
| MPI finalize | .05533 | .06176 | .05370 | .00158 | .00157 |
| results.json write | .02934 | .02718 | .02723 | .00492 | .00440 |

The supplied comparison says Si64 MPI8 outside-process cost is already about
QE's 0.12 s. No new QE run or cross-code ratio is claimed. This replaces the
older approximate 0.84 s attribution. Work should target **inside-process** work.

## Source attribution and floors

Let `D = launchio.sizes.density_bytes = 8 * Nx * Ny * Nz` and `N` be the actual
executable file length, both in bytes. The raw counts must accompany the next
relay; the following formulas are exact workload counts, not assumed bandwidth.

- `session.rs` materializes an owned terminal density through
  `ScfProblem::materialize_owned` / `grid_dist::Plan::gather_field` before returning
  the outcome. This is an **allgatherv to every rank**, not a gather to root.
  Logical delivered-density receive payload per rank is D; off-rank delivered
  bytes aggregate to `(P-1)*D` for one full-density gather, ignoring wire protocol
  and implementation routing. Existing handoff/verification gathers are separate.
- CLI converts the already-replicated density into a contiguous little-endian
  Vec (D input bytes read, D output bytes stored), hashes it (D read), then
  `density_bin_write` encloses only `fs::write(path, density_bytes)` and the moved
  Vec's drop. There is **no MPI collective or layout conversion inside it**.
  Rust's write covers open/create/truncate, write_all and close; no fsync is used.
- Every rank writes its own complete density.bin (`rankN/` except root): **P*D**
  logical file bytes in total, D per file. The OS write copies D bytes into the
  page cache, so at least D user reads plus D destination writes per rank occur.
  Physical disk bytes and persistence latency are unknown without filesystem/
  page-cache telemetry; a successful non-fsync write is not a durable-media rate.
- Define measured per-rank payload throughput `D/t` and aggregate throughput
  `P*D/t_max` for a synchronized output interval (the latter is only a rate proxy
  if rank starts differ). Si64 rates are `D/.20749` and `8D/.20749` bytes/s versus
  `D/.03186` serial. Si32: `D/.10612`, `8D/.10612` versus `D/.01699`.
  At equal D the slowdowns are **6.5126x** and **6.2454x**. Eight full-file writers
  are consistent with contention; these data do not establish whether metadata,
  writeback throttling, storage or allocator/drop time dominates. Do not claim
  rank scaling is a numerical/MPI-gather failure from this timer.
- Given independently measured sustained memory bandwidth Bm and shared output
  bandwidth Bw (bytes/s), ideal write lower bounds are `2*P*D/Bm` for the explicit
  copy traffic and `P*D/Bw` for aggregate payload drainage (if the latter lies on
  the timed path). Metadata/latency adds an unknown floor. Do not divide by a
  guessed SSD/DRAM bandwidth to produce a purported measured saving.
- The fresh hash reads **N** bytes once per process (P*N across P ranks), with
  SHA-256 block count `ceil((N+9)/64)` and 64 rounds/block. With hash throughput
  Bh and effective read bandwidth Br, its lower bound is `max(N/Bh, N/Br)` plus
  file/metadata overhead. The observed per-rank effective rate is `N/.10446`
  bytes/s on Si64 MPI8; without N no numerical GB/s is defensible. The existing
  64 KiB streaming implementation has already removed the executable-sized
  allocation; shrinking scratch again does not remove the digest scan.
- A valid persistent hash-cache hit reads at most **4096** sidecar bytes and
  **zero executable payload bytes** for digest computation. It still opens the
  executable and validates metadata twice. All ranks can use the same immutable
  cache record; no barrier or communication is necessary.

`results_write` has 41–69% asymmetry, density_json_write up to 92.7%, and MPI
finalize 18.9–30.6%: waits and output scheduling can transfer cost between spans.
Treat them as correlated diagnostic observations; **external whole-run wall**
remains the optimization verdict. Construction is not optimized in this lane.

## Candidates and decision

1. **Selected: default-off persistent executable SHA-256 cache.** Keep the exact
   path and digest on every run, including measurement jobs. Key by Unix device,
   inode, length, mtime+ns, ctime+ns, uid and mode. Validate an already-open file
   and the path identity before accepting/publishing. Cold/malformed/stale/
   unavailable caches compute the original streaming digest. Expected warm saving
   is a substantial part of 0.094–0.105 s minus metadata/cache costs; planning
   range **0.05–0.10 s**, not measured. Supplied wall shares are 10.56% Si16 MPI8,
   4.53% Si32 MPI8, 1.01% Si64 MPI8 and 1.2% Si64 one core. Cold runs may be slower.
   F-LAUNCH-IO previously deferred caching because it could stale provenance.
   New premise: real dev08 fixed-cost evidence and explicit user authorization
   to cache, with metadata invalidation and independent master checksum gates.
2. **Deferred: root-write plus hard-link/reflink publication of identical density
   files to rank directories, or staggered writes.** Possible Si64 ceiling is
   roughly 0.20749 s, and a serial-write-like cost would save at most about
   0.17563 s in that span before coordination/metadata. This is **not** a whole-run
   prediction. Filesystem locality and link support must be proven; hard links
   change later mutation semantics, rank-zero failure becomes shared, and ranks
   need synchronization. Merely root-only output removes existing files and
   violates the requested byte-preserving output contract. Not implemented.
3. **Deferred: hardware-accelerated fresh SHA.** Retains content scanning and
   avoids persistent trust but needs confirmed CPU/library support and dependency
   ownership. No source or build-feature change is justified without that evidence.

Do not omit provenance, move hashing outside timed B only, or change the core
count to hide it. Only candidate 1 is implemented. No Mac profiling or timing
conclusion. The exact validation and trust limits are in dev08-request.md.
