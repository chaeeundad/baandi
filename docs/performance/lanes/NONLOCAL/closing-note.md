# NONLOCAL Phase 2 closing note

2026-11-20 — recommend closing this lane for the current optimization round.

The relayed dev08 `bad02ee` Si64 np8 profile covers 0.7825 s of 9.6702 s.
Projection/back wrappers reach about 34/28 GFLOP/s; the 64–108 GFLOP/s
reference leaves theoretical headroom, but does not identify a change that
captures it. The source-derived 256-panel model already uses real DGEMMs.
Complex-to-real packing is constructor work, and the MPI route does not use
the serial MKL_NONLOCAL switch.

Two concrete shape-preserving candidates fail the requested 0.15 s expected
gain screen: direct output/add fusion has an optimistic bandwidth-scale
estimate around 0.050 s; eliminating coupling entirely saves at most its
0.0474 s attributed stage. Persistent backend packing has no measured cost
attribution supporting the threshold. Reduction is wait-dominated and its
payload experiment already failed. No rejected traversal, tile narrowing or
column grouping is repeated.

See [diagnose.md](diagnose.md) for arithmetic, bandwidth, shape assumptions,
candidate estimates and the distinction between reference floors and expected
gain. No optimization, switch, default, numerical path or instrumentation was
changed in Phase 2. No dev08 A/B is requested. No new Mac correctness test is
needed for these documentation-only changes; Phase 0 evidence is not presented
as a Phase 2 test. Changes are left uncommitted.

Reopen if same-shape dev08 evidence identifies a removable cost sufficient for
>=0.15 s expected Si64 np8 wall saving, or if the orchestrator changes this
round's threshold. A future implementation would require the prescribed
default-off switch, numerical equivalence and crossed A/B |mean|/SE >=3.
