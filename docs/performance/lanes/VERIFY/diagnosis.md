# LANE-VERIFY Phase 0: code map only

The revised user instruction of 2026-11-20 supersedes the original lane charter.
Status: PROFILE_REQUEST. No optimization or candidate selection is retained.
Mac runs validate correctness and instrumentation plumbing only. Earlier scratch
timing experiments predate the instruction and are excluded from all performance
conclusions and from the profile request.

`session::verify_fixed_point` temporarily installs an owned input (cloning the
owned density), invokes the inner verification, then restores the original input.
The inner routine clones warm matrices/lazy handles, runs one ordinary map,
evaluates the Coulomb metric and residual/density/eigenvalue metrics, makes the
existing decision, and constructs certificate/attempt records. A passing
verification delivers the preceding SCF state. Failed-attempt reuse and scheduling
are unchanged. The inner map already has potential/solve/density/energy phase
buckets; this lane adds per-attempt, per-rank attribution around that map.

Ledger read: VERIFY-REUSE variants and adopted HALF, VERIFY-SCHEDULE and
G-SCHED-OVERHEAD. No rejected policy is repeated. Judgment rule, tolerances,
map count and numerical work are unchanged. See profile-request.md for the
requested dev08 measurements; only after that relay may Phase 2 begin.
