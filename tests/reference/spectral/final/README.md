# Final-policy follow-up evidence

Current outcome: explicit runtime `[scf] spectral_guard = "final"`, default off
because si8/k444s H work exceeds 3%. See
[report](../../../../docs/validation/spectral-completeness.md) and
[dev08 request](../../../../docs/validation/spectral-dev08-request.md).

- `summary.json`: generated comparisons, cost, material spectra, dense skew trajectory.
- `latest-t1/`, `t2/`, `t4/`, `t8/`: current production executable, seven cases,
  both policies; full results and commands. `t1/` is an earlier implementation snapshot.
- `materials/`: all six P13D cutoff fixtures, ordinary/fresh final dense comparisons.
- `dense-scf/`: diagnostic dense-every-step skew trajectories, including exhaustion.
- `default-bench.json`: original five-case numerical-equivalence gate.
- `robustness.json`: ordinary default robustness gate (16 pass, 3 nonconverged).
- `verification.json`, `logs/`: current verification and reviewed Rust source hashes.

Root-level historical candidate records and BLOCKED manifest are preserved;
this directory supersedes their implementation-status conclusions. Mac timing
samples are observations, not isolated performance acceptance measurements.
