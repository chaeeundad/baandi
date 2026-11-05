# P15A: frozen-Hamiltonian stopping measurements

Source baseline: `ce64b57437e2cb7b14804a4bdd15fb252dab6286`.
The measurement sources are the uncommitted P15A changes to that commit;
`verification.json` records source and executable hashes. Native Mac, pure Rust,
one thread for replay. Timings are incidental observations, never acceptance.
No QE source was accessed. No UPF contents are stored.

## Reproduce

Save a default CLI from the baseline commit as `/tmp/p15a-before-cli`, using
`cargo build --release -p dft-cli --offline`. Then in the P15A checkout:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
P15A_REPLAY=1 P15A_DENSE_IDS=0,7,14,128,288 cargo build --release -p dft-cli --offline
cp target/release/dft-cli /tmp/p15a-oracle-cli
python3 -I tests/p15a_replay.py --binary /tmp/p15a-oracle-cli --verify-repeat
python3 -I tests/p15a_analyze.py
python3 -I tests/p15a_analysis_test.py
python3 -I tests/p15a_evidence_test.py
cargo build --release -p dft-cli --offline
python3 -I tests/p15a_verify.py
scripts/check.sh
python3 -I scripts/tests/test_readme_figures.py
```

Pass `--upf /absolute/path/Si.pz-vbc.UPF` to both runners if necessary. Synthetic
UPFs are generated in a temporary directory by the existing P10G generator;
input text, generator commit and UPF SHA-256 are recorded. The four original
fixtures retain P14A's inputs. The added slab has four atoms in a plane in an
8x8x20 Bohr cell; the 16-atom cell is 12x8x8 Bohr at Gamma and a 2x1x1 mesh.
All synthetic fixtures use valence 2, channels [0,0,1], soft shape, cutoff 3 Ry.
The slab and 16-atom mesh exhaust 100 SCF iterations before and after; they
are diagnostic trajectories, not certified ground-state material results.

## Files and semantics

- `tables.md` and `analysis.json`: generated stop scans, matched accuracy,
  per-root rate curves, schedules, first redo flags, and selected snapshots.
- `<case>.json`: exact input, commands, hashes, both full SCF histories,
  dr2 in Ry, numerical equivalence and repeat checks.
- `<case>-traces.json.gz`: deterministic gzip of a JSON array indexed by solve
  completion order. Each solve contains variants 0 (native), 1..5 (residual
  thresholds 1e-2..1e-6 Ha), and 6 (all roots 1e-10 Ha, true-residual check).
  `dense_eigenvalues_ha`, when present, is an independent full-space oracle.
  Decompress with Python's standard `gzip` module or `gzip -dc`.
- `default-before/after-t*.json`, `identity-commands.json`: unchanged bench_scf
  identity checks for five cases at T=1,2,4,8,12. Earlier `before.json` and
  `after.json` retain the initial T=1 check. Timings overlap other work.

Every replay borrows the same Hamiltonian object, including its frozen local
potential, kinetic/nonlocal data and cell, and the identical starting columns.
It runs synchronously before SCF resumes. Only variant 0 is returned to SCF.
There is no approximate serialization of a Hamiltonian, external restart file,
or change of the trajectory feeding the next SCF iteration. Each trace reserves
its step/root arrays before the eigensolver loop. No trace printing, timers or
file I/O occur inside the numerical loops. Output drains through the existing
post-SCF CLI accounting boundary. An ordinary build monomorphizes the solver
with tracing off and omits all measurement replays and density diagnostics.

Ritz-change rules are **offline prefix replays**: the underlying correction,
preconditioner and restart choices stay fixed on each measured trajectory.
This isolates stopping from convergence rate; it is not an implementation of
QE's different direction selection or its generalized nonorthonormal basis.
Independent residual-tolerance replays do change correction selection through
the existing threshold logic; they are recorded as separate trajectories.

Ritz thresholds scan 1e-2 through 1e-13 Ha; the same ladder is scanned for
residual prefix stops. Allowed unconverged roots scan 0,1,2,5. The early-allowance
experiment differs from QE's allowance at its 20-expansion limit, so both are
reported separately. Initial Ritz change is unavailable, not zero. Sorted root
indices can change character across rotations; first-below is not locking.

A missing stop is null (censored), not the terminal step. Work is cumulative
H vector applications, H blocks, restarts and basis width, independent of wall
time. Dense-oracle work and discarded replay work are excluded from native SCF
counts. Native steps include the initial reduced solve; QE's expansion counter
does not. The ordinary/certificate partition is verified against the history.

The density diagnostic uses the normalized Fourier coefficients of the signed
output-minus-input density, the density sphere and no G=0 term:
`dr2_Ry = Omega * sum(4*pi*|delta_n_G|^2/G^2)`. The offline schedule starts at
0.005 Ha, ratchets with `0.05*dr2_Ry/max(1,N_e)`, and floors at 5e-14 Ha.
First-iteration redo is flagged using the native stored residual; its changed
vectors/residual and subsequent SCF trajectory are not manufactured. Schedules
on these fixed histories are conditional diagnostics, not predictions of a
new self-consistent run.

See [the interpretation](../../../../docs/research/qe-study/replay-stopping-results.md)
for accuracy limitations, spectral completeness findings and ranked experiments.
