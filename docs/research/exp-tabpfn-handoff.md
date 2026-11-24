# Handoff: EXP-TABPFN-SCF (active SCF with a TabPFN controller)

Written 2026-11-25 by the orchestrating session (main performance work stays there). This experiment is
handed to a separate session at the user's request. It is research only: nothing here is merged into
`main` or used as a default without a separate user decision. Not a speed measurement, so it runs on the
Mac (user, 2026-11-25); keep load modest (the Mac also controls dev08 work).

## Goal (user, 2026-11-25)
Observe the reference SCF trajectory map by map with TabPFN and let it choose the next mixing action
("active SCF"); report how the physics (delivered fixed point, energy, density, verification) and the
computation (maps, H columns, wall incl. controller time) change against the fixed rule (modified
Broyden, beta 0.7, history 8).

## Where everything is
| What | Path |
|---|---|
| Rust code (hook, trace, exploration) | branch `task/EXP-TABPFN-SCF`, worktree `/Volumes/990-2TB/Develop/rust-dft-wt/EXP-TABPFN-SCF` (committed on the branch, not pushed, not merged) |
| Protocol and trace schema | `docs/research/exp-tabpfn-scf.md` in that worktree |
| Binary used | `/Volumes/990-2TB/Develop/rust-dft-wt/EXP-TABPFN-SCF/target/release/dft-cli` (Mac release, pure build: no MKL, no MPI) |
| Experiment directory | `/Volumes/990-2TB/Develop/rust-dft-wt/exp-tabpfn/` |
| Python env (3.12, tabpfn 10.2.0, torch 2.15) | `/Volumes/990-2TB/Develop/rust-dft-wt/tabpfn-env/venv` |
| Weights cache | `tabpfn-env/models` (v2: Apache-2.0 with attribution; v3.5: non-commercial, token in `tabpfn-env/token.sh`, mode 600, never print or commit it) |
| UPF | `/Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF` (never commit UPF files) |

Environment for any TabPFN process:
`source tabpfn-env/token.sh; export HF_HOME=tabpfn-env/hf TABPFN_MODEL_CACHE_DIR=tabpfn-env/models TMPDIR=tabpfn-env/tmp TABPFN_DISABLE_TELEMETRY=1`

## Switches on the experiment branch (all default off)
- `EXP_SCF_TRACE=on`: writes `<out>/scf_trace.jsonl` (one completed line per outer map; features must be
  taken from the embedded `observation` object to avoid leaking post-action values).
- `EXP_SCF_CONTROLLER=<executable>`: one child per run; per map one observation line on stdin, one reply
  `{"beta": (0,1.5], "reset_history": bool}` on stdout; 30 s exchange deadline; failure falls back to the
  fixed rule. Stopping estimator, ethr ratchet and the fixed-point verification criterion are unchanged.
- `EXP_SCF_BETA_EXPLORE=<seed>`: beta drawn uniformly in [0.2, 1.0] per map (training data).
Requires `driver='reference'`, one rank, `verify='fixed_point'`.

## Tools (`exp-tabpfn/tools/`)
- `cases.py`: input generator (si8, si16 = 2x1x1, skew = sheared 2-atom, elong = 3x1x1 conventional),
  seeded random displacements (Bohr amplitude), 24 Ry, LDA-PZ, occupied + 4 bands.
- `run.py <outdir> <kind> <amp> <seed> [K=V ...]`: one run; writes `in.toml`, `out/`, `run.json`, `stdout.txt`;
  sets `CTL_LOG=<outdir>/ctl.jsonl`. Needs `SI_UPF` and `DFT_CLI` in the environment.
- `features.py`: 10 features per (observation, candidate beta): iteration, log10 dr2, previous log ratio,
  previous beta, history occupancy, log ethr, log max Ritz change, log |dE_total|, Davidson iterations, beta.
- `table.py <runs dir> <out.json>`: training table, target y = log10(dr2_{k+1} / dr2_k); converged runs only.
- `cv.py <train.json> <v2|v3.5>`: grouped 4-fold CV (by geometry) against linear and constant baselines.
- `server.py <train.json> <v2|v3.5|linear> <socket>`: persistent model (fit once); greedy one-step policy:
  pick beta in {0.2, 0.3, ..., 1.2} minimizing predicted y. Logs predictions per map to `CTL_LOG`.
- `client.py` + `ctl.sh`: thin stdin/stdout bridge to the server (`CTL_SOCK=<socket>`); this is what
  `EXP_SCF_CONTROLLER` points to. The server accepts one connection at a time: run controller arms with
  `xargs -P 1` per server.
- Servers may still be running from the last session (sockets `/tmp/claude-501/tabpfn-{v2,v3.5,linear}.sock`,
  logs `server-*.log`). Check with `pgrep -fl tools/server.py`; stop them when done.

## Results so far
Training: 88 runs (si8 x10, si16 x6, skew x6 geometries; 1 default + 3 exploration seeds each);
amplitude 0.2 Bohr never converged (gap closes; closed-shell fixed occupation), so 54 converged runs,
857 rows (`train.json`).

Offline CV, RMSE of log10(dr2 ratio): TabPFN v3.5 0.528, v2 0.542, linear 0.617, constant 0.645.

Held-out test (`runs/test/`, summary `SUMMARY.txt`), 15 geometries x 4 arms (default, linear, v2, v3.5).
skew a0.12 failed in every arm (gap). Over the 14 converged geometries:

| arm | maps | H columns | wall (s) | controller (s) |
|---|---|---|---|---|
| default | 194 | 19,925 | 23.6 | 0 |
| linear | 261 | 20,671 | 24.5 | 0.6 |
| TabPFN v2 | 195 | 18,057 (-9%) | 1037.5 | 1008.4 |
| TabPFN v3.5 | 197 | 19,296 (-3%) | 1276.6 | 1253.1 |

Physics: every converged arm reaches the same fixed point; max |E - E_default| = 3.7e-10 Ha (standard
class 5e-9 Ha), all verified. The path changes, the answer does not.
Computation: TabPFN roughly matches the fixed rule in maps (wins and losses per geometry, e.g. si16 a0.02
14 -> 12, si8 a0.02 s203 14 -> 21 with v2) and lowers H columns modestly; CPU inference costs ~40x the SCF.
The unseen elongated cell: v2 helped on two of three geometries (16 -> 13, 17 -> 14 maps).

## Open next steps (proposals, not decisions)
1. Better policy than greedy one-step: target maps-to-convergence or a multi-step objective; add an
   "uncertainty" guard (fall back to 0.7 when the prediction interval is wide).
2. Harder systems where mixing matters: Si vacancy / divacancy / interstitial, strongly disordered si64,
   6x1x1 long cells (charge sloshing). Multi-element systems need either LDA pseudopotentials (only Si is
   available) or PBE support in the reference driver (`refstart.rs` rejects PBE); the user deferred that
   choice: do not run PBE-generated files under LDA without asking.
3. Cost: batch the beta grid in one predict call (already done), reduce context size, or distill the
   TabPFN policy into a cheap rule; report controller time separately from SCF time in every table.
4. Write results into `docs/research/exp-tabpfn-scf.md` on the experiment branch (English), and give the
   user Korean plain-language summaries.

## Rules that still apply here
English for files; Korean plain language to the user. Never commit UPF files or tokens. Do not touch
`main`, `scripts/delegate.sh`, or dev08 queues (owned by the orchestrating session). Builds and scratch on
`/Volumes/990-2TB` only.
