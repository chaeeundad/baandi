# ADR 0053: Single entry (`[task]`), progress stream and job monitor (UX1)

- Status: **Proposed** (ticket UX1-SINGLE-ENTRY-MONITOR, 2026-11-27; orchestrator review pending).
- Scope: `crates/dft-cli/src/{task.rs, progress.rs, monitor.rs, monitor/tui.rs}`,
  the command-line front end in `main.rs`, the `Config.task` field, boundary
  hooks in `run.rs`, `relax.rs`, `bands.rs`, `telemetry.rs`; user documentation
  in [CLI usage](../cli.md#running-and-monitoring).
- User request (2026-11-27): one executable, VASP-style, with the job chosen
  by the input; a modern terminal view of a running job's state.

## Decision

1. **One entry.** `baandi` with no subcommand reads `input.toml` in the
   working directory (or `baandi -i FILE`) and runs the job of its `[task]`
   table: `type = "scf"` (default when `[task]` is absent), `"relax"`,
   `"cell-relax"`, `"bands"`, `"dos"`, `"nscf"`; `resume = true` (relax,
   cell-relax) and `checkpoint_dir` (nscf) are the only other keys. The bare
   entry is translated into the existing alias subcommand (`run`, `relax`,
   `bands`, `dos`, `nscf`) before anything else happens, so the job runs
   exactly the code it ran before. The subcommands stay, unchanged, for
   scripts and tests; an explicit subcommand ignores `[task]`.
2. **Selection rule.** `task.type` is the only selector. A `[relax]` table
   does not imply relaxation: on the bare entry, `[relax]` without
   `type = "relax"`/`"cell-relax"` is an input error, as is `relax` with
   `relax.mode` other than `ions`, `cell-relax` with `relax.mode = "ions"`,
   `bands`/`nscf` without their table, and a key that does not apply to the
   type. Nothing is silently ignored. The output directory is `[output]
   directory` (default `results`, relative to the working directory), the
   same rule every driver uses; relax, which requires `--output-dir` as a
   subcommand, receives it from there.
3. **Progress stream `baandi-progress-1`.** Every job command (bare or
   alias) writes `progress.jsonl` in its output directory: append-only JSON
   lines, one per boundary (job start, phase change, SCF start, each
   completed SCF map, SCF end, each relaxation frame, warnings, job end),
   written and flushed directly, never inside a map, and only on MPI rank 0.
   Envelope: `v` (=1), `seq` (1, 2, ...), `t` (Unix seconds), `elapsed_s`
   (since job start), `event`. Event fields:
   - `start`: `schema`, `task`, `input`, `output_dir`, `pid`, `host`,
     `version`, `system` {`name`, `formula`, `atoms`, `kpoints`, `nbands`
     (null = automatic), `ecutwfc_ry`, `volume_bohr3`, `spin`, `symmetry`,
     `threads`, `ranks`, `relax` {`mode`, `force_tol_ha_bohr`, `energy_tol_ha`,
     `stress_tol_ha_bohr3`, `max_steps`} or null} (null when the input does
     not parse; the job then fails with its usual error).
   - `phase`: `phase` = `setup` | `scf` | `post-scf` (forces, stress, outputs).
   - `scf_start`: `cycle` (1 for a single SCF, one per relaxation step),
     `threshold_ry`, `threshold_ha`.
   - `scf_map`: `cycle`, `iteration`, `kind` (`ordinary` | `verification` |
     `retry`), `stage`, `energy_ha` (null while the driver has not evaluated
     it), `de_ha` (change from the previous map with an energy), `dr2_ry`,
     `estimator_ha`, `ethr_ry`, `h_columns`, `reused`, `seconds` (wall time
     of the map).
   - `scf_end`: `cycle`, `maps`, `converged` (null when the driver result
     does not carry it), `error`.
   - `relax_step`: `step`, `energy_ha`, `max_force_ha_bohr`,
     `step_length_bohr`, `scf_iterations`, `provenance`, `pressure_ha_bohr3`
     (cell modes), `status`.
   - `warning`: `source`, `message` (stall diagnostics, failed relaxation
     attempts, a replaced stale lock).
   - `end`: `status` (`completed` | `not_converged` | `failed` | `aborted`),
     `exit_code`, `message`.
   Readers ignore unknown events and fields; a new field is a compatible
   change, a changed meaning bumps `v` and the schema name. A previous stream
   in the same directory is kept as `progress.prev.jsonl` (one generation).
   `--no-progress` turns the stream and the lock off. A job rejected before it
   wrote anything else (an input error) removes its stream again, and the
   output directory when the session created it, so the existing contract
   "an invalid input creates no output directory" holds.
4. **Lock.** `baandi.lock` (`baandi-lock-1`: `pid`, `host`, `started_unix`,
   `task`, `input`) is created exclusively before the job and removed at its
   end. A lock whose process runs on this host, or a lock from another host
   (cannot be checked), refuses a second job in the directory (exit 5) and
   names `baandi watch`; a lock whose process is gone (checked with `ps -p`,
   no signals, no FFI) is stale, is replaced, and the replacement is
   recorded as a warning.
5. **Monitor.** `baandi watch [DIR]`, and a bare `baandi` whose output
   directory holds a live lock (or a bare `baandi` run inside an output
   directory), tails `progress.jsonl` (no IPC) and shows task and system,
   phase, SCF convergence (log10 dr2 and |dE| per map with the threshold
   line), relaxation trajectory (energy and log10 max force per step with
   the force tolerance), per-map and per-step times, elapsed time, events
   and warnings, and the final status. On a terminal it is a ratatui
   full-screen view (q quits, Tab switches Overview / SCF maps / Relaxation
   / Events); with `--plain`, or when stdout or stdin is not a terminal, it
   prints a text view and reprints it when the stream changes until the job
   ends (`--once` prints once). A finished directory shows its final summary.
6. **Numerical identity.** The reference SCF calls go through
   `progress::solve`, which calls `refscf::api::solve` when no stream is
   active and otherwise performs the identical sequence `stepper`,
   `initialize`, `step` until done, `finish` (the body of
   `Stepper::run_to_end`) with a record after each map. No record feeds
   back; no write happens inside a map; a failed write disables the stream
   with one warning and the job continues.

## Consequences

- New dependency ratatui 0.30 (MIT, default features off, `crossterm`
  backend; crossterm 0.29 MIT through ratatui-crossterm), see
  [dependencies](../dependencies.md). The UI code lives only in `dft-cli`.
- The strict-oracle and finite-temperature/spin oracle drivers (not the
  reference path) produce `start`/`phase`/`relax_step`/`end` records but no
  per-map records. `run --observe`/`--agent` add per-map records from the
  runtime loop in addition to their unchanged `events.jsonl`.
- Lock liveness by process id can be fooled by pid reuse after a crash; the
  monitor then shows `RUNNING` until the reused process ends. Remove the lock
  by hand in that case.
- The `[task]` table becomes part of the input text and thus of the input
  SHA-256 recorded in results (it is not part of any model identity).

## Evidence

`cargo test -p dft-cli --lib -- task:: progress:: monitor::` (task mapping
and rejection table, lock live/remote/stale, envelope and dedupe, recorded
relaxation stream folded and rendered against a committed plain-view
snapshot, torn last line, every TUI panel rendered at 120x40 and 30x10 on
ratatui's test backend) and `cargo test -p dft-cli --test all single_entry::`
(bare entry = `run`: density bytes and total-energy bits identical with the
stream, without it (`--no-progress`) and through the alias; live lock refuses
a second job (exit 5) and turns the bare entry into the monitor; input
errors).
