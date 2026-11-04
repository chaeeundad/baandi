# CLI execution modes

The existing executable is `dft-cli`. `baandi` is a thin POSIX shell launcher
(`crates/dft-cli/src/baandi`) that executes the same binary; both retain `run`,
`validate`, `grid`, `inspect-upf`, `pseudo`, `bands`, `dos`, `eval`, and `convert`.
Build with `cargo build -p dft-cli`. For local use:

```
export PATH="$PWD/target/debug:$PWD/crates/dft-cli/src:$PATH"
baandi --help
```

For installation, copy the launcher alongside `dft-cli` in the chosen executable
directory. It first uses a sibling executable, then searches PATH. There is still
one Cargo binary, preserving existing `cargo run -p dft-cli -- ...` calls. A second
Cargo target would require `default-run = "dft-cli"` in Cargo.toml, which is outside
this ticket's ownership. On platforms without a POSIX shell, use `dft-cli` directly.

## Classic

```
baandi run input.toml --output-dir results
```

No observation or control is enabled. The existing `ScfProblem::run_outcome`
path, final outputs, console text and exit codes are retained. The original
`dft-cli run` spelling continues to work. Execution mode does not select physics,
change tolerances, or provide a new convergence criterion.

## Observed

```
baandi run input.toml --observe --output-dir observed
baandi run input.toml --observe --events - --output-dir observed-stdout
baandi run input.toml --observe --events telemetry.jsonl --output-dir observed-file
baandi run input.toml --observe --events unix:/tmp/listener.sock --output-dir observed-socket
```

Each observed calculation creates `OUTPUT/events.jsonl` exclusively. An existing
log is an error, never truncated or silently appended to another calculation.
Use a fresh output directory per calculation. `--events` adds a mirror: `-` means
JSON-lines stdout, a path means another exclusively created file, and `unix:PATH`
means a pre-bound Unix **datagram** listener (not a stream socket). With stdout
telemetry, preflight and final human output go to stderr. Human telemetry always
goes to stderr, reporting iteration, energy in Ha, normalized density residual
in Bohr^-3/2 and `solver=iteration_completed`. This status means an iteration
finished; only a `converged` event establishes engine certification. Certificate
pass/fail and residuals are separate events.

Every record has schema version `1.0`, calculation ID, increasing `sequence`,
Unix timestamp in seconds, iteration, state version, `type`, and a `data` object.
Event types are `initialized`, `iteration_finished`, `certificate_evaluated`,
`mixer_reset`, `diverged`, `converged`, `not_converged`, `failed`, `events_dropped`.
Field names encode units; relative charge error is dimensionless. JSON contains
scalars, never density/orbital arrays. IDs identify executions, not physical models.

The driver drains engine events after each complete step. A bounded 256-record
worker queue moves formatting and I/O outside SCF; publication never waits for a
subscriber. Queue exhaustion or sink failure stops driving at a boundary with
exit code 5 and an explicit stderr error. Queued records are drained before
return. File records are flushed per line and synced when the worker finishes.
The socket mirror is nonblocking; disconnection/backpressure is an I/O error.
These are local CLI sinks, not a durable multi-subscriber broker. A blocked stdout
consumer can delay final worker shutdown; it cannot hold an SCF step. On I/O or
process failure the final log may remain nonterminal; status never claims liveness
from that log. No overhead/performance claim is made here.

## Agent-Controlled (Unix)

```
baandi run input.toml --agent --output-dir controlled
```

This mode uses P11B `Runtime`, with supervised approval and no autonomy policy.
It adds `OUTPUT/control/control.sock`, a nonblocking Unix datagram proposal socket,
and an exclusively created `OUTPUT/control-audit.jsonl`. The default socket
parent is created with mode 0700; the socket is 0600. `--control-socket PATH`
selects an existing private parent directory (no group/other permissions). Unix
socket path-length limits apply. Filesystem permissions define the local trust
boundary; this is not remote authentication, nor isolation from hostile processes
with the same user credentials. The endpoint is removed on normal completion.

A proposal is a JSON datagram, at most 4096 bytes, with no unknown fields:

```json
{"id":"proposal-1","calc_id":"ID_FROM_EVENT_LOG","observed_state_version":1,"observed_policy_revision":0,"max_iteration":20,"action":"mixing_beta","beta":0.2}
```

Actions are `pause`, `resume`, `cancel`, and `mixing_beta`; `beta` is only allowed
for the latter and must be finite in (0, 1]. Versions and expiry are mandatory;
IDs and the calculation identity are validated by the runtime. The initial policy
revision is 0; applied revisions appear in the audit. Late/stale requests are
rejected, never silently retargeted. A bound client socket receives a best-effort
admission response. `pending_approval` is not acceptance or application.

The launch terminal displays the immutable proposal. The operator types
`approve ID` or `deny ID` on that terminal. Only terminal stdin is an approval
source; piped stdin and proposal fields cannot grant approval. Approval expires
after 60 seconds and does not authorize stale revalidation. No terminal, absent
agent, or absent approval means continue the current numerical policy. A pause
requires approval; an approved pause waits for an approved resume or cancel.
Intake is limited to 16 datagrams and 16 terminal commands per boundary; runtime
pending/audit quotas remain in force. A run may finish before an agent replies.

The P11B audit sink writes and syncs each batch before committing a control
mutation, rolling back partial writes on reported failure. Audit rows include
request ID, action, observed versions, expiry iteration, requester, approver,
outcome, resulting state and policy revision. Debug strings for action/outcome
are display evidence, not a replay format. Cancellation has no fabricated final
numerical result. Checkpoint/control mode switches, general policy edits, remote
MCP access and durable replay are not advertised by this CLI adapter.

## Status and explain

```
baandi status observed
baandi status observed/events.jsonl
baandi status results/results.json
baandi status boundary.cp
baandi explain input.toml
baandi explain observed
baandi explain results/results.json
baandi explain evaluated.cp
```

Status prints structured evidence. Directories prefer events over final results
to avoid hiding an interrupted newer run behind an older result. The reader
streams bounded records, checks schema/calculation/sequence consistency, ignores
an incomplete final line with an explicit flag, and rejects malformed complete
lines. A nonterminal record means **process liveness unknown**, not “running”.
Final results report recorded convergence and energies.

Checkpoint inspection uses the merged P11C `checkpoint_io::read` API and the
[checkpoint v1 format](../checkpoint-format.md). Pass the file explicitly:
`.cp` routes damaged files through that reader, and `BAANDICP` magic recognizes
valid containers with any extension (including `.toml`). Directory lookup keeps
its established event-log/results priority and never guesses which checkpoint
belongs to a calculation. Schema, size limits, checksum and scalar/array decoding
are checked before displaying any evidence. No density or orbital arrays are
printed; the report gives dimensions/counts, model/execution identity and recorded
scalar history with units.

- `ExactContinuation`: saved iteration/state version/mode and latest recorded
  energy/residuals; convergence is unknown. Runtime envelopes additionally report
  saved runtime state, policy revision and audit count. Saved state is not proof
  that a process is running or paused now.
- `WarmStart`: an initial guess, with no fabricated iteration, convergence flag,
  history or energy. Missing evidence is JSON null.
- `ResultExport`: the recorded convergence flag (including false), final energy
  in Ha, charge in electrons, gap in Ha and latest iteration diagnostics. This is
  inspection evidence, not a restart capability.

Inspection neither restores a session nor executes SCF, and requires no original
input/UPF files. It does not establish restart compatibility, physical accuracy or
source authenticity. Kind/payload contradictions fail explicitly. Explain labels
these limits and reports that preflight was not rerun; only the original TOML
configuration route performs P8G preflight.

Explain runs the existing P8G preflight for TOML inputs without running SCF or
creating output files. It reports geometry, basis/grid, XC, resolved numerical
controls, electron/band counts, memory estimates, warnings and overrides. Invalid
fields retain the shared preflight's field name and allowed range. For run output
it reports recorded numbers and convergence evidence. It never interprets a large
residual or budget exhaustion as a physical diagnosis. Hypotheses are explicitly
labelled; none are automatically invented. Recorded engine diagnostics can be
extended when a stable serialized diagnostic contract is adopted.

## Natural language to specification (design only)

A future translator must emit a separate, inspectable TOML specification and a
list of unresolved fields. It must cite the provenance of user-supplied choices,
show units and explicit defaults, and run the same parser and preflight before
execution. Geometry, species/pseudopotentials, charge, XC, cutoffs, k sampling and
physical tolerances must never be silently guessed from ambiguous prose. Missing
physical intent requires clarification. Approval applies to the displayed,
validated spec, not to the original sentence. No translator is implemented here.

## Verification and limits

`cargo test -p dft-cli --test run execution_modes -- --nocapture` compares a
three-iteration synthetic small-cell fixture across the existing Classic entry,
the alias and Observed. `agent_socket_mode` extends this to AgentControlled
without intervention. It checks exact
density bytes and numerical JSON (excluding timings/timestamps), and unchanged
Classic console content (excluding measured wall times). This is adapter equivalence, not a QE comparison, a
multi-system performance result, or a claim about controlled trajectories.
Additional tests exercise stdout logs, partial status tails, malformed sequences
and field/range errors. A Unix socket-pair test exercises proposal intake, separate
operator approval, pause/resume, stale rejection and the durable audit. Path-bound
socket tests are ignored in the normal gate because the Codex sandbox rejects
Unix bind with EPERM; run them in a host environment that permits local sockets:

```
cargo test -p dft-cli --test run agent_socket_mode -- --ignored --nocapture
cargo test -p dft-cli --test run observed_socket -- --ignored --nocapture
```

These tests cover AgentControlled end-to-end equivalence, socket mirroring and
existing-log protection; they are not claimed as passing in this sandbox. The crate and workspace
gates are `scripts/check.sh dft-cli` and `scripts/check.sh`.

### Recorded verification (2026-11-04)

- `scripts/check.sh dft-cli`: exit 0; build, test, fmt and clippy passed.
- `scripts/check.sh`: exit 0; all four workspace gates passed after replacing the
  second Cargo binary with the compatibility launcher. This preserves the
  existing UPF test's unqualified `cargo run -p dft-cli` invocation.
- `cargo test -p dft-cli --test run execution_modes -- --nocapture`: exit 0;
  three SCF iterations, Classic/launcher/Observed density bytes and numerical JSON
  identical; energy difference 0 Ha. Agent socket-mode equivalence remains a
  host-only test and is not included in this result.
- `python3 -I scripts/make_readme_figures.py`: exit 0; existing generated README
  and figures unchanged. `python3 -I scripts/tests/test_readme_figures.py`: exit 0,
  three tests passed. No new QE comparison or physical validation was produced.
- Filesystem Unix datagram bind was denied with EPERM both under `/tmp` and this
  worktree. Socket-pair control tests passed; the two explicit path-binding tests
  remain ignored in the default gate, with host commands listed above.

Only ticket-owned files are changed. Session/dashboard JSON and publication are
orchestrator integration work; this checkout does not contain a running P11H
session record, and this task does not launch a nested delegate run. Adding this
adapter evidence to the generated README requires an evidence/script extension
outside this ticket's owned paths. After P11C integration, checkpoint inspection
uses the shared adapter; there is no substitute checkpoint format here.

### Checkpoint inspection follow-up

`cargo test -p dft-cli --test checkpoint_inspect -- --nocapture` exercises all
three checkpoint kinds, exact boundaries 0/1, converged and unconverged results,
and checksum/schema/truncation/kind-payload failures through both CLI commands.
Files remain byte-identical and inspection creates no outputs. Displayed scalar
energy agrees with the decoded engine record within 1e-14 Ha; the JSON display
is not a bitwise continuation format. Classic execution code is unchanged by this
follow-up; `execution_modes` retains its exact density/numerical comparison.
The main.rs merge resolution preserves both the checkpoint/restart commands and
the status/explain commands.

Follow-up verification: `scripts/check.sh dft-cli` and `scripts/check.sh` both
exited 0 (build/test/fmt/clippy). The two checkpoint inspection tests passed in
0.75 s; the Classic/launcher/Observed comparison again reported identical density
bytes and numerical JSON across three iterations, energy difference 0 Ha.
README generation and its three consistency tests exited 0 without additional
generated changes. No new numerical reference or performance claim is made.
Conflict markers were removed from main.rs; the Git index remains unmerged until
the reviewing orchestrator stages that resolution (this worktree's Git metadata
is read-only to the implementer). Incoming staged P11C changes are preserved.
