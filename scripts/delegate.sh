#!/usr/bin/env bash
# Claude Code -> Codex delegation wrapper (codex exec). Every run is recorded in artifacts/sessions.json
# (PID, role, session id, status, exit code). Never kill Codex by name; use the `stop` mode.
#
#   ROLE=implementer scripts/delegate.sh new    <TASK_ID> <ticket.md>
#   scripts/delegate.sh resume <TASK_ID> <message.md>
#   scripts/delegate.sh fork   <TASK_ID> <ticket.md> <BASE_SESSION_ID>
#   scripts/delegate.sh stop   <TASK_ID>          # kills only the recorded Codex PID of the latest running run
# Outputs: .delegate/<TASK_ID>/{events.jsonl,report.json,session_id}
# Launch this script itself with the harness' run_in_background, not with a trailing '&'.
set -uo pipefail

mode="${1:?mode}"; task="${2:?TASK_ID}"; prompt_file="${3:-}"; base="${4:-}"
role="${ROLE:-implementer}"
repo="$(git rev-parse --show-toplevel)"
wt="$repo/../rust-dft-wt/$task"
out="$repo/.delegate/$task"
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
reg() { python3 "$repo/scripts/session_registry.py" "$@"; }
mkdir -p "$out"

if [[ "$mode" == "stop" ]]; then
  run="$(python3 - "$task" <<'PY'
import json,sys
runs=[r for r in json.load(open("artifacts/sessions.json")) if r["task_id"]==sys.argv[1] and r["status"]=="running"]
print(runs[-1]["run_id"] if runs else "")
PY
)"
  [[ -n "$run" ]] || { echo "no running run for $task"; exit 1; }
  pid="$(reg get "$run" codex_pid)"
  [[ -n "$pid" ]] || { echo "no codex_pid recorded for $run"; exit 1; }
  if ps -p "$pid" -o command= | grep -q "codex"; then
    reg set "$run" status=killed note="stopped by Claude via delegate.sh stop"
    kill "$pid" && echo "stopped $run (pid $pid)"
  else
    echo "pid $pid is not a codex process; refusing to kill"; exit 1
  fi
  exit 0
fi

[[ -n "$prompt_file" ]] || { echo "prompt file required" >&2; exit 2; }
if [[ "$mode" == "new" || "$mode" == "fork" ]]; then
  [[ -d "$wt" ]] || git -C "$repo" worktree add -b "task/$task" "$wt" HEAD
fi

common=(-C "$wt" -s workspace-write --json
        --output-schema "$repo/docs/schemas/report.schema.json"
        -o "$out/report.json")
case "$mode" in
  new)    cmd=(codex exec "${common[@]}" -) ;;
  fork)   cmd=(codex exec fork "$base" "${common[@]}" -) ;;
  resume) cmd=(codex exec resume "$(cat "$out/session_id")" "${common[@]}" -) ;;
  *) echo "unknown mode $mode" >&2; exit 2 ;;
esac

run="$(reg start --task "$task" --mode "$mode" --role "$role" --ticket "$prompt_file" \
       --worktree "$wt" --branch "task/$task" --wrapper-pid $$)"
[[ "$mode" == "resume" ]] && redir=">>" || redir=">"
if [[ "$mode" == "resume" ]]; then "${cmd[@]}" < "$prompt_file" >> "$out/events.jsonl" & else "${cmd[@]}" < "$prompt_file" > "$out/events.jsonl" & fi
cpid=$!
reg set "$run" codex_pid="$cpid"
wait "$cpid"; rc=$?

sid="$(grep -o -E '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}' "$out/events.jsonl" | head -1)"
[[ -n "$sid" && ! -f "$out/session_id" ]] && echo "$sid" > "$out/session_id"
status="$(python3 - "$out/report.json" "$rc" <<'PY'
import json,sys
try: s=json.load(open(sys.argv[1]))["status"].lower()
except Exception: s="failed" if sys.argv[2]!="0" else "unknown"
print(s)
PY
)"
[[ "$(reg get "$run" status)" == "killed" ]] && status=killed
reg set "$run" status="$status" exit_code="$rc" session_id="${sid:-null}" report="$out/report.json"
echo "run=$run status=$status exit=$rc pid=$cpid session=${sid:-?} worktree=$wt"
