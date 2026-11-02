#!/usr/bin/env bash
# Claude Code -> Codex 위임 래퍼 (codex exec 기반).
# 사용법:
#   scripts/delegate.sh new    <TASK_ID> <ticket.md>   # worktree+브랜치 생성 후 새 세션
#   scripts/delegate.sh resume <TASK_ID> <message.md>  # 같은 세션 이어서 수정 요청
#   scripts/delegate.sh fork   <TASK_ID> <ticket.md> <BASE_SESSION_ID>  # 프라임된 기준 세션에서 분기
# 산출물: .delegate/<TASK_ID>/{events.jsonl,report.json,session_id}
set -euo pipefail

mode="${1:?mode}"; task="${2:?TASK_ID}"; prompt_file="${3:?prompt file}"; base="${4:-}"
repo="$(git rev-parse --show-toplevel)"
wt="$repo/../rust-dft-wt/$task"
out="$repo/.delegate/$task"
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
mkdir -p "$out"

if [[ "$mode" == "new" || "$mode" == "fork" ]]; then
  [[ -d "$wt" ]] || git -C "$repo" worktree add -b "task/$task" "$wt" HEAD
fi

common=(-C "$wt" -s workspace-write --json
        --output-schema "$repo/docs/schemas/report.schema.json"
        -o "$out/report.json")

case "$mode" in
  new)    codex exec "${common[@]}" - < "$prompt_file" > "$out/events.jsonl" ;;
  fork)   codex exec fork "$base" "${common[@]}" - < "$prompt_file" > "$out/events.jsonl" ;;
  resume) codex exec resume "$(cat "$out/session_id")" "${common[@]}" - < "$prompt_file" >> "$out/events.jsonl" ;;
  *) echo "unknown mode $mode" >&2; exit 2 ;;
esac

# 세션 ID 기록 (이벤트 스트림에서 첫 thread/session id 추출; 스키마는 실제 출력으로 확인 필요)
if [[ ! -f "$out/session_id" ]]; then
  grep -o -E '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}' "$out/events.jsonl" | head -1 > "$out/session_id" || true
fi
echo "report: $out/report.json  session: $(cat "$out/session_id" 2>/dev/null)  worktree: $wt"
