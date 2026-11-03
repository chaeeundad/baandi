#!/usr/bin/env python3
"""Build artifacts/agents.json (dashboard Agents card) from artifacts/sessions.json plus fixed non-Codex agents."""
import json, os
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "artifacts")
runs = json.load(open(os.path.join(root, "sessions.json")))
state = {"running": "active", "done": "done", "blocked": "blocked", "failed": "failed", "killed": "killed", "unknown": "unknown"}
rows = [
    {"id": "A1", "agent": "Claude Code", "role": "Architect, orchestrator, reviewer", "state": "active",
     "task": "총괄 설계와 검토: 성능 개선 지휘(목표: QE 같은 코어 수에서 1.5배 이내), 새 기능 조사 결과를 마스터 로드맵으로 통합", "ref": "this session", "updated": "2026-11-02"},
    {"id": "A2", "agent": "dft5090", "role": "Remote QE reference runner (via dev08)", "state": "idle",
     "task": "dev08에서 QE 기준 계산과 성능 측정을 대신 실행하는 세션(요청이 오면 실행, 지금은 대기)", "ref": "Remote Control session", "updated": "2026-11-02"},
]
titles = json.load(open(os.path.join(root, "task_titles.json")))
# one row per task (latest run), running tasks first, Korean description in the task column
latest = {}
for r in runs:
    latest[r["task_id"]] = r
order = sorted(latest.values(), key=lambda r: (r["status"] != "running", r["started_at"]))
for i, r in enumerate(order, start=3):
    pid = f"pid {r['codex_pid']}" if r.get("codex_pid") else "pid not recorded"
    sid = (r.get("session_id") or "no session")[:8]
    note = (r.get("note") or "").split(";")[0]
    rows.append({"id": f"A{i}", "agent": "Codex", "role": r["role"], "state": state.get(r["status"], r["status"]),
                 "task": f"{r['task_id']}: {titles.get(r['task_id'], note or r['task_id'])}",
                 "ref": f"{pid}, session {sid}, {r['branch']}", "updated": (r.get("ended_at") or r["started_at"])[:10]})
json.dump(rows, open(os.path.join(root, "agents.json"), "w"), indent=2, ensure_ascii=False)
open(os.path.join(root, "agents.json"), "a").write("\n")
print(len(rows), "agent rows")
