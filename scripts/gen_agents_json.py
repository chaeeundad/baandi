#!/usr/bin/env python3
"""Build artifacts/agents.json (dashboard Agents card) from artifacts/sessions.json plus fixed non-Codex agents."""
import json, os
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "artifacts")
runs = json.load(open(os.path.join(root, "sessions.json")))
state = {"running": "active", "done": "done", "blocked": "blocked", "failed": "failed", "killed": "killed", "unknown": "unknown"}
rows = [
    {"id": "A1", "agent": "Claude Code", "role": "Architect, orchestrator, reviewer", "state": "active",
     "task": "Library survey (ADR 0003), AI-first repo structure; next: backend spike P0C", "ref": "this session", "updated": "2026-11-02"},
    {"id": "A2", "agent": "dft5090", "role": "Remote QE reference runner (via dev08)", "state": "idle",
     "task": "Readiness check done: pw.x 7.6 on dev08 reachable (reported). Waiting for P6 inputs with SHA-256", "ref": "Remote Control session", "updated": "2026-11-02"},
]
for i, r in enumerate(runs, start=3):
    pid = f"pid {r['codex_pid']}" if r.get("codex_pid") else "pid not recorded"
    sid = (r.get("session_id") or "no session")[:8]
    note = (r.get("note") or "").split(";")[0]
    rows.append({"id": f"A{i}", "agent": "Codex", "role": r["role"], "state": state.get(r["status"], r["status"]),
                 "task": f"{r['task_id']}: {note}" if note else r["task_id"],
                 "ref": f"{pid}, session {sid}, {r['branch']}", "updated": (r.get("ended_at") or r["started_at"])[:10]})
json.dump(rows, open(os.path.join(root, "agents.json"), "w"), indent=2, ensure_ascii=False)
open(os.path.join(root, "agents.json"), "a").write("\n")
print(len(rows), "agent rows")
