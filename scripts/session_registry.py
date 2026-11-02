#!/usr/bin/env python3
"""Registry of Codex delegation runs: artifacts/sessions.json (tracked in git).

Commands:
  start --task T --mode M --role R --ticket PATH --worktree PATH --branch B --wrapper-pid N   -> prints run_id
  set RUN_ID key=value ...      (ints for *pid / exit_code are converted; "null" clears)
  get RUN_ID key                prints one field
  list [--status S]

Fields: run_id, task_id, mode, role, status, codex_pid, wrapper_pid, session_id, branch, worktree,
ticket, report, started_at, ended_at, exit_code, note.
status: running | done | blocked | failed | killed | unknown
"""
import argparse, datetime, fcntl, json, os, sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "artifacts", "sessions.json")
INT_KEYS = {"codex_pid", "wrapper_pid", "exit_code"}


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Registry:
    def __enter__(self):
        os.makedirs(os.path.dirname(PATH), exist_ok=True)
        self.lock = open(PATH + ".lock", "w")
        fcntl.flock(self.lock, fcntl.LOCK_EX)
        try:
            self.runs = json.load(open(PATH))
        except (FileNotFoundError, json.JSONDecodeError):
            self.runs = []
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            tmp = PATH + ".tmp"
            with open(tmp, "w") as f:
                json.dump(self.runs, f, indent=2, ensure_ascii=False)
                f.write("\n")
            os.replace(tmp, PATH)
        fcntl.flock(self.lock, fcntl.LOCK_UN)
        self.lock.close()

    def find(self, run_id):
        for r in self.runs:
            if r["run_id"] == run_id:
                return r
        sys.exit(f"unknown run_id {run_id}")


def coerce(k, v):
    if v == "null":
        return None
    return int(v) if k in INT_KEYS else v


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("start")
    for a in ("task", "mode", "role", "ticket", "worktree", "branch"):
        s.add_argument("--" + a, default="")
    s.add_argument("--wrapper-pid", type=int, default=None)
    t = sub.add_parser("set"); t.add_argument("run_id"); t.add_argument("kv", nargs="+")
    g = sub.add_parser("get"); g.add_argument("run_id"); g.add_argument("key")
    l = sub.add_parser("list"); l.add_argument("--status")
    a = p.parse_args()

    with Registry() as reg:
        if a.cmd == "start":
            n = sum(1 for r in reg.runs if r["task_id"] == a.task) + 1
            run = {"run_id": f"{a.task}#{n}", "task_id": a.task, "mode": a.mode, "role": a.role,
                   "status": "running", "codex_pid": None, "wrapper_pid": a.wrapper_pid,
                   "session_id": None, "branch": a.branch, "worktree": a.worktree, "ticket": a.ticket,
                   "report": None, "started_at": now(), "ended_at": None, "exit_code": None, "note": ""}
            reg.runs.append(run)
            print(run["run_id"])
        elif a.cmd == "set":
            r = reg.find(a.run_id)
            for kv in a.kv:
                k, v = kv.split("=", 1)
                r[k] = coerce(k, v)
            if r["status"] != "running" and not r.get("ended_at"):
                r["ended_at"] = now()
        elif a.cmd == "get":
            v = reg.find(a.run_id).get(a.key)
            print("" if v is None else v)
        elif a.cmd == "list":
            for r in reg.runs:
                if a.status and r["status"] != a.status:
                    continue
                print(f'{r["run_id"]:22} {r["status"]:8} pid={r["codex_pid"]} role={r["role"]} session={r["session_id"]}')


if __name__ == "__main__":
    main()
