#!/usr/bin/env python3 -I
"""Remove the build directories (`target/`) of finished worktrees to free disk space.

Two worktree families are scanned:

1. Codex task worktrees `../rust-dft-wt/<TASK>` (relative to the main checkout):
   `<TASK>/target` is removed only when the latest run of <TASK> recorded in
   artifacts/sessions.json is not "running".
2. Claude agent worktrees `<main checkout>/.claude/worktrees/<name>`: `<name>/target`
   is removed only when the worktree is finished: its HEAD is merged into `main`
   (an ancestor of main), its working tree is clean (no uncommitted or untracked
   files) and no running process has the worktree path on its command line.
   `--include-unmerged` also takes clean, idle worktrees whose HEAD is not merged
   (rejected or abandoned branches; their commits stay in git).

`--check-dirs` removes only the gate's side directories `target/check-nested` (the
nested `cargo run -p dft-cli` build, unused since check.sh passes DFT_CLI_BIN) and
`target/check-clippy` (clippy's own build), in every worktree and the main checkout,
whatever the task status, unless a process is using that worktree. They are rebuilt
by the next gate.

Every action prints its size; `--dry-run` lists without removing. Sources,
uncommitted changes and everything else stay. Standard library only.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main_checkout():
    """The main checkout (the common git dir's parent), also when run from a worktree."""
    p = subprocess.run(["git", "-C", HERE, "rev-parse", "--path-format=absolute", "--git-common-dir"],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    if p.returncode == 0 and p.stdout.strip():
        return os.path.dirname(p.stdout.strip())
    return HERE


ROOT = main_checkout()
WT = os.path.realpath(os.path.join(ROOT, "..", "rust-dft-wt"))
CLAUDE_WT = os.path.join(ROOT, ".claude", "worktrees")


def size(path):
    total = 0
    for base, _, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(base, f)).st_blocks * 512
            except OSError:
                pass
    return total


def busy_paths():
    """Command lines and working directories of running processes (a worktree named in
    one, or containing a process's working directory, is in use)."""
    p = subprocess.run(["ps", "-A", "-o", "command="], stdout=subprocess.PIPE, text=True)
    procs = [line for line in p.stdout.splitlines() if "wt_cleanup" not in line]
    q = subprocess.run(["lsof", "-w", "-a", "-d", "cwd", "-Fn"], stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL, text=True)
    mine, pid = {os.getpid(), os.getppid()}, None
    for line in q.stdout.splitlines():
        if line.startswith("p"):
            pid = int(line[1:])
        elif line.startswith("n") and pid not in mine:
            procs.append(line[1:] + "/")
    return procs


def in_use(path, procs):
    real = os.path.realpath(path)
    nested = (CLAUDE_WT + "/", os.path.realpath(CLAUDE_WT) + "/")
    hits = (c for c in procs if path + "/" in c or real + "/" in c or c.endswith(path) or c.endswith(real))
    if real == os.path.realpath(ROOT):  # agent worktrees live inside the main checkout
        hits = (c for c in hits if not any(n in c for n in nested))
    return any(True for _ in hits)


def git(wt, *args):
    return subprocess.run(["git", "-C", wt, *args], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, text=True)


def remove(label, path, dry, freed):
    n = size(path)
    freed[0] += n
    print(f"{'would remove' if dry else 'remove'} {label} {n / 1e9:.2f} GB")
    if not dry:
        shutil.rmtree(path, ignore_errors=True)


def codex_targets(a, procs, freed):
    if not os.path.isdir(WT):
        return
    try:
        runs = json.load(open(os.path.join(ROOT, "artifacts", "sessions.json")))
    except (OSError, ValueError):
        runs = []
    latest = {}
    for r in runs:
        latest[r["task_id"]] = r["status"]
    for task in sorted(os.listdir(WT)):
        target = os.path.realpath(os.path.join(WT, task, "target"))
        if not os.path.isdir(target) or os.path.basename(target) != "target" or os.path.dirname(os.path.dirname(target)) != WT:
            continue  # only <WT>/<task>/target, nothing else
        status = latest.get(task)
        if status is None or status == "running" or in_use(os.path.join(WT, task), procs):
            print(f"keep   {task} (status {status})")
            continue
        remove(f"{task}/target (status {status})", target, a.dry_run, freed)


def claude_targets(a, procs, freed):
    if not os.path.isdir(CLAUDE_WT):
        return
    for name in sorted(os.listdir(CLAUDE_WT)):
        wt = os.path.join(CLAUDE_WT, name)
        target = os.path.join(wt, "target")
        if not os.path.isdir(target) or os.path.islink(target):
            continue
        if os.path.realpath(wt) == os.path.realpath(HERE):
            print(f"keep   .claude/worktrees/{name} (this checkout)")
            continue
        if in_use(wt, procs):
            print(f"keep   .claude/worktrees/{name} (a running process uses it)")
            continue
        status = git(wt, "status", "--porcelain")
        if status.returncode != 0 or status.stdout.strip():
            print(f"keep   .claude/worktrees/{name} (uncommitted changes or not a git worktree)")
            continue
        merged = git(wt, "merge-base", "--is-ancestor", "HEAD", "main").returncode == 0
        if not merged and not a.include_unmerged:
            print(f"keep   .claude/worktrees/{name} (HEAD not merged into main; --include-unmerged takes it)")
            continue
        remove(f".claude/worktrees/{name}/target ({'merged' if merged else 'unmerged, clean'})",
               target, a.dry_run, freed)


def check_dirs(a, procs, freed):
    roots = [ROOT]
    for base in (WT, CLAUDE_WT):
        if os.path.isdir(base):
            roots += [os.path.join(base, n) for n in sorted(os.listdir(base))]
    for wt in roots:
        if in_use(wt, procs):
            print(f"keep   {os.path.relpath(wt, os.path.dirname(ROOT))} check dirs (a running process uses it)")
            continue
        for sub in ("check-nested", "check-clippy"):
            d = os.path.join(wt, "target", sub)
            if os.path.isdir(d) and not os.path.islink(d):
                remove(os.path.relpath(d, os.path.dirname(ROOT)), d, a.dry_run, freed)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check-dirs", action="store_true",
                    help="remove only target/check-nested and target/check-clippy (all worktrees)")
    ap.add_argument("--include-unmerged", action="store_true",
                    help="also remove targets of clean, idle Claude worktrees whose HEAD is not in main")
    a = ap.parse_args(argv)
    procs = busy_paths()
    freed = [0]
    if a.check_dirs:
        check_dirs(a, procs, freed)
    else:
        codex_targets(a, procs, freed)
        claude_targets(a, procs, freed)
    print(f"{'would free' if a.dry_run else 'freed'} {freed[0] / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
