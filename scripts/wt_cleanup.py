#!/usr/bin/env python3 -I
"""Remove the build directories (`target/`) of finished task worktrees to free disk space.

Reads artifacts/sessions.json: a worktree ../rust-dft-wt/<TASK>/target is removed only when the latest
recorded run of <TASK> is not "running". Sources, uncommitted changes and everything else stay.
`--dry-run` lists what would be removed. Standard library only.
"""
import argparse, json, os, shutil, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WT = os.path.realpath(os.path.join(ROOT, "..", "rust-dft-wt"))


def size(path):
    total = 0
    for base, _, files in os.walk(path):
        for f in files:
            try:
                total += os.lstat(os.path.join(base, f)).st_size
            except OSError:
                pass
    return total


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    runs = json.load(open(os.path.join(ROOT, "artifacts", "sessions.json")))
    latest = {}
    for r in runs:
        latest[r["task_id"]] = r["status"]
    freed = 0
    for task in sorted(os.listdir(WT)):
        target = os.path.realpath(os.path.join(WT, task, "target"))
        if not os.path.isdir(target) or os.path.basename(target) != "target" or os.path.dirname(os.path.dirname(target)) != WT:
            continue  # only <WT>/<task>/target, nothing else
        status = latest.get(task)
        if status is None or status == "running":
            print(f"keep   {task} (status {status})")
            continue
        n = size(target)
        freed += n
        print(f"{'would remove' if a.dry_run else 'remove'} {task}/target {n / 1e9:.2f} GB (status {status})")
        if not a.dry_run:
            shutil.rmtree(target, ignore_errors=True)
    print(f"{'would free' if a.dry_run else 'freed'} {freed / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
