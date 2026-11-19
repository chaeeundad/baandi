#!/usr/bin/env python3
"""Select the workspace crates that the `fast` tier of scripts/check.sh must check.

Prints one package name per line: the crates touched by the change and their direct
dependents (normal, dev and build dependencies within the workspace), or the single
word `WORKSPACE` when the change touches the whole workspace. Standard library only.

The change = `git diff --name-only <base>` (committed and uncommitted changes against
the base) plus untracked files. The base is `CHECK_BASE` when set, else the nearer of
`git merge-base HEAD main` and `git merge-base HEAD origin/main` (agent branches fork
from a local main that may be ahead of origin/main; taking the nearer base keeps the
diff to the agent's own work).

Mapping of a changed path to crates:
- `crates/<dir>/...`: the package whose manifest lives in `crates/<dir>`.
- `Cargo.toml`, `Cargo.lock`, `rust-toolchain.toml`, `.cargo/...`, `clippy.toml`,
  `rustfmt.toml`: WORKSPACE (every crate is rebuilt).
- any other path (docs, scripts, tests/reference, AGENTS.md, ...): the crates whose
  sources mention the file's basename (`git grep -F`), because some crate tests read
  root files (AGENTS.md, ownership.toml, docs, scripts, tests/reference). Only `.rs`
  files are searched. Conservative: a basename match may add a crate that does not
  read the file.

`--explain` prints the base, the changed paths and the reason per crate on stderr.
"""

import argparse
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE_FILES = {"Cargo.toml", "Cargo.lock", "rust-toolchain.toml", "clippy.toml",
                   ".clippy.toml", "rustfmt.toml", ".rustfmt.toml"}


def git(*args, check=True):
    p = subprocess.run(["git", "-C", ROOT, *args], stdout=subprocess.PIPE,
                       stderr=subprocess.DEVNULL, text=True)
    if check and p.returncode != 0:
        return None
    return p.stdout


def is_ancestor(a, b):
    return subprocess.run(["git", "-C", ROOT, "merge-base", "--is-ancestor", a, b],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def pick_base():
    env = os.environ.get("CHECK_BASE")
    if env:
        return env
    bases = []
    for ref in ("main", "origin/main"):
        out = git("merge-base", "HEAD", ref)
        if out and out.strip():
            bases.append(out.strip())
    if not bases:
        return None
    best = bases[0]
    for b in bases[1:]:
        if is_ancestor(best, b):  # b is nearer to HEAD
            best = b
    return best


def metadata():
    cargo = os.environ.get("CARGO", "cargo")
    p = subprocess.run([cargo, "metadata", "--format-version", "1", "--no-deps", "--offline"],
                       cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    if p.returncode != 0:
        return None
    return json.loads(p.stdout)


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--explain", action="store_true")
    ap.add_argument("--no-dependents", action="store_true", help="changed crates only")
    a = ap.parse_args(argv)
    say = (lambda *x: print(*x, file=sys.stderr)) if a.explain else (lambda *x: None)

    base = pick_base()
    if base is None:
        say("no base found (no main / origin/main): WORKSPACE")
        print("WORKSPACE")
        return 0
    diff = git("diff", "--name-only", base)
    untracked = git("ls-files", "--others", "--exclude-standard")
    if diff is None:
        say(f"git diff against {base} failed: WORKSPACE")
        print("WORKSPACE")
        return 0
    paths = sorted({p for p in (diff + (untracked or "")).splitlines() if p.strip()})
    say(f"base {base[:12]}; {len(paths)} changed paths")

    m = metadata()
    if m is None:
        say("cargo metadata failed: WORKSPACE")
        print("WORKSPACE")
        return 0
    by_dir, names = {}, set()
    for p in m["packages"]:
        rel = os.path.relpath(os.path.dirname(p["manifest_path"]), m["workspace_root"])
        by_dir[rel] = p["name"]
        names.add(p["name"])

    changed, why = set(), {}
    for path in paths:
        if path in WORKSPACE_FILES or path.startswith(".cargo/"):
            say(f"  {path}: workspace-wide file")
            print("WORKSPACE")
            return 0
        parts = path.split("/")
        if parts[0] == "crates" and len(parts) > 2:
            name = by_dir.get("/".join(parts[:2]))
            if name:
                changed.add(name)
                why.setdefault(name, f"changed ({path})")
                continue
        base_name = os.path.basename(path)
        hits = git("grep", "-l", "-F", "-e", base_name, "--", ":(glob)crates/**/*.rs",
                   check=False) or ""
        for hit in hits.splitlines():
            hp = hit.split("/")
            name = by_dir.get("/".join(hp[:2])) if len(hp) > 2 else None
            if name:
                changed.add(name)
                why.setdefault(name, f"reads {path} ({hit})")

    selected = set(changed)
    if not a.no_dependents:
        for p in m["packages"]:
            deps = {d["name"] for d in p["dependencies"] if d["name"] in names}
            hit = deps & changed
            if hit and p["name"] not in selected:
                selected.add(p["name"])
                why[p["name"]] = "direct dependent of " + ", ".join(sorted(hit))
    for name in sorted(selected):
        say(f"  {name}: {why[name]}")
        print(name)
    if not selected:
        say("  no crate affected")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
