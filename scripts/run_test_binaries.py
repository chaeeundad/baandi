#!/usr/bin/env python3
"""Run already-built cargo test binaries concurrently (used by scripts/check.sh).

`cargo test` runs its test binaries one after another; on this workspace (180+
binaries, most under a second, a few near 30 s in debug) that serial tail is most
of the gate. This runner reads the JSON messages of
`cargo test --no-run --message-format=json` on stdin, runs every test executable
in its package directory (the working directory cargo uses) with a pool of
`--jobs` workers, longest first by the wall times of the previous run (cached in
the target directory; unknown binaries go first), and exits nonzero if any
binary fails. The failing binary's output tail is printed; everything else is
summarised in one line. Arguments after `--` are passed to every binary.

Coverage is identical to `cargo test --no-fail-fast` without doctests (run those
with `cargo test --doc`): the same binaries, the same libtest harness and
threads, the same ignored set. Only the order and overlap of binaries change.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--times", help="JSON cache of per-binary wall seconds")
    ap.add_argument("rest", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    extra = a.rest[1:] if a.rest[:1] == ["--"] else a.rest

    units = []
    for line in sys.stdin:
        try:
            m = json.loads(line)
        except ValueError:
            continue
        if m.get("reason") == "compiler-artifact" and m.get("executable") and m["profile"]["test"]:
            t = m["target"]
            cwd = os.path.dirname(m["manifest_path"])
            key = f"{os.path.basename(cwd)} {t['kind'][0]} {t['name']}"
            units.append((key, cwd, m["executable"]))
        elif m.get("reason") == "build-finished" and not m.get("success"):
            print("build failed")
            return 1
    if not units:
        print("no test binaries found")
        return 1

    known = {}
    if a.times:
        try:
            with open(a.times) as f:
                known = json.load(f)
        except (OSError, ValueError):
            known = {}
    units.sort(key=lambda u: -known.get(u[0], 1e9))

    def run(unit):
        key, cwd, exe = unit
        env = dict(os.environ, CARGO_MANIFEST_DIR=cwd)
        t0 = time.time()
        p = subprocess.run([exe, *extra], cwd=cwd, env=env, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, errors="replace")
        return key, p.returncode, time.time() - t0, p.stdout

    t0 = time.time()
    failed, totals, times = [], [0, 0, 0], {}
    with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        for key, rc, wall, out in pool.map(run, units):
            times[key] = round(wall, 2)
            for r in re.finditer(r"test result: \w+\. (\d+) passed; (\d+) failed; (\d+) ignored", out):
                for i in range(3):
                    totals[i] += int(r.group(i + 1))
            if rc != 0:
                failed.append(key)
                print(f"---- {key} (exit {rc}) ----")
                print("\n".join(out.splitlines()[-40:]))
    if a.times:
        try:
            with open(a.times, "w") as f:
                json.dump({**known, **times}, f, indent=0, sort_keys=True)
        except OSError:
            pass
    slow = sorted(times.items(), key=lambda kv: -kv[1])[:3]
    print(f"{len(units)} binaries, {totals[0]} passed, {totals[1]} failed, {totals[2]} ignored, "
          f"{time.time() - t0:.0f} s wall with {a.jobs} jobs; slowest: "
          + ", ".join(f"{k} {v:.0f} s" for k, v in slow))
    if failed:
        print("failed binaries: " + ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
