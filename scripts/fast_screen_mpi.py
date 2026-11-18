#!/usr/bin/env python3
"""Stage-1 MPI screen (docs/performance/fast-loop-protocol.md, section "Harness (MPI)").

Runs the M2 profile harness (crates/dft-engine/tests/m2_entry_cost.rs, simulated
ranks = threads of one process) on ONE binary with an environment A/B, N
alternating pairs, and prints the median wall per arm, the spread, the
deterministic counters (collective calls and bytes, parallel.reductions per site,
H columns, Davidson iterations, rotation FLOPs, reduced-solve sum of m^3) and a
verdict REJECT / NEUTRAL / CANDIDATE / UNUSABLE.

Mac wall on simulated ranks is an OBSERVATION (label O), never adoption
evidence; the counters are exact. Mandatory self-check: a counter that is zero on
both arms is reported as "unusable (all zero)" and is never read as "no
difference".

Example (python3 -I compatible, no network):
  python3 -I scripts/fast_screen_mpi.py --upf /path/Si.pz-vbc.UPF \
      --b BAANDI_M2_PANEL_GS=fused --expect collectives.tree_sum.calls:down
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLASS_LIMIT_HA = {"strict": 1e-10, "standard": 5e-9}
TIMING_KEYS = ("wall_s",)


def flatten(record: dict) -> tuple[dict[str, float], dict[str, float]]:
    """(deterministic counters, timings) of one M2SCREEN_JSON record."""
    counters: dict[str, float] = {}
    timings: dict[str, float] = {"wall_s": float(record["wall_s"])}
    for key in ("scf", "davidson", "h_columns", "rotation_flops", "reduced_solves", "reduced_sum_m3"):
        if key in record:
            counters[key] = float(record[key])
    for site, (calls, payload) in record.get("reductions", {}).items():
        counters[f"reductions.{site}.calls"] = float(calls)
        counters[f"reductions.{site}.doubles"] = float(payload)
    for op, values in record.get("collectives", {}).items():
        counters[f"collectives.{op}.calls"] = float(values[0])
        counters[f"collectives.{op}.bytes"] = float(values[1])
        timings[f"collectives.{op}.s"] = float(values[2])
    for name, value in record.get("m2_counters", {}).items():
        counters[f"m2.{name}"] = float(value)
    for name, value in record.get("grid_counters", {}).items():
        counters[f"grid.{name}"] = float(value)
    # M2-OWNED: counts are deterministic counters, *_s entries are timings.
    for name, value in record.get("owned_counters", {}).items():
        if name.endswith("_s"):
            timings[f"owned.{name}"] = float(value)
        else:
            counters[f"owned.{name}"] = float(value)
    for name, value in record.get("phases", {}).items():
        timings[f"phase.{name}"] = float(value)
    for name, value in record.get("stages", {}).items():
        timings[f"stage.{name}"] = float(value)
    return counters, timings


def find_binary(build: bool) -> Path:
    """The release test binary of m2_entry_cost (built once; both arms share it)."""
    cmd = [
        "cargo", "test", "--release", "-p", "dft-engine", "--test", "m2_entry_cost",
        "--no-run", "--message-format=json",
    ]
    if not build:
        cmd.append("--offline")
    out = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"build failed (exit {out.returncode}):\n{out.stderr[-4000:]}")
    exe = None
    for line in out.stdout.splitlines():
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue
        if msg.get("reason") == "compiler-artifact" and msg.get("target", {}).get("name") == "m2_entry_cost":
            exe = msg.get("executable") or exe
    if not exe:
        sys.exit("could not locate the m2_entry_cost test executable")
    return Path(exe)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_env(items: list[str]) -> dict[str, str]:
    env = {}
    for item in items:
        for part in item.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                sys.exit(f"bad environment item {part!r} (want NAME=value)")
            k, v = part.split("=", 1)
            env[k] = v
    return env


def run_once(binary: Path, base_env: dict, arm_env: dict, timeout: float) -> dict:
    env = dict(os.environ)
    env.update(base_env)
    env.update(arm_env)
    load = os.getloadavg()[0]
    t0 = time.time()
    out = subprocess.run(
        [str(binary), "--ignored", "--nocapture", "--exact", "m2_entry_cost_profile"],
        env=env, capture_output=True, text=True, timeout=timeout,
    )
    elapsed = time.time() - t0
    if out.returncode != 0:
        raise RuntimeError(f"run failed (exit {out.returncode}):\n{out.stdout[-3000:]}\n{out.stderr[-3000:]}")
    records = [json.loads(l.split(" ", 1)[1]) for l in out.stdout.splitlines() if l.startswith("M2SCREEN_JSON ")]
    if len(records) != 1:
        raise RuntimeError(f"expected one M2SCREEN_JSON line, got {len(records)} (INVALID run)")
    rec = records[0]
    if rec["wall_s"] < 0.1:
        raise RuntimeError("wall below 0.1 s: INVALID run")
    rec["_load"] = load
    rec["_process_s"] = elapsed
    return rec


def spread(values: list[float]) -> float:
    return max(values) - min(values) if values else 0.0


def verdict(arms: dict, expect: list[tuple[str, str]], class_limit: float, min_drop: float) -> dict:
    """Pure decision logic (unit-tested in scripts/tests/test_fast_screen_mpi.py).

    arms: {"A": [record, ...], "B": [record, ...]} of M2SCREEN_JSON records.
    """
    res: dict = {"notes": []}
    flat = {k: [flatten(r) for r in v] for k, v in arms.items()}
    # Counters must be identical within an arm (they are trajectory properties).
    for arm, rows in flat.items():
        keys = rows[0][0].keys()
        bad = [k for k in keys if len({r[0].get(k) for r in rows}) > 1]
        if bad:
            res["notes"].append(f"arm {arm}: counters vary between repeats: {', '.join(sorted(bad)[:6])}")
    ca, cb = flat["A"][0][0], flat["B"][0][0]
    counters = {}
    unusable = []
    for k in sorted(set(ca) | set(cb)):
        a, b = ca.get(k, 0.0), cb.get(k, 0.0)
        if a == 0.0 and b == 0.0:
            unusable.append(k)
            continue
        counters[k] = (a, b, (b - a) / a if a else float("inf"))
    res["counters"] = counters
    res["unusable"] = unusable
    wa = [r[1]["wall_s"] for r in flat["A"]]
    wb = [r[1]["wall_s"] for r in flat["B"]]
    ma, mb = statistics.median(wa), statistics.median(wb)
    sp = max(spread(wa), spread(wb))
    res["wall"] = {"A": ma, "B": mb, "delta": (mb - ma) / ma, "spread_s": sp, "A_runs": wa, "B_runs": wb}
    # Correctness: energy class against arm A, bit identity reported.
    ea = arms["A"][0]["energy_ha"]
    de = max(abs(r["energy_ha"] - ea) for r in arms["B"])
    same_bits = all(r["energy_bits"] == arms["A"][0]["energy_bits"] for r in arms["B"])
    same_density = all(r["density_fnv"] == arms["A"][0]["density_fnv"] for r in arms["B"])
    reproducible = len({r["density_fnv"] for r in arms["B"]}) == 1 and len({r["density_fnv"] for r in arms["A"]}) == 1
    res["correctness"] = {
        "dE_ha": de, "class_limit_ha": class_limit, "energy_bits_equal": same_bits,
        "density_hash_equal": same_density, "reproducible_per_arm": reproducible,
        "identity": "bit-identical" if (same_bits and same_density) else "rounding class",
    }
    reasons = []
    if de > class_limit:
        reasons.append(f"|dE| {de:.3e} Ha above the class limit {class_limit:.0e}")
    if not reproducible:
        reasons.append("density hash differs between repeats of one arm (determinism contract)")
    moved_ok = []
    expect_unusable = []
    for key, direction in expect:
        if key in unusable or key not in counters:
            expect_unusable.append(key)
            continue
        a, b, rel = counters[key]
        good = (b < a) if direction == "down" else (b > a)
        if not good:
            reasons.append(f"expected counter {key} did not move {direction} ({a:.6g} -> {b:.6g})")
        elif abs(rel) >= min_drop:
            moved_ok.append(key)
    if reasons:
        res["verdict"] = "REJECT"
    elif expect and len(expect_unusable) == len(expect):
        res["verdict"] = "UNUSABLE"
        reasons.append("every expected counter is unusable (all zero on both arms)")
    elif mb > ma + sp:
        res["verdict"] = "REJECT"
        reasons.append(f"median wall worse beyond the spread ({ma:.3f} -> {mb:.3f} s, spread {sp:.3f} s)")
    elif mb < ma - sp:
        res["verdict"] = "CANDIDATE"
        reasons.append(f"median wall better beyond the spread ({ma:.3f} -> {mb:.3f} s, spread {sp:.3f} s)")
    elif moved_ok:
        res["verdict"] = "CANDIDATE"
        reasons.append(f"large counter drop ({', '.join(moved_ok)}), wall within the spread")
    else:
        res["verdict"] = "NEUTRAL"
        reasons.append("no counter moved beyond the minimum and the wall is within the spread")
    res["reasons"] = reasons
    res["expect_unusable"] = expect_unusable
    return res


def report(res: dict, show_all: bool) -> str:
    lines = []
    w = res["wall"]
    lines.append(
        f"wall (O, simulated ranks): A median {w['A']:.3f} s {['%.3f' % x for x in w['A_runs']]}, "
        f"B median {w['B']:.3f} s {['%.3f' % x for x in w['B_runs']]}, delta {100 * w['delta']:+.2f} %, spread {w['spread_s']:.3f} s"
    )
    c = res["correctness"]
    lines.append(
        f"correctness: {c['identity']}, |dE| {c['dE_ha']:.3e} Ha (limit {c['class_limit_ha']:.0e}), "
        f"energy bits equal {c['energy_bits_equal']}, density hash equal {c['density_hash_equal']}, reproducible {c['reproducible_per_arm']}"
    )
    lines.append("counters (per rank, exact):")
    for k, (a, b, rel) in res["counters"].items():
        if show_all or a != b:
            lines.append(f"  {k:48s} {a:16.6g} -> {b:16.6g}  ({100 * rel:+.2f} %)" if a else f"  {k:48s} {a:16.6g} -> {b:16.6g}")
    if res["unusable"]:
        lines.append(f"unusable (all zero): {', '.join(res['unusable'])}")
    for n in res["notes"]:
        lines.append(f"note: {n}")
    lines.append(f"VERDICT: {res['verdict']} ({'; '.join(res['reasons'])})")
    return "\n".join(lines)


def stage_table(arms: dict) -> str:
    rows = []
    keys = sorted({k for r in arms["A"] for k in flatten(r)[1]})
    for k in keys:
        if k == "wall_s":
            continue
        a = statistics.median(flatten(r)[1].get(k, 0.0) for r in arms["A"])
        b = statistics.median(flatten(r)[1].get(k, 0.0) for r in arms["B"])
        if a == 0.0 and b == 0.0:
            continue
        rows.append(f"  {k:40s} {a:8.3f} -> {b:8.3f} s")
    return "stage seconds per rank (O, medians):\n" + "\n".join(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--upf", default=os.environ.get("DFT_SI_UPF"), help="Si UPF (default $DFT_SI_UPF)")
    ap.add_argument("--case", default="si32", choices=["si8", "si16", "si32", "si64"])
    ap.add_argument("--ranks", type=int, default=4)
    ap.add_argument("--mode", default="production", choices=["production", "standard", "strict"])
    ap.add_argument("--ladder", default="on", choices=["on", "off"])
    ap.add_argument("--pairs", type=int, default=3)
    ap.add_argument("--a", action="append", default=[], help="arm A environment NAME=value[,NAME=value]")
    ap.add_argument("--b", action="append", default=[], help="arm B environment NAME=value[,NAME=value]")
    ap.add_argument("--expect", action="append", default=[], help="counter:down|up that the idea must move")
    ap.add_argument("--min-drop", type=float, default=0.10, help="relative counter move that counts as large")
    ap.add_argument("--binary", help="use this test executable instead of building")
    ap.add_argument("--no-build", action="store_true", help="locate the binary with --offline only")
    ap.add_argument("--timeout", type=float, default=1800.0)
    ap.add_argument("--json-out", help="write the full result as JSON")
    ap.add_argument("--all-counters", action="store_true")
    ap.add_argument("--level", default="scf", choices=["scf", "solve"],
                    help="scf: whole SCF (stage 1, this harness); solve: one Davidson solve at a fixed "
                         "Hamiltonian on the same simulated ranks (stage 0, scripts/fast_screen_dav.py)")
    args = ap.parse_args()
    if not args.upf or not Path(args.upf).is_file():
        sys.exit("--upf (or $DFT_SI_UPF) must name the Si UPF file")
    if args.level == "solve":
        # Stage 0: delegate to the solve-level harness (counter names there: coll.<op>.calls,
        # dist.<stage>.calls; the SCF-level prefix "collectives." is mapped to "coll.").
        sys.path.insert(0, str(ROOT / "scripts"))
        import fast_screen_dav
        argv = ["--upf", args.upf, "--case", args.case, "--ranks", str(args.ranks),
                "--pairs", str(args.pairs)]
        for e in args.a:
            argv += ["--a", e]
        for e in args.b:
            argv += ["--b", e]
        for e in args.expect:
            argv += ["--expect", e.replace("collectives.", "coll.", 1)]
        if args.json_out:
            argv += ["--json-out", args.json_out]
        if args.all_counters:
            argv.append("--all-counters")
        return fast_screen_dav.main(argv)
    expect = []
    for e in args.expect:
        key, _, direction = e.rpartition(":")
        if direction not in ("down", "up") or not key:
            sys.exit(f"bad --expect {e!r} (want counter:down|up)")
        expect.append((key, direction))
    binary = Path(args.binary) if args.binary else find_binary(not args.no_build)
    digest = sha256(binary)
    base = {
        "DFT_SI_UPF": args.upf,
        "M2_PROFILE_CASE": args.case,
        "M2_PROFILE_RANKS": str(args.ranks),
        "M2_PROFILE_SERIAL": "off",
        "M2_PROFILE_LADDER": args.ladder,
    }
    if args.mode == "production":
        base["M2_PROFILE_MODE"] = "production"
    elif args.mode == "standard":
        base["M2_PROFILE_STANDARD"] = "on"
    env_a, env_b = parse_env(args.a), parse_env(args.b)
    print(f"binary {binary} sha256 {digest}")
    print(f"case {args.case} np={args.ranks} (simulated) mode {args.mode} ladder {args.ladder} pairs {args.pairs}")
    print(f"A: {env_a or '(defaults)'}  B: {env_b or '(defaults)'}")
    arms: dict[str, list] = {"A": [], "B": []}
    for i in range(args.pairs):
        order = ("A", "B") if i % 2 == 0 else ("B", "A")
        for arm in order:
            rec = run_once(binary, base, env_a if arm == "A" else env_b, args.timeout)
            arms[arm].append(rec)
            print(f"  pair {i + 1} {arm}: wall {rec['wall_s']:.3f} s load {rec['_load']:.2f} E {rec['energy_ha']:.12f} "
                  f"(SCF, Dav, H) ({rec['scf']}, {rec['davidson']}, {rec['h_columns']}) tree_sum {rec['collectives']['tree_sum'][0]}",
                  flush=True)
    limit = CLASS_LIMIT_HA["strict" if args.mode == "strict" else "standard"]
    res = verdict(arms, expect, limit, args.min_drop)
    print(report(res, args.all_counters))
    print(stage_table(arms))
    if args.json_out:
        res["binary_sha256"] = digest
        res["setup"] = {"case": args.case, "ranks": args.ranks, "mode": args.mode, "ladder": args.ladder,
                        "pairs": args.pairs, "A": env_a, "B": env_b}
        res["records"] = arms
        Path(args.json_out).write_text(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
