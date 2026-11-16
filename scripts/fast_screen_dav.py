#!/usr/bin/env python3
"""Stage-0 screen: ONE Davidson solve at a fixed Hamiltonian, environment A/B (FASTBENCH-DAVIDSON).

Protocol: docs/performance/fast-loop-protocol.md (section "Stage 0"). Standard library only,
python3 -I compatible, no network. Mac timings are observations (label O), never adoption
evidence; the counters are exact.

It runs crates/dft-cli/examples/davidson_bench.rs (built once, sha256 printed) on N alternating
pairs (A B, B A, ...). Each process builds the case once, then times `--repeats` cold solves
(atomic start, first-step tolerance 5e-3 Ha) and warm solves (previous orbitals kept, ratchet
tolerance) after one discarded warm-up each, and prints exact counters (H columns, Davidson
iterations, GEMM calls and FLOP per site/route, lean Davidson columns, REUSE-PSI and projection
cache counts, MPI collectives and stage calls on ranks) plus the eigenvalues and residuals.
Switches are read once per process, so every arm is its own process.

  python3 -I scripts/fast_screen_dav.py --upf /path/Si.pz-vbc.UPF --b BAANDI_DAVIDSON=classic
  python3 -I scripts/fast_screen_dav.py --upf ... --case si8 --b BAANDI_LEAN_PROJ_CACHE=on \
      --expect lean.proj_cache_columns:up
  python3 -I scripts/fast_screen_dav.py --upf ... --ranks 4 --b BAANDI_M2_PANEL_GS=fused \
      --expect coll.tree_sum.calls:down                      # simulated ranks (counters)
  python3 -I scripts/fast_screen_dav.py --upf ... --emit-input in.toml   # input for mpirun

Verdict per solve kind, then overall (REJECT > CANDIDATE > UNUSABLE > NEUTRAL):
  REJECT     a run failed; eigenvalues not reproducible between repeats or processes of one arm;
             an occupied eigenvalue differs from arm A by more than max(10 tol, 1e-10) Ha (tol =
             the solve's tolerance); fewer converged roots than A; an --expect counter (summed over
             kinds) moved the wrong way; median solve-time ratio B/A above 1 + spread
  CANDIDATE  median ratio below 1 - spread, or H columns / Davidson iterations / GEMM FLOP (or an
             --expect counter) moved by at least --min-drop with the time inside the spread
  UNUSABLE   every --expect counter is zero on both arms ("counter zero on both arms = unusable":
             a zero is never read as "no difference")
  NEUTRAL    otherwise
spread = max(--min-spread, half the range of the paired ratios). Identity: "bit-identical"
(eigenvalue FNV equal) or "rounding class" with max |d eps| over the occupied bands.
Not testable here (they depend on the SCF trajectory): VERIFY_REUSE, mixing, the cutoff ladder,
stop rules; use stage 1 (scripts/fast_screen.py, scripts/fast_screen_mpi.py).
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
sys.path.insert(0, str(ROOT / "scripts"))
import fast_screen  # noqa: E402  (case loader and the stage-1 displacement stream)

PRIMARY = ("h_columns", "davidson_iterations", "gemm_flop")


def find_binary(features: str | None, build: bool) -> Path:
    cmd = ["cargo", "build", "--release", "-p", "dft-cli", "--example", "davidson_bench",
           "--message-format=json"]
    if features:
        cmd += ["--features", features]
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
        if msg.get("reason") == "compiler-artifact" and msg.get("target", {}).get("name") == "davidson_bench":
            exe = msg.get("executable") or exe
    if not exe:
        sys.exit("could not locate the davidson_bench executable")
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


def case_path(name: str) -> Path:
    p = Path(name)
    if p.suffix == ".toml":
        return p if p.is_absolute() else ROOT / p
    return ROOT / "benchmarks" / "cases" / f"{name}-gamma.toml"


def run_once(cmd: list[str], env_extra: dict, timeout: float) -> dict:
    """One process; peak RSS from wait4 (bytes on macOS, KiB on Linux, as reported)."""
    env = dict(os.environ)
    env.update(env_extra)
    load = os.getloadavg()[0]
    import tempfile
    with tempfile.TemporaryFile("w+") as out, tempfile.TemporaryFile("w+") as err:
        t0 = time.perf_counter()
        proc = subprocess.Popen(cmd, env=env, stdout=out, stderr=err, text=True)
        deadline = t0 + timeout
        while True:
            pid, status, usage = os.wait4(proc.pid, os.WNOHANG)
            if pid:
                break
            if time.perf_counter() > deadline:
                proc.kill()
                os.wait4(proc.pid, 0)
                raise RuntimeError("timeout")
            time.sleep(0.01)
        wall = time.perf_counter() - t0
        code = os.waitstatus_to_exitcode(status) if hasattr(os, "waitstatus_to_exitcode") else status >> 8
        proc.returncode = code
        out.seek(0)
        err.seek(0)
        stdout, stderr = out.read(), err.read()
    if code != 0:
        raise RuntimeError(f"run failed (exit {code}):\n{stdout[-2000:]}\n{stderr[-2000:]}")
    recs = [json.loads(l.split(" ", 1)[1]) for l in stdout.splitlines() if l.startswith("DAVBENCH_JSON ")]
    if len(recs) != 1:
        raise RuntimeError(f"expected one DAVBENCH_JSON line, got {len(recs)}")
    rec = recs[0]
    rec["_load"] = load
    rec["_process_s"] = wall
    rec["_max_rss"] = usage.ru_maxrss  # bytes on macOS, KiB on Linux (OS convention)
    return rec


def kind_counters(k: dict) -> dict[str, float]:
    c = {"h_columns": float(k["h_columns"]), "davidson_iterations": float(k["davidson_iterations"])}
    for name, v in k.get("counters", {}).items():
        c[name] = float(v)
    return c


def verdict(arms: dict, expect: list[tuple[str, str]], min_drop: float, min_spread: float) -> dict:
    """Pure decision logic (unit-tested in scripts/tests/test_fast_screen_dav.py).

    arms: {"A": [DAVBENCH_JSON record, ...], "B": [...]}, same pair count, same kinds.
    """
    res: dict = {"kinds": {}, "notes": [], "reasons": []}
    kinds = [k["kind"] for k in arms["A"][0]["kinds"]]
    occupied = arms["A"][0]["occupied"]
    sums = {"A": {}, "B": {}}
    reject = []
    any_candidate = False
    for i, kind in enumerate(kinds):
        per = {arm: [r["kinds"][i] for r in recs] for arm, recs in arms.items()}
        kr: dict = {"reasons": []}
        # Determinism inside each arm: repeats (in-process) and processes.
        for arm, ks in per.items():
            if not all(k["repeats_identical"] for k in ks) or len({k["eigen_fnv"] for k in ks}) > 1:
                kr["reasons"].append(f"arm {arm}: eigenvalues not reproducible")
            cs = [kind_counters(k) for k in ks]
            bad = sorted(n for n in cs[0] if len({c.get(n) for c in cs}) > 1)
            if bad or not all(k["counters_deterministic"] for k in ks):
                res["notes"].append(f"{kind} arm {arm}: counters vary between repeats/processes: {', '.join(bad[:6])}")
        ca, cb = kind_counters(per["A"][0]), kind_counters(per["B"][0])
        for arm, c in (("A", ca), ("B", cb)):
            for n, v in c.items():
                sums[arm][n] = sums[arm].get(n, 0.0) + v
        counters, unusable = {}, []
        for n in sorted(set(ca) | set(cb)):
            a, b = ca.get(n, 0.0), cb.get(n, 0.0)
            if a == 0.0 and b == 0.0:
                unusable.append(n)
                continue
            counters[n] = (a, b, (b - a) / a if a else float("inf"))
        kr["counters"], kr["unusable"] = counters, unusable
        # Identity A vs B.
        ea, eb = per["A"][0], per["B"][0]
        tol = ea["tolerance_ha"]
        deps = max(abs(x - y) for x, y in zip(ea["eigenvalues"][:occupied], eb["eigenvalues"][:occupied]))
        gate = max(10.0 * tol, 1e-10)
        kr["identity"] = {
            "bit_identical": ea["eigen_fnv"] == eb["eigen_fnv"],
            "max_deps_occupied_ha": deps,
            "gate_ha": gate,
            "d_band_energy_ha": eb["band_energy_occupied"] - ea["band_energy_occupied"],
            "residual_max": (ea["residual_max"], eb["residual_max"]),
            "converged": (ea["converged"], eb["converged"]),
        }
        if deps > gate:
            kr["reasons"].append(f"occupied eigenvalue moved {deps:.3e} Ha > gate {gate:.1e}")
        if eb["converged"] < ea["converged"]:
            kr["reasons"].append(f"converged roots {ea['converged']} -> {eb['converged']}")
        # Time: per-process medians, paired in run order.
        ta = [k["solve_ms_median"] for k in per["A"]]
        tb = [k["solve_ms_median"] for k in per["B"]]
        ratios = [b / a for a, b in zip(ta, tb)]
        med = statistics.median(ratios)
        sp = max(min_spread, (max(ratios) - min(ratios)) / 2)
        inproc = max((max(k["solve_ms"]) - min(k["solve_ms"])) / k["solve_ms_median"]
                     for ks in per.values() for k in ks)
        kr["time"] = {"A_ms": ta, "B_ms": tb, "ratios": ratios, "median_ratio": med, "spread": sp,
                      "in_process_range": inproc}
        drops = [n for n in PRIMARY if n in counters and counters[n][2] <= -min_drop]
        if kr["reasons"]:
            kr["verdict"] = "REJECT"
        elif len(ratios) < 3:
            kr["verdict"] = "CANDIDATE" if drops else "NEUTRAL"
            kr["reasons"].append((f"counter drop {', '.join(drops)}; " if drops else "")
                                 + "time not judged with fewer than 3 pairs")
        elif med > 1 + sp:
            kr["verdict"] = "REJECT"
            kr["reasons"].append(f"slower: median ratio {med:.4f} > 1 + {sp:.4f}")
        elif med < 1 - sp:
            kr["verdict"] = "CANDIDATE"
            kr["reasons"].append(f"faster: median ratio {med:.4f} < 1 - {sp:.4f}")
        elif drops:
            kr["verdict"] = "CANDIDATE"
            kr["reasons"].append(f"counter drop {', '.join(drops)}; time within +-{sp:.4f}")
        else:
            kr["verdict"] = "NEUTRAL"
            kr["reasons"].append(f"median ratio {med:.4f} within +-{sp:.4f}, no primary counter drop")
        if kr["verdict"] == "REJECT":
            reject.append(kind)
        any_candidate |= kr["verdict"] == "CANDIDATE"
        res["kinds"][kind] = kr
    # Expected counters, summed over the kinds.
    moved, expect_unusable = [], []
    for key, direction in expect:
        a, b = sums["A"].get(key, 0.0), sums["B"].get(key, 0.0)
        if a == 0.0 and b == 0.0:
            expect_unusable.append(key)
            continue
        good = (b < a) if direction == "down" else (b > a)
        if not good:
            res["reasons"].append(f"expected counter {key} did not move {direction} ({a:.6g} -> {b:.6g})")
            reject.append("expect")
        elif a == 0.0 or abs(b - a) / a >= min_drop:
            moved.append(key)
    res["expect_unusable"] = expect_unusable
    if reject:
        res["verdict"] = "REJECT"
        res["reasons"] += [f"{k}: {'; '.join(res['kinds'][k]['reasons'])}" for k in reject if k in res["kinds"]]
    elif expect and len(expect_unusable) == len(expect):
        res["verdict"] = "UNUSABLE"
        res["reasons"].append("every expected counter is zero on both arms (unusable, not 'no difference')")
    elif any_candidate or moved:
        res["verdict"] = "CANDIDATE"
        res["reasons"] += [f"{k}: {res['kinds'][k]['reasons'][-1]}" for k in kinds
                           if res["kinds"][k]["verdict"] == "CANDIDATE"]
        if moved:
            res["reasons"].append(f"expected counters moved: {', '.join(moved)}")
    else:
        res["verdict"] = "NEUTRAL"
        res["reasons"] += [f"{k}: {res['kinds'][k]['reasons'][-1]}" for k in kinds]
    return res


def report(res: dict, show_all: bool) -> str:
    lines = []
    for kind, kr in res["kinds"].items():
        t, idn = kr["time"], kr["identity"]
        lines.append(f"[{kind}] solve ms (O, per-process medians) A {['%.1f' % x for x in t['A_ms']]} "
                     f"B {['%.1f' % x for x in t['B_ms']]}; ratios {['%.4f' % x for x in t['ratios']]} "
                     f"median {t['median_ratio']:.4f} spread +-{t['spread']:.4f} "
                     f"(in-process range up to {100 * t['in_process_range']:.1f} %)")
        lines.append(f"  identity: {'bit-identical' if idn['bit_identical'] else 'rounding class'}, "
                     f"max |d eps_occ| {idn['max_deps_occupied_ha']:.3e} Ha (gate {idn['gate_ha']:.1e}), "
                     f"d band energy {idn['d_band_energy_ha']:+.3e} Ha, residual max {idn['residual_max'][0]:.2e} -> "
                     f"{idn['residual_max'][1]:.2e}, converged {idn['converged'][0]} -> {idn['converged'][1]}")
        lines.append("  | counter | A | B | B/A |")
        lines.append("  |---|---|---|---|")
        for n, (a, b, rel) in kr["counters"].items():
            if show_all or a != b or n in PRIMARY or n in ("gemm_calls",):
                lines.append(f"  | {n} | {a:.6g} | {b:.6g} | {(1 + rel):.4f} |" if a else f"  | {n} | 0 | {b:.6g} | new |")
        if kr["unusable"]:
            groups: dict[str, list[str]] = {}
            for n in kr["unusable"]:
                groups.setdefault(n.split(".")[0] if n.startswith(("coll.", "dist.")) else n, []).append(n)
            shown = [g if v == [g] else f"{g}.* ({len(v)})" for g, v in groups.items()]
            lines.append(f"  unusable (all zero on both arms): {', '.join(shown)}")
        lines.append(f"  kind verdict: {kr['verdict']} ({'; '.join(kr['reasons'])})")
    for n in res["notes"]:
        lines.append(f"note: {n}")
    lines.append(f"VERDICT: {res['verdict']} ({'; '.join(res['reasons'])})")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--upf", default=os.environ.get("DFT_SI_UPF"), help="Si UPF (default $DFT_SI_UPF)")
    ap.add_argument("--case", default="si32", help="si8|si16|si32|si64 or a benchmarks/cases/*.toml path")
    ap.add_argument("--seed", type=int, default=601)
    ap.add_argument("--amp", type=float, default=1e-3, help="fractional displacement amplitude (stage-1 stream)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--ranks", type=int, default=1, help="simulated ranks (threads of one process, M2 path)")
    ap.add_argument("--pairs", type=int, default=3)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--kinds", default="cold,warm")
    ap.add_argument("--a", action="append", default=[], help="arm A environment NAME=value[,NAME=value]")
    ap.add_argument("--b", action="append", default=[], help="arm B environment NAME=value[,NAME=value]")
    ap.add_argument("--expect", action="append", default=[], help="counter:down|up (summed over kinds)")
    ap.add_argument("--min-drop", type=float, default=0.05)
    ap.add_argument("--min-spread", type=float, default=0.01)
    ap.add_argument("--density", help="base density file (density.bin of a converged dft-cli run)")
    ap.add_argument("--cli", help="dft-cli binary: run the case once (cached in --work) and use its density.bin")
    ap.add_argument("--bench-args", default="", help="extra davidson_bench arguments (space separated)")
    ap.add_argument("--binary", help="davidson_bench executable (default: build it)")
    ap.add_argument("--features", help="cargo features for the build (mkl)")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--work", help="work directory for the input and cached density")
    ap.add_argument("--emit-input", help="write the CLI input of the case to this path and exit")
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--json-out")
    ap.add_argument("--all-counters", action="store_true")
    a = ap.parse_args(argv)
    if not a.upf or not Path(a.upf).is_file():
        sys.exit("--upf (or $DFT_SI_UPF) must name the Si UPF file")
    case = fast_screen.load_case(case_path(a.case))
    text = fast_screen.case_input(case, Path(a.upf).resolve(), a.seed, a.amp, a.threads, [])
    if a.emit_input:
        Path(a.emit_input).write_text(text)
        print(f"wrote {a.emit_input} ({case['id']} seed {a.seed} amp {a.amp} T={a.threads})")
        return 0
    if a.work:
        work = Path(a.work)
        work.mkdir(parents=True, exist_ok=True)
    else:
        import tempfile
        work = Path(tempfile.mkdtemp(prefix="fast_screen_dav_"))
    src = work / f"{case['id']}-s{a.seed}-a{a.amp}-T{a.threads}.toml"
    src.write_text(text)
    density = a.density
    if a.cli and not density:
        key = hashlib.sha256(text.encode()).hexdigest()[:12]
        out = work / f"cli-{key}"
        if not (out / "density.bin").exists():
            print(f"base density: running {a.cli} once on {src.name} ...", flush=True)
            r = subprocess.run([a.cli, "run", str(src), "--output-dir", str(out)], capture_output=True, text=True)
            if r.returncode != 0:
                sys.exit(f"dft-cli failed (exit {r.returncode}): {r.stderr[-1000:]}")
        density = str(out / "density.bin")
    expect = []
    for e in a.expect:
        key, _, direction = e.rpartition(":")
        if direction not in ("down", "up") or not key:
            sys.exit(f"bad --expect {e!r} (want counter:down|up)")
        expect.append((key, direction))
    binary = Path(a.binary) if a.binary else find_binary(a.features, not a.no_build)
    digest = sha256(binary)
    cmd = [str(binary), "--input", str(src), "--repeats", str(a.repeats), "--kinds", a.kinds]
    if a.ranks > 1:
        cmd += ["--sim-ranks", str(a.ranks)]
    if density:
        cmd += ["--density", density]
    cmd += a.bench_args.split()
    env_a, env_b = parse_env(a.a), parse_env(a.b)
    print(f"binary {binary} sha256 {digest[:16]}")
    print(f"case {case['id']} seed {a.seed} amp {a.amp} T={a.threads} ranks {a.ranks} pairs {a.pairs} "
          f"repeats {a.repeats} kinds {a.kinds}; base density {density or 'atomic superposition'}")
    print(f"A: {env_a or '(defaults)'}  B: {env_b or '(defaults)'}")
    t_screen = time.perf_counter()
    arms: dict[str, list] = {"A": [], "B": []}
    try:
        for i in range(a.pairs):
            for arm in (("A", "B") if i % 2 == 0 else ("B", "A")):
                rec = run_once(cmd, env_a if arm == "A" else env_b, a.timeout)
                arms[arm].append(rec)
                ks = " ".join(f"{k['kind']} {k['solve_ms_median']:.1f} ms (Dav {k['davidson_iterations']}, "
                              f"H {k['h_columns']})" for k in rec["kinds"])
                print(f"  pair {i + 1} {arm}: process {rec['_process_s']:.2f} s (setup {rec['setup_s']:.2f} s) "
                      f"load {rec['_load']:.1f} peak RSS {rec['_max_rss']}: {ks}", flush=True)
    except RuntimeError as e:
        print(f"VERDICT: REJECT (a run failed: {str(e)[:600]})")
        return 0
    r0 = arms["A"][0]
    print(f"setup: bands {r0['bands']} occupied {r0['occupied']} plane waves {r0['plane_waves']} grid {r0['grid']} "
          f"gamma_real {r0['gamma_real']} base {r0['base_density']}; E_H change {r0['e_h_change_ha']:.3e} Ha, "
          f"warm tol {r0['warm_tol_ha']:.3e} Ha; screen wall {time.perf_counter() - t_screen:.1f} s")
    res = verdict(arms, expect, a.min_drop, a.min_spread)
    print(report(res, a.all_counters))
    if a.json_out:
        res["binary_sha256"] = digest
        res["setup"] = {"case": case["id"], "seed": a.seed, "amp": a.amp, "threads": a.threads, "ranks": a.ranks,
                        "pairs": a.pairs, "repeats": a.repeats, "A": env_a, "B": env_b, "density": density}
        res["records"] = arms
        Path(a.json_out).write_text(json.dumps(res, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
