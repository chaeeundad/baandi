#!/usr/bin/env python3 -I
"""Fast-loop stage-1 screen: one case, environment (or binary) A/B, alternating pairs.

Protocol: docs/performance/fast-loop-protocol.md (section "Harness"). Standard library only,
no network. Screening numbers are never adoption evidence (Mac timings are observations, O).

  python3 -I scripts/fast_screen.py --binary target/release/dft-cli --upf /path/Si.pz-vbc.UPF \
      --env-b BAANDI_PROJECT_WIDE=on                       # A = environment as is, B = with the switch
  python3 -I scripts/fast_screen.py --binary old/dft-cli --binary-b new/dft-cli --upf ...   # two binaries
  ... --case benchmarks/cases/si16-gamma.toml --seed 601 --amp 1e-3 --threads 1 --pairs 3

The input is built from a benchmarks/cases/*.toml file like scripts/bench_suite.py does, with
the single-path defaults (no scf block beyond seed/threads/max_iterations unless --scf adds
lines). Displacement: every atom, every fractional component += random.Random(seed).uniform(
-amp, amp) in atom order (the documented MPI-M1-GATE3 stream; byte equality with dev08's
p48lib geometry is not verified, so screen numbers are not comparable with dev08 tables).
amp = 0 keeps the case unperturbed.

Pairs alternate the arm order (AB, BA, AB, ...) to cancel drift. The verdict:
  REJECT     a run failed, did not converge, failed verification, |E_B - E_A| > --energy-gate
             (default 5e-9 Ha, the standard class), or the median paired wall ratio is above
             1 + spread (slower beyond noise)
  CANDIDATE  the median paired ratio is below 1 - spread, or a work counter (H vectors,
             Davidson iterations, GEMM GFLOP) drops by at least --counter-drop while the
             wall is within the spread
  NEUTRAL    otherwise
spread = half the range of the paired ratios (at least --min-spread, default 0.01).
Exit code 0 always when the runs completed (the verdict is printed), 1 on usage/IO error.
"""
import argparse, hashlib, json, math, os, random, statistics, subprocess, sys, time
from pathlib import Path

KEYS_REQUIRED = ("id", "lattice", "atoms", "cutoff_ry")


def load_case(path):
    case = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        key, sep, value = line.partition("=")
        if sep:
            case[key.strip()] = json.loads(value.strip())
    for k in KEYS_REQUIRED:
        if k not in case:
            raise ValueError(f"{path}: missing {k}")
    return case


def displaced(atoms, seed, amp):
    if amp == 0:
        return [list(a) for a in atoms]
    rng = random.Random(seed)
    return [[(x + rng.uniform(-amp, amp)) % 1.0 for x in a] for a in atoms]


def case_input(case, upf, seed, amp, threads, scf_extra, xc="lda-pz"):
    lines = ["[system]", f"name = {json.dumps(case['id'])}", 'lattice_unit = "bohr"',
             f"lattice_vectors = {json.dumps(case['lattice'])}"]
    for site in displaced(case["atoms"], seed, amp):
        lines += ["[[system.atoms]]", 'species = "Si"', f"fractional = {json.dumps(site)}"]
    lines += ["[pseudopotentials]", f"Si = {json.dumps(str(upf))}", "[basis]",
              f"ecutwfc_ry = {case['cutoff_ry']}"]
    if "kgrid" in case:
        lines += [f"kgrid = {json.dumps(case['kgrid'])}",
                  f"kshift = {json.dumps(case.get('kshift', [0, 0, 0]))}"]
    else:
        lines += ["kpoints = [[0.0, 0.0, 0.0]]"]
    if "fft_grid" in case:
        lines += [f"fft_grid = {json.dumps(case['fft_grid'])}"]
    lines += ["[xc]", f"functional = {json.dumps(xc)}", "[scf]", f"seed = {seed}", f"threads = {threads}",
              f"max_iterations = {case.get('max_iterations', 100)}"]
    if "nbands" in case:
        lines += [f"nbands = {case['nbands']}"]
    lines += list(scf_extra)
    return "\n".join(lines) + "\n"


def parse_env(items):
    env = {}
    for item in items or []:
        k, sep, v = item.partition("=")
        if not sep or not k:
            raise ValueError(f"bad --env item {item!r} (KEY=VALUE)")
        env[k] = v
    return env


def get(d, *path, default=None):
    for p in path:
        if not isinstance(d, dict) or p not in d:
            return default
        d = d[p]
    return d


def counters(r):
    w = get(r, "convergence", "work", default={}) or {}
    rows = get(w, "gemm_dispatch", "rows", default=[]) or []
    lad = get(w, "cutoff_ladder", default={}) or {}
    ver = r.get("verification") or {}
    c = {
        "h_vectors": (w.get("h_vectors_recorded") or 0) + (w.get("certificate_h_vectors") or 0),
        "h_vectors_recorded": w.get("h_vectors_recorded") or 0,
        "certificate_h_vectors": w.get("certificate_h_vectors") or 0,
        "coarse_h_vectors": lad.get("coarse_h_vectors") or 0,
        "davidson_iterations": (w.get("davidson_iterations_recorded") or 0)
        + (w.get("certificate_davidson_iterations") or 0),
        "coarse_davidson_iterations": lad.get("coarse_davidson_iterations") or 0,
        "scf_steps": get(r, "convergence", "iterations", default=0) or 0,
        "coarse_scf_steps": lad.get("coarse_iterations") or 0,
        "verification_attempts": len(ver.get("attempts") or []),
        "gemm_calls": sum(x.get("calls", 0) for x in rows),
        "gemm_gflop": round(sum(x.get("gflop", 0.0) for x in rows), 6),
    }
    return c


def run_once(binary, src, outdir, env_extra, timeout):
    env = dict(os.environ)
    env.update(env_extra)
    load = os.getloadavg()[0]
    t0 = time.perf_counter()
    p = subprocess.run([binary, "run", str(src), "--output-dir", str(outdir)], env=env,
                       capture_output=True, text=True, timeout=timeout)
    wall = time.perf_counter() - t0
    rec = {"exit": p.returncode, "wall_s": wall, "load": load}
    if p.returncode != 0:
        rec["stderr"] = p.stderr[-600:]
        return rec
    r = json.loads((Path(outdir) / "results.json").read_text())
    rec["total_s"] = get(r, "timings", "total_s")
    rec["verification_s"] = get(r, "timings", "verification_s")
    rec["converged"] = get(r, "convergence", "converged")
    rec["verified"] = r.get("verified")
    rec["energy_ha"] = get(r, "total_energy", "ha")
    dens = Path(outdir) / "density.bin"
    rec["density_sha256"] = hashlib.sha256(dens.read_bytes()).hexdigest() if dens.exists() else None
    rec["counters"] = counters(r)
    rec["timings"] = r.get("timings")
    return rec


def fmt(x, nd=3):
    return "-" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--binary", required=True)
    p.add_argument("--binary-b", help="second binary for arm B (default: --binary)")
    p.add_argument("--upf", required=True)
    p.add_argument("--xc", default="lda-pz", help="xc.functional of the input (lda-pz, pbe); must match the UPF")
    p.add_argument("--case", default="benchmarks/cases/si32-gamma.toml")
    p.add_argument("--seed", type=int, default=601)
    p.add_argument("--amp", type=float, default=1e-3, help="fractional displacement amplitude")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--pairs", type=int, default=3)
    p.add_argument("--env-a", action="append", default=[], help="KEY=VALUE for arm A (repeatable)")
    p.add_argument("--env-b", action="append", default=[], help="KEY=VALUE for arm B (repeatable)")
    p.add_argument("--scf", action="append", default=[], help="extra line in [scf], both arms (repeatable)")
    p.add_argument("--scf-a", action="append", default=[], help="extra line in [scf], arm A only")
    p.add_argument("--scf-b", action="append", default=[], help="extra line in [scf], arm B only")
    p.add_argument("--out", default=None, help="work directory (default: a fresh temp dir)")
    p.add_argument("--timeout", type=float, default=1800)
    p.add_argument("--energy-gate", type=float, default=5e-9)
    p.add_argument("--counter-drop", type=float, default=0.05)
    p.add_argument("--min-spread", type=float, default=0.01)
    p.add_argument("--max-load", type=float, default=None,
                   help="wait (up to 120 s) until the 1-min load is below this before each run")
    p.add_argument("--json", default=None, help="write all records to this file")
    a = p.parse_args(argv)
    try:
        env_a, env_b = parse_env(a.env_a), parse_env(a.env_b)
        case = load_case(a.case)
    except (ValueError, OSError) as e:
        print("error:", e, file=sys.stderr)
        return 1
    bin_a = str(Path(a.binary).resolve())
    bin_b = str(Path(a.binary_b).resolve()) if a.binary_b else bin_a
    if a.out:
        work = Path(a.out)
        work.mkdir(parents=True, exist_ok=True)
    else:
        import tempfile
        work = Path(tempfile.mkdtemp(prefix="fast_screen_"))
    srcs = {}
    for arm, extra in (("A", a.scf_a), ("B", a.scf_b)):
        srcs[arm] = work / f"input_{arm}.toml"
        srcs[arm].write_text(case_input(case, Path(a.upf).resolve(), a.seed, a.amp, a.threads,
                                        a.scf + extra, a.xc))
    src = srcs["A"]
    sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()[:16]
    print(f"case {case['id']} seed {a.seed} amp {a.amp} T={a.threads} pairs {a.pairs}; input {src}")
    print(f"A: {bin_a} ({sha(bin_a)}) env {env_a or '{}'} scf {a.scf_a or '[]'}")
    print(f"B: {bin_b} ({sha(bin_b)}) env {env_b or '{}'} scf {a.scf_b or '[]'}")
    arms = {"A": (bin_a, env_a), "B": (bin_b, env_b)}
    recs = {"A": [], "B": []}
    for i in range(a.pairs):
        order = ("A", "B") if i % 2 == 0 else ("B", "A")
        for arm in order:
            if a.max_load is not None:
                t_wait = time.time()
                while os.getloadavg()[0] > a.max_load and time.time() - t_wait < 120:
                    time.sleep(5)
            binary, env = arms[arm]
            rec = run_once(binary, srcs[arm], work / f"{arm}{i}", env, a.timeout)
            recs[arm].append(rec)
            print(f"  pair {i} {arm}: exit {rec['exit']} wall {rec['wall_s']:.3f} s "
                  f"load {rec['load']:.2f} E {rec.get('energy_ha')!r}", flush=True)
            if rec["exit"] != 0:
                print("    stderr:", rec.get("stderr", "").strip()[-300:])
    if a.json:
        Path(a.json).write_text(json.dumps({"args": vars(a), "records": recs}, indent=1))
    return report(recs, a)


def report(recs, a):
    reasons = []
    ok = all(r["exit"] == 0 for arm in recs.values() for r in arm)
    if not ok:
        print("VERDICT: REJECT (a run failed)")
        return 0
    for arm, rs in recs.items():
        for r in rs:
            if not r.get("converged"):
                reasons.append(f"{arm} not converged")
            if r.get("verified") is False:
                reasons.append(f"{arm} verification failed")
    ea = [r["energy_ha"] for r in recs["A"]]
    eb = [r["energy_ha"] for r in recs["B"]]
    de = max(abs(x - y) for x in ea for y in eb)
    shas = {arm: sorted({r["density_sha256"] for r in rs}) for arm, rs in recs.items()}
    if de > a.energy_gate:
        reasons.append(f"|dE| {de:.3e} Ha above gate {a.energy_gate:.1e}")
    wa = [r["wall_s"] for r in recs["A"]]
    wb = [r["wall_s"] for r in recs["B"]]
    ratios = [b / x for x, b in zip(wa, wb)]
    med_ratio = statistics.median(ratios)
    spread = max(a.min_spread, (max(ratios) - min(ratios)) / 2)
    print()
    print("| arm | median wall s | min | max | range/median |")
    print("|---|---|---|---|---|")
    for arm, w in (("A", wa), ("B", wb)):
        m = statistics.median(w)
        print(f"| {arm} | {m:.3f} | {min(w):.3f} | {max(w):.3f} | {(max(w) - min(w)) / m:.3f} |")
    print(f"paired ratios B/A: {', '.join(f'{x:.4f}' for x in ratios)}; median {med_ratio:.4f}; "
          f"spread +-{spread:.4f}")
    print(f"energy A {ea[0]!r} B {eb[0]!r} max|dE| {de:.3e} Ha")
    same = shas["A"] == shas["B"]
    print(f"density sha256 A {[s[:12] for s in shas['A']]} B {[s[:12] for s in shas['B']]} "
          f"({'identical' if same else 'DIFFERENT'})")
    ca = recs["A"][0]["counters"]
    cb = recs["B"][0]["counters"]
    det = all(r["counters"] == ca for r in recs["A"]) and all(r["counters"] == cb for r in recs["B"])
    print(f"counters deterministic within each arm: {det}")
    print("| counter | A | B | B/A |")
    print("|---|---|---|---|")
    drops = []
    for k in ca:
        x, y = ca[k], cb.get(k, 0)
        if x == 0 and y == 0:
            print(f"| {k} | 0 | 0 | unusable (all zero) |")
            continue
        rel = (y / x) if x else math.inf
        print(f"| {k} | {x} | {y} | {rel:.4f} |")
        if k in ("h_vectors", "davidson_iterations", "gemm_gflop") and x and rel <= 1 - a.counter_drop:
            drops.append(f"{k} {rel:.3f}")
    ta = recs["A"][0].get("timings") or {}
    print("timings (first run per arm, s): " + ", ".join(
        f"{k} {fmt(get(recs['A'][0], 'timings', *path))}->{fmt(get(recs['B'][0], 'timings', *path))}"
        for k, path in (("total", ("total_s",)), ("setup", ("setup_s",)),
                        ("ladder_coarse", ("closure", "ladder_coarse_s")),
                        ("diag", ("phases", "diagonalization_s")),
                        ("verification", ("verification_s",)))) if ta else "")
    if reasons:
        verdict = "REJECT (" + "; ".join(reasons) + ")"
    elif len(ratios) < 3:
        verdict = ("CANDIDATE (counter drop: " + ", ".join(drops) + "; wall not judged, pairs < 3)"
                   if drops else "NEUTRAL (wall not judged with fewer than 3 pairs; counters only)")
    elif med_ratio > 1 + spread:
        verdict = f"REJECT (slower: median ratio {med_ratio:.4f} > 1 + {spread:.4f})"
    elif med_ratio < 1 - spread:
        verdict = f"CANDIDATE (wall: median ratio {med_ratio:.4f} < 1 - {spread:.4f})"
    elif drops:
        verdict = "CANDIDATE (counter drop: " + ", ".join(drops) + "; wall within spread)"
    else:
        verdict = f"NEUTRAL (median ratio {med_ratio:.4f} within +-{spread:.4f}, no counter drop)"
    print("VERDICT:", verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
