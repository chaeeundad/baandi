#!/usr/bin/env python3
"""Paired summary of Phase 1 raw results (relax-precon-run-1 JSON files). Stdlib only.

  python3 -I summarize.py RESULTS_DIR [--baseline baandi-bfgs:warm] [--json OUT.json]

Per case and arm (method:start): success count, medians of attempted evaluations
(to 0.1 / 0.05 / 0.02 eV/A), failed/rejected/auxiliary evaluations, SCF
iterations, H vectors and wall. Paired comparisons use the same (case, seed):
median ratio of attempted evaluations against the baseline arm, wins/ties/losses.
Wall ratios are printed ONLY between arms with the same start class (cold vs cold,
warm vs warm), never across (spec section 6). Final energies are compared per
seed against the baseline; |dE| above 1e-4 Ha or a final-structure difference above
0.05 Bohr (after removing the rigid translation) is flagged as a possibly different minimum.
"""
import argparse
import glob
import json
import math
import os
import statistics as st
from collections import defaultdict


def med(xs):
    xs = [x for x in xs if x is not None]
    return st.median(xs) if xs else None


def fmt(x, nd=1):
    if x is None:
        return "-"
    return ("%%.%df" % nd) % x if isinstance(x, float) else str(x)


def max_disp(a, b):
    """Max atomic displacement between two final structures after removing the
    mean (rigid) translation, which is a zero mode of a periodic cell."""
    if not a or not b:
        return None
    d = [[x - y for x, y in zip(p, q)] for p, q in zip(a, b)]
    m = [sum(v[k] for v in d) / len(d) for k in range(3)]
    return max(math.sqrt(sum((v[k] - m[k]) ** 2 for k in range(3))) for v in d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--baseline", default="baandi-bfgs:warm")
    ap.add_argument("--json")
    args = ap.parse_args()
    runs = []
    for f in sorted(glob.glob(os.path.join(args.results, "*.json"))):
        with open(f) as fh:
            r = json.load(fh)
        if r.get("schema") == "relax-precon-run-1":
            runs.append(r)
    by = defaultdict(dict)  # (case, arm) -> seed -> run
    for r in runs:
        by[(r["case"], "%s:%s" % (r["method"], r["start"]))][r["seed"]] = r
    out = {"per_arm": [], "paired": []}
    for case in sorted({c for c, _ in by}):
        arms = sorted(a for c, a in by if c == case)
        print("\n== %s" % case)
        print("%-24s %4s %6s %6s %6s %6s %5s %5s %4s %7s %8s %8s" % (
            "arm", "ok", "e0.1", "e0.05", "e0.02", "evals", "fail", "rej", "aux", "scf", "hvec", "wall_s"))
        for arm in arms:
            rs = list(by[(case, arm)].values())
            ss = [r["summary"] for r in rs]
            ok = sum(1 for s in ss if s["converged_at_evaluation"])
            row = {
                "case": case, "arm": arm, "runs": len(rs), "converged": ok,
                "median_evals_to": {t: med([s["evaluations_to_threshold_ev_a"][t] for s in ss])
                                    for t in ("0.1", "0.05", "0.02")},
                "median_attempted": med([s["attempted_evaluations"] for s in ss]),
                "failed": sum(s["failed_evaluations"] for s in ss),
                "rejected": sum(s["rejected_evaluations"] for s in ss),
                "auxiliary": sum(s["auxiliary_evaluations"] for s in ss),
                "median_scf": med([s["total_scf_iterations"] for s in ss]),
                "median_h_vectors": med([s["total_h_vectors"] for s in ss]),
                "median_wall_s": med([s["end_to_end_wall_s"] for s in ss]),
            }
            out["per_arm"].append(row)
            m = row["median_evals_to"]
            print("%-24s %2d/%-2d %6s %6s %6s %6s %5d %5d %4d %7s %8s %8s" % (
                arm, ok, len(rs), fmt(m["0.1"]), fmt(m["0.05"]), fmt(m["0.02"]), fmt(row["median_attempted"]),
                row["failed"], row["rejected"], row["auxiliary"], fmt(row["median_scf"]),
                fmt(row["median_h_vectors"]), fmt(row["median_wall_s"])))
        base = by.get((case, args.baseline))
        if not base:
            continue
        bstart = args.baseline.split(":")[1]
        print("paired vs %s (evaluations to 0.02 eV/A; wall only within the same start class)" % args.baseline)
        for arm in arms:
            if arm == args.baseline:
                continue
            ratios, wall_ratios, w, t, l, flags = [], [], 0, 0, 0, []
            for seed, r in by[(case, arm)].items():
                b = base.get(seed)
                if not b:
                    continue
                a_n, b_n = r["summary"]["converged_at_evaluation"], b["summary"]["converged_at_evaluation"]
                if a_n and b_n:
                    ratios.append(a_n / b_n)
                    w, t, l = w + (a_n < b_n), t + (a_n == b_n), l + (a_n > b_n)
                    if arm.split(":")[1] == bstart:
                        wall_ratios.append(r["summary"]["end_to_end_wall_s"] / b["summary"]["end_to_end_wall_s"])
                    de = abs(r["summary"]["final_energy_ha"] - b["summary"]["final_energy_ha"])
                    dx = max_disp(r["summary"]["final_positions_bohr"], b["summary"]["final_positions_bohr"])
                    if de > 1e-4 or (dx is not None and dx > 0.05):
                        flags.append({"seed": seed, "dE_ha": de, "max_dx_bohr": dx})
            row = {"case": case, "arm": arm, "baseline": args.baseline, "pairs": len(ratios),
                   "median_eval_ratio": med(ratios), "wins_ties_losses": [w, t, l],
                   "median_wall_ratio_same_start": med(wall_ratios), "different_minimum_flags": flags}
            out["paired"].append(row)
            print("  %-24s pairs %d  eval ratio %s  W/T/L %d/%d/%d  wall ratio %s%s" % (
                arm, len(ratios), fmt(row["median_eval_ratio"], 2), w, t, l,
                fmt(row["median_wall_ratio_same_start"], 2) if wall_ratios else "n/a (other start class)",
                "  DIFFERENT-MINIMUM? %s" % flags if flags else ""))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
