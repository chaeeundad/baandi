#!/usr/bin/env python3 -I
"""Compare a rust-dft k-point SCF (results.json) with a QE pw.x run.

QE side: either a full pw.x output (--qe-out, verbosity high) or a summary JSON
(--qe-summary) with keys energy_ry, nk, iterations, kpoints_cryst [[x,y,z,wk]...],
bands_ev [[8 values]...] in QE k order. Standard library only.
Exit codes: 0 all tolerances met, 5 comparison done and a tolerance exceeded,
3 inconsistent inputs, 1 usage or I/O error.
Defaults (printed; any override is flagged): total energy 1e-4 Ha/atom,
eigenvalues 1e-3 Ha, k weights 1e-6. Never loosened silently.
"""
import argparse
import json
import re
import sys

HA_PER_RY = 0.5
EV_PER_HA = 27.211386245988  # unit conversion for the QE eV printout; CODATA 2018 value used by uom (docs/decisions/0012)
DEFAULTS = {"energy_ha_per_atom": 1e-4, "eigenvalue_ha": 1e-3, "weight": 1e-6}


class Inconsistent(Exception):
    pass


def parse_qe_out(path):
    txt = open(path).read()
    m = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry", txt, re.M)
    if len(m) != 1:
        raise Inconsistent("expected exactly one '!' total energy line")
    nk = int(re.search(r"number of k points=\s*(\d+)", txt).group(1))
    it = re.search(r"convergence has been achieved in\s+(\d+) iterations", txt)
    if not it:
        raise Inconsistent("QE run did not report convergence")
    sec = re.split(r"cryst\. coord\.\n", txt)
    if len(sec) < 2:
        raise Inconsistent("no crystal k list (needs verbosity high)")
    ks = re.findall(r"k\(\s*\d+\) = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\), wk =\s*([\d.]+)", sec[1])[:nk]
    blocks = re.split(r"\n\s*k =", txt.split("End of self-consistent calculation")[1])[1:]
    bands = []
    for b in blocks:
        rest = b.split("bands (ev):")[1].split("occupation numbers")[0]
        bands.append([float(x) for x in re.findall(r"-?\d+\.\d+", rest)])
    return {"energy_ry": float(m[0]), "nk": nk, "iterations": int(it.group(1)),
            "kpoints_cryst": [[float(a), float(b), float(c), float(w)] for a, b, c, w in ks],
            "bands_ev": bands}


def same_k(a, b, sign):
    return all(abs(((sign * x - y) + 0.5) % 1.0 - 0.5) < 1e-5 for x, y in zip(a, b))


def compare(ours, qe, tol, natoms):
    if "kpoints" in ours:
        pts, count = ours["kpoints"]["points"], ours["kpoints"]["count"]
    else:
        # legacy Gamma-only results.json has no k-point block: one Gamma point of weight 1
        pts = [{"fractional": [0.0, 0.0, 0.0], "weight": 1.0, "eigenvalues": ours["eigenvalues"]}]
        count = 1
    if count != qe["nk"] or len(pts) != qe["nk"]:
        raise Inconsistent(f"k count differs: ours {len(pts)} QE {qe['nk']}")
    if len(qe["bands_ev"]) != qe["nk"] or len(qe["kpoints_cryst"]) != qe["nk"]:
        raise Inconsistent("QE summary lengths inconsistent")
    # QE prints wk normalized to a sum of 2 (spin degeneracy); ours sum to 1: compare normalized weights.
    wsum = sum(k[3] for k in qe["kpoints_cryst"])
    if wsum <= 0:
        raise Inconsistent("QE weights do not sum to a positive number")
    used = set()
    max_eig = 0.0
    max_w = 0.0
    for p in pts:
        hits = [i for i, k in enumerate(qe["kpoints_cryst"])
                if i not in used and (same_k(p["fractional"], k[:3], 1) or same_k(p["fractional"], k[:3], -1))]
        if len(hits) != 1:
            raise Inconsistent(f"k {p['fractional']} matches {len(hits)} QE points")
        i = hits[0]
        used.add(i)
        q = qe["bands_ev"][i]
        e = p["eigenvalues"]["ev"]
        if len(q) > len(e):
            raise Inconsistent("QE has more bands than ours")
        max_eig = max(max_eig, max(abs(a - b) for a, b in zip(q, e)) / EV_PER_HA)
        max_w = max(max_w, abs(p["weight"] - qe["kpoints_cryst"][i][3] / wsum))
    e_ours_ry = None
    te = ours["total_energy"]
    ha = te["ha"] if isinstance(te, dict) else te
    e_diff = (ha - qe["energy_ry"] * HA_PER_RY) / natoms
    res = {"k_points": qe["nk"], "max_weight_diff": max_w, "max_eigenvalue_diff_ha": max_eig,
           "total_energy_diff_ha_per_atom": e_diff, "ours_total_ha": ha, "qe_total_ry": qe["energy_ry"],
           "tolerances": tol}
    ok = (abs(e_diff) <= tol["energy_ha_per_atom"] and max_eig <= tol["eigenvalue_ha"] and max_w <= tol["weight"])
    return ok, res


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ours", required=True, help="rust-dft results.json (k-point run)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--qe-out")
    g.add_argument("--qe-summary")
    ap.add_argument("--natoms", type=int, default=2)
    ap.add_argument("--energy-tol", type=float, default=DEFAULTS["energy_ha_per_atom"])
    ap.add_argument("--eig-tol", type=float, default=DEFAULTS["eigenvalue_ha"])
    ap.add_argument("--weight-tol", type=float, default=DEFAULTS["weight"])
    a = ap.parse_args(argv)
    tol = {"energy_ha_per_atom": a.energy_tol, "eigenvalue_ha": a.eig_tol, "weight": a.weight_tol}
    for k, v in tol.items():
        if v != DEFAULTS[k]:
            print(f"NON-DEFAULT tolerance {k}={v} (default {DEFAULTS[k]})")
    try:
        ours = json.load(open(a.ours))
        qe = parse_qe_out(a.qe_out) if a.qe_out else json.load(open(a.qe_summary))
        ok, res = compare(ours, qe, tol, a.natoms)
    except Inconsistent as e:
        print("INCONSISTENT:", e)
        return 3
    except OSError as e:
        print("I/O error:", e)
        return 1
    print(json.dumps(res, indent=2))
    print("PASS" if ok else "FAIL (tolerance exceeded)")
    return 0 if ok else 5


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
