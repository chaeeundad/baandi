#!/usr/bin/env python3
"""R3 gate evaluation: BAANDI rpa.json against an ABINIT gwrpacorr run at
matched settings. The ABINIT values are supplied in a small JSON record
written by the operator from the ABINIT output (the lines quoted verbatim, so
the reader can check them):

  {"code": "ABINIT", "version": "...", "command": "...", "host": "dev08",
   "input": "si2_444.abi", "pseudopotential_sha256": {"Si": "..."},
   "ecut_chi_hartree": [4.0, 5.0, 6.0, 8.0],
   "e_c_hartree_per_cell": [...],
   "quoted_lines": ["..."]}

Gates (plan R3): |E_c(BAANDI) - E_c(ABINIT)| <= 1 meV/atom per matched
cutoff; extrapolated values (same (E_cut^chi)^-3/2 fit on both) within
2 meV/atom. Writes <out>.json with the per-point differences and pass/fail.
Python 3.9+, standard library only.

  python3 -I compare_r3.py OUT/rpa.json abinit-si2_444.json --out si2_444-compare.json
"""
import argparse
import json
from pathlib import Path

MEV = 27211.386245988


def fit(c, e):
    x = [v ** -1.5 for v in c]
    n = len(x)
    xm, ym = sum(x) / n, sum(e) / n
    sxx = sum((a - xm) ** 2 for a in x)
    slope = sum((a - xm) * (b - ym) for a, b in zip(x, e)) / sxx
    return ym - slope * xm, slope


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rpa_json")
    ap.add_argument("oracle_json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    b = json.loads(Path(a.rpa_json).read_text())
    o = json.loads(Path(a.oracle_json).read_text())
    atoms = b["atoms"]
    bc = b["settings"]["ecut_chi_series_hartree"]
    be = b["e_c_hartree_per_cell"]
    oc, oe = o["ecut_chi_hartree"], o["e_c_hartree_per_cell"]
    for sp, h in o.get("pseudopotential_sha256", {}).items():
        mine = {p["species"]: p["sha256"] for p in b["pseudopotentials"]}.get(sp)
        if mine != h:
            print(f"note: {sp} hash differs (BAANDI {mine}, oracle {h}); psp8/UPF conversions must record both")
    rows, ok = [], True
    for c, e in zip(bc, be):
        if c not in oc:
            continue
        d = (e - oe[oc.index(c)]) * MEV / atoms
        rows.append({"ecut_chi_hartree": c, "baandi_ha": e, "oracle_ha": oe[oc.index(c)],
                     "diff_mev_per_atom": d, "pass_1mev": abs(d) <= 1.0})
        ok &= abs(d) <= 1.0
    result = {"matched": rows, "matched_pass": ok and bool(rows)}
    common = [c for c in bc if c in oc]
    if len(common) >= 2:
        bi, _ = fit(common, [be[bc.index(c)] for c in common])
        oi, _ = fit(common, [oe[oc.index(c)] for c in common])
        d = (bi - oi) * MEV / atoms
        result["extrapolated"] = {"baandi_ha": bi, "oracle_ha": oi,
                                  "diff_mev_per_atom": d, "pass_2mev": abs(d) <= 2.0}
    result["oracle"] = {k: o.get(k) for k in ("code", "version", "command", "host", "input")}
    Path(a.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
