#!/usr/bin/env python3
"""R1 oracle comparison: BAANDI nscf.json against QE pw.x nscf
data-file-schema.xml (eigenvalues in Hartree, k order of the explicit list).
Reports per run: plane-wave counts per k (must match), max |dE| over all bands
up to N_b absolute, and after removing the constant offset of the highest
occupied eigenvalue at the first k (the ground-state/convention offset,
recorded separately per the R1 gate). Exit 1 when the aligned maximum exceeds
--gate (default 2e-5 Ha) or the band/k/npw layout differs.

  python3 -I compare_r1.py OUT/nscf.json qe-work/PREFIX.save/data-file-schema.xml [--json summary.json]
"""
import argparse
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def qe_bands(path):
    root = ET.parse(path).getroot()
    out = []
    for ks in root.iter("ks_energies"):
        npw = int(ks.find("npw").text)
        values = [float(x) for x in ks.find("eigenvalues").text.split()]
        out.append((npw, values))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("nscf_json")
    ap.add_argument("qe_xml")
    ap.add_argument("--gate", type=float, default=2e-5)
    ap.add_argument("--json")
    a = ap.parse_args()
    meta = json.loads(Path(a.nscf_json).read_text())
    ours = meta["band_data"]["eigenvalues_ha"][0]
    npw = [s["basis_size"] for s in meta["stats"][0]]
    qe = qe_bands(a.qe_xml)
    if len(qe) != len(ours):
        print(f"k count differs: BAANDI {len(ours)}, QE {len(qe)}")
        return 1
    bad_npw = [k for k, (n, _) in enumerate(qe) if n != npw[k]]
    occ = max(meta["occupied_bands"])
    offset = qe[0][1][occ - 1] - ours[0][occ - 1]
    worst_abs = worst_aligned = 0.0
    worst_at = None
    per_k = []
    for k, ((_, q), b) in enumerate(zip(qe, ours)):
        if len(q) < len(b):
            print(f"QE has fewer bands at k {k}: BAANDI {len(b)}, QE {len(q)}")
            return 1
        # QE solves extra top bands (make_qe_r1.py --extra-bands); compare the first len(b).
        q = q[: len(b)]
        d_abs = max(abs(x - y) for x, y in zip(q, b))
        d_al = [abs(x - offset - y) for x, y in zip(q, b)]
        m = max(d_al)
        per_k.append({"k_index": k, "npw": npw[k], "max_abs_ha": d_abs, "max_aligned_ha": m,
                      "worst_band": d_al.index(m) + 1})
        worst_abs = max(worst_abs, d_abs)
        if m > worst_aligned:
            worst_aligned, worst_at = m, (k, d_al.index(m) + 1)
    summary = {"bands": len(ours[0]), "k_points": len(ours), "npw_mismatch_k": bad_npw,
               "offset_ha": offset, "max_abs_ha": worst_abs, "max_aligned_ha": worst_aligned,
               "worst_k_band": worst_at, "gate_ha": a.gate, "per_k": per_k}
    print(json.dumps({k: v for k, v in summary.items() if k != "per_k"}, indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps(summary, indent=1) + "\n")
    return 0 if not bad_npw and worst_aligned <= a.gate else 1


if __name__ == "__main__":
    sys.exit(main())
