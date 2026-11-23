#!/usr/bin/env python3
"""ISSUE-25: tables and density checks for the QE reference runs (python3 -I, stdlib only).

charge-density.dat (QE 7.6, no HDF5, Fortran unformatted, little-endian, 4-byte record markers):
  rec 1: gamma_only (int32 logical), ngm_g (int32), nspin (int32)
  rec 2: b1, b2, b3 reciprocal vectors (9 float64)
  rec 3: Miller indices, int32 (3, ngm_g)
  rec 4: rho(G), complex128 (ngm_g) per spin; gamma_only -> half sphere, G0 = N_e / Omega
"""
import os, re, sys, struct, array, math

D = sys.argv[1] if len(sys.argv) > 1 else "/data/sig-dft/i25"


def records(p):
    b = open(p, "rb").read(); o = 0; out = []
    while o < len(b):
        n = struct.unpack_from("<i", b, o)[0]; out.append(b[o + 4:o + 4 + n]); o += 8 + n
    return out


def rho(p):
    r = records(p); records_cache[p] = r
    gam, ngm, nspin = struct.unpack("<iii", r[0][:12])
    a = array.array("d"); a.frombytes(r[3])
    mill = array.array("i"); mill.frombytes(r[2])
    return dict(gamma=gam, ngm=ngm, nspin=nspin, rho=a, mill=mill, path=p)


def dr2(r0, rin, rout):
    """QE's printed 'estimated scf accuracy' for nspin=1, gamma_only (Hartree-like metric):
    dr2 = e2 * 4 pi * Omega/2 * 2 * sum_{G != 0} |rho_out(G) - rho_in(G)|^2 / |G|^2, e2 = 2 (Ry).
    b1..b3 in the file are reciprocal vectors in bohr^-1."""
    bv = struct.unpack("<9d", records_cache[r0["path"]][1])
    b = [bv[0:3], bv[3:6], bv[6:9]]
    det = (b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
           + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]))
    omega = (2 * math.pi) ** 3 / abs(det)
    m = r0["mill"]; s = 0.0
    for i in range(r0["ngm"]):
        h, k, l = m[3 * i], m[3 * i + 1], m[3 * i + 2]
        g2 = sum((h * b[0][j] + k * b[1][j] + l * b[2][j]) ** 2 for j in range(3))
        if g2 < 1e-12:
            continue
        dr = rout[2 * i] - rin[2 * i]; di = rout[2 * i + 1] - rin[2 * i + 1]
        s += (dr * dr + di * di) / g2
    return 2.0 * 4 * math.pi * omega * s


records_cache = {}


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def diff(a, b):
    return [x - y for x, y in zip(a, b)]


RX = dict(
    it=re.compile(r"iteration #\s*(\d+)"),
    ethr=re.compile(r"ethr =\s*([0-9.Ee+-]+),\s*avg # of iterations =\s*([0-9.]+)"),
    etot=re.compile(r"^\s*!?\s*total energy\s*=\s*([-0-9.]+) Ry", re.M),
    acc=re.compile(r"estimated scf accuracy\s*<\s*([0-9.Ee+-]+) Ry"),
)


def parse_stdout(p):
    txt = open(p, errors="replace").read()
    its = []
    blocks = re.split(r"(?=\n\s*iteration #)", txt)
    for b in blocks:
        m = RX["it"].search(b)
        if not m:
            continue
        e = RX["ethr"].findall(b)
        et = RX["etot"].search(b)
        ac = RX["acc"].search(b)
        its.append(dict(it=int(m.group(1)),
                        ethr=";".join(x[0] for x in e), avg=";".join(x[1] for x in e),
                        etot=et.group(1) if et else "", acc=ac.group(1) if ac else "",
                        relowered="Diagonalizing with lowered threshold" in b))
    conv = re.search(r"convergence (NOT achieved after\s+\d+ iterations: stopping|has been achieved in\s+\d+ iterations)", txt)
    final = re.search(r"^!\s*total energy\s*=\s*([-0-9.]+) Ry", txt, re.M)
    return its, (conv.group(0) if conv else "none"), (final.group(1) if final else "")


def cell_report(cell):
    C = os.path.join(D, cell)
    if not os.path.isdir(C):
        return
    print(f"\n## {cell}\n")
    # 1. per-k table: the last iteration of each run (= SCF step k)
    print("| k | mode | stop line | step-k estimated scf accuracy (Ry) | step-k ethr | step-k avg # Davidson its | step-k total energy (Ry) | '!' energy | wall s |")
    print("|---|---|---|---|---|---|---|---|---|")
    runs = {}
    for k in range(1, 7):
        for mode in ("mix", "out"):
            R = os.path.join(C, f"k{k}-{mode}")
            if not os.path.isdir(R):
                continue
            its, conv, final = parse_stdout(os.path.join(R, "stdout.txt"))
            runs[(k, mode)] = its
            last = its[-1] if its else {}
            wall = open(os.path.join(R, "time.txt")).read().split()[0] if os.path.exists(os.path.join(R, "time.txt")) else ""
            print(f"| {k} | {mode} | {conv} | {last.get('acc','')} | {last.get('ethr','')} | {last.get('avg','')} | {last.get('etot','')} | {final} | {wall} |")
    # 2. trajectory identity across runs (every run must reproduce the same first iterations)
    # reference trajectory = k6 mix run (an out run's last iteration prints the energy of rho_out, not the mixed one)
    ref = runs.get((6, "mix")) or []
    bad = 0
    for (k, mode), its in runs.items():
        for j, it in enumerate(its[:-1] if mode == "out" else its):
            if j < len(ref) and (it["acc"], it["etot"], it["ethr"]) != (ref[j]["acc"], ref[j]["etot"], ref[j]["ethr"]):
                bad += 1
    print(f"\nTrajectory prefix identical across all runs (acc, energy, ethr per iteration): {'yes' if bad == 0 else 'NO (%d mismatches)' % bad}")
    if any(it["relowered"] for its in runs.values() for it in its):
        print("NOTE: a run printed 'Diagonalizing with lowered threshold' (extra first-step diagonalisation, same iteration number).")
    print("\nFull trajectory (from the k=6 mix run):\n")
    print("| iteration | estimated scf accuracy (Ry) | ethr | avg # Davidson its | total energy (Ry) |")
    print("|---|---|---|---|---|")
    for it in ref:
        print(f"| {it['it']} | {it['acc']} | {it['ethr']} | {it['avg']} | {it['etot']} |")

    # 3. density checks
    def load(name):
        p = os.path.join(C, name, "charge-density.dat")
        return rho(p) if os.path.exists(p) else None
    out = {k: load(f"k{k}-out") for k in range(1, 7)}
    mix = {k: load(f"k{k}-mix") for k in range(1, 7)}
    m03 = load("k1-mix-beta0.3"); m10 = load("k1-mix-beta1.0")
    r0 = out[1] or mix[1]
    if not r0:
        print("\n(no density files)"); return
    print(f"\nDensity files: gamma_only={r0['gamma']} ngm_g={r0['ngm']} nspin={r0['nspin']} G0={r0['rho'][0]:.12g}")
    same_mill = all(d is None or d["mill"] == r0["mill"] for d in list(out.values()) + list(mix.values()) + [m03, m10])
    print(f"Miller index list identical in every file: {'yes' if same_mill else 'NO'}")
    print("\n### Relative differences ||a-b|| / ||b|| (2-norm over the stored half-sphere coefficients)\n")
    print("| check | value | meaning |")
    print("|---|---|---|")
    def rel(a, b):
        return norm(diff(a["rho"], b["rho"])) / norm(b["rho"])
    if out[1] and m10:
        print(f"| k1 mix beta=1.0 vs k1 out | {rel(m10, out[1]):.3e} | first-step mix at beta=1 equals rho_out(1) up to rounding |")
    if out[1] and mix[1] and m03:
        b07, b03, o = mix[1]["rho"], m03["rho"], out[1]["rho"]
        rin_a = [(x - 0.7 * y) / 0.3 for x, y in zip(b07, o)]
        rin_b = [(x - 0.3 * y) / 0.7 for x, y in zip(b03, o)]
        print(f"| rho_in(1) from beta 0.7 vs from beta 0.3 | {norm(diff(rin_a, rin_b)) / norm(rin_b):.3e} | reconstruction rho_in(1) = (rho_in(2) - beta rho_out(1))/(1-beta) is consistent |")
        print(f"| ||rho_out(1) - rho_in(1)|| / ||rho_in(1)|| | {norm(diff(o, rin_a)) / norm(rin_a):.3e} | size of the first update |")
        print(f"| rho_in(1) G0 | {rin_a[0]:.12g} | electron count / Omega |")
    for k in range(1, 7):
        if out[k] and mix[k]:
            print(f"| k{k}: mix(k) vs out(k) | {rel(mix[k], out[k]):.3e} | rho_in(k+1) vs rho_out(k) |")
    for k in range(2, 7):
        if out[k] and mix[k] and mix[k - 1]:
            rin = mix[k - 1]["rho"]
            lin = [a + 0.7 * (b - a) for a, b in zip(rin, out[k]["rho"])]
            dev = norm(diff(mix[k]["rho"], lin)) / norm(mix[k]["rho"])
            print(f"| k{k}: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | {dev:.3e} | 0 would mean plain linear mixing; nonzero = Broyden history |")
    print("\n### Recomputed 'estimated scf accuracy' from the stored files (decides what each file holds)\n")
    print("dr2(k) recomputed with QE's own metric from rho_in(k) and the k-out file; rho_in(1) = reconstruction above, rho_in(k>=2) = the (k-1)-mix file.")
    print("Agreement to every printed digit means: k-out holds rho_out(k) and (k-1)-mix holds rho_in(k), for every k.\n")
    print("| k | recomputed dr2 (Ry) | QE printed (Ry) | agree to printed digits |")
    print("|---|---|---|---|")
    for k in range(1, 7):
        if not out[k]:
            continue
        if k == 1:
            if not (mix[1] and out[1]):
                continue
            rin = [(x - 0.7 * y) / 0.3 for x, y in zip(mix[1]["rho"], out[1]["rho"])]
        else:
            if not mix[k - 1]:
                continue
            rin = mix[k - 1]["rho"]
        v = dr2(out[k], rin, out[k]["rho"])
        printed = runs.get((k, "out"), [{}])[-1].get("acc", "")
        ok = printed and ("%.8f" % v) == ("%.8f" % float(printed)) if printed else False
        print(f"| {k} | {v:.10g} | {printed} | {'yes' if ok else 'NO'} |")
    print("\n| check | value | meaning |")
    print("|---|---|---|")
    for k in range(2, 7):
        if out[k] and out[k - 1]:
            print(f"| k{k}: out(k) vs out(k-1) | {rel(out[k], out[k-1]):.3e} | successive rho_out |")


print("# ISSUE-25 QE reference: tables and density checks")
print(f"\nsource: {D}; generated by i25_qe_table.py")
for cell in ("si128", "si8", "si8-np1"):
    cell_report(cell)
