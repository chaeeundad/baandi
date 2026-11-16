"""F11: run si16/si32/si64 (one atom displaced 1e-3 fractional, given profile) and
dump the GEMM dispatch rows with bands, plane waves and capacity.
Usage: python3 -I scripts/gemm_dispatch_survey.py BINARY CASE_DIR UPF OUTDIR [profile]
See docs/performance/gemm-dispatch-analysis-2026-11-16.md (MKL orthogonalization
defaults: set BAANDI_MGS=blocked BAANDI_ORTHO_BLOCK=off)."""
import json, os, subprocess, sys

binary, case_dir, upf, outdir = sys.argv[1:5]
profile = sys.argv[5] if len(sys.argv) > 5 else "standard"
for cid in ("si16-gamma", "si32-gamma", "si64-gamma"):
    case = {}
    for line in open(os.path.join(case_dir, cid + ".toml")):
        line = line.strip()
        if line and not line.startswith("#"):
            k, v = line.split("=", 1)
            case[k.strip()] = json.loads(v.strip())
    atoms = [list(a) for a in case["atoms"]]
    atoms[-1][2] += 1e-3
    lines = ["[system]", f'name = "{cid}"', 'lattice_unit = "bohr"',
             f'lattice_vectors = {json.dumps(case["lattice"])}']
    for a in atoms:
        lines += ["[[system.atoms]]", 'species = "Si"', f"fractional = {json.dumps(a)}"]
    lines += ["[pseudopotentials]", f"Si = {json.dumps(upf)}", "[basis]",
              f'ecutwfc_ry = {case["cutoff_ry"]}', "kpoints = [[0.0, 0.0, 0.0]]",
              "[xc]", 'functional = "lda-pz"', "[scf]", "threads = 8", "max_iterations = 100",
              f'convergence_profile = "{profile}"']
    src = os.path.join(outdir, f"{cid}-{profile}.toml")
    open(src, "w").write("\n".join(lines) + "\n")
    out = os.path.join(outdir, f"{cid}-{profile}")
    p = subprocess.run([binary, "run", src, "--output-dir", out], capture_output=True, text=True)
    if p.returncode != 0:
        print(cid, "exit", p.returncode, p.stderr[-400:])
        continue
    r = json.load(open(os.path.join(out, "results.json")))
    g = r["convergence"]["work"]["gemm_dispatch"]
    print(json.dumps({"case": cid, "profile": profile, "bands": r["bands"],
                      "n_plane_waves": r["basis"]["n_plane_waves"],
                      "capacity_per_band": r["scf"].get("davidson_capacity_per_band"),
                      "scf": len(r["convergence"]["history"]),
                      "wide_requests": g["wide_requests"]}))
    for row in g["rows"]:
        print("   ", json.dumps(row))
    print("  as_if_mkl:")
    for row in g.get("rows_as_if_mkl", []):
        print("   ", json.dumps(row))
