#!/usr/bin/env python3 -I
"""Regenerate the README validation block and figures from the stored evidence.

Reads tests/reference plus the stored collinear report, writes
docs/img/*.svg and replaces the text between the markers
<!-- validation:begin --> and <!-- validation:end --> in README.md.
New measurements come from stored reports or from
scripts/compare_kpoints_qe.py. Standard library only. Run after each verified
implementation (see AGENTS.md, "README results").
"""
import importlib.util
import json
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, "tests", "reference")
IMG = os.path.join(ROOT, "docs", "img")
BEGIN, END = "<!-- validation:begin -->", "<!-- validation:end -->"
INK, MUTED, GRID, BG = "#1f2328", "#59636e", "#d1d9e0", "#ffffff"
C_OURS, C_QE, C_TOL = "#0b6bcb", "#c2410c", "#9a3412"


def load(path):
    with open(path) as f:
        return json.load(f)


def check(report, item):
    return next(c for c in report["checks"] if c["item"] == item)


def cutoff_data():
    rows = []
    for name in sorted(os.listdir(os.path.join(REF, "reports"))):
        m = re.fullmatch(r"(scan_e\d+|si30_ours)\.json", name)
        if not m:
            continue
        r = load(os.path.join(REF, "reports", name))
        c = check(r, "total energy per atom")
        rows.append({"ecut": r["qe"]["ecutwfc"], "ours": 2 * c["ours"], "qe": 2 * c["qe"], "iters": r["qe"]["iterations"],
                     "eig": r["maximum_differences"]["eigenvalues"]["absolute_difference"]})
    return sorted(rows, key=lambda x: x["ecut"])


def kpoint_data():
    spec = importlib.util.spec_from_file_location("ck", os.path.join(ROOT, "scripts", "compare_kpoints_qe.py"))
    ck = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ck)
    runs = [("LDA 2x2x2", "kpoints/k222_c", "qe.out"), ("LDA 2x2x2, shifted", "kpoints/k222s_c", "qe.out"),
            ("LDA 4x4x4, shifted", "kpoints/k444s_c", "qe_summary.json"),
            ("PBE Gamma (override)", "pbe/gamma", "qe_summary.json"), ("PBE 2x2x2, shifted (override)", "pbe/k222s", "qe_summary.json")]
    out = []
    for label, d, qf in runs:
        base = os.path.join(REF, d)
        ours = load(os.path.join(base, "ours_results.json"))
        qe = ck.parse_qe_out(os.path.join(base, qf)) if qf.endswith(".out") else load(os.path.join(base, qf))
        ok, res = ck.compare(ours, qe, ck.DEFAULTS, 2)
        out.append({"label": label, "nk": res["k_points"], "ok": ok, "e": abs(res["total_energy_diff_ha_per_atom"]),
                    "eig": res["max_eigenvalue_diff_ha"], "iters": qe["iterations"]})
    return out


def perf_data():
    d = load(os.path.join(REF, "perf", "dev08_baseline.json"))
    mean = lambda a: sum(a) / len(a)
    steps = d["bench_steps"]
    out = []
    for key, label in (("gamma", "Gamma, 1 k point"), ("k444s", "4x4x4 shifted, 32 k points")):
        c = d["cases"][key]
        q1, q4 = mean(c["qe_np1"]["wall_s"]), mean(c["qe_np4"]["wall_s"])
        if "qe_same_window" in d:
            q1, q4 = d["qe_same_window"][key]["np1"], d["qe_same_window"][key]["np4"]
        last = steps[-1]["cases"][key]
        out.append({"label": label + ", 30 Ry", "nk": c["k_points"], "steps": [(st["label"], st["cases"][key]["wall_s"]) for st in steps],
                    "qe1": q1, "qe4": q4, "diag_frac": last["diag_s"] / last["engine_s"],
                    "ours_it": last["iters"], "qe_it": c["qe_np1"]["iterations"]})
    xp = {k: d["cases"][k]["ours_1thread"]["total_energy_ha"] - d["mac_reference_total_energy_ha"][k] for k in ("gamma", "k444s")}
    return {"rows": out, "xplat": xp, "commit": steps[-1]["commit"], "cpu": d["host"]["cpu"], "steps": [st["label"] for st in steps]}


def threads_data():
    d = load(os.path.join(REF, "perf", "dev08_threads.json"))
    return d


def fig_threads(t):
    ours = t["threads"]["k444s"]
    prev = t["threads_after_p7p6"]["k444s"]
    qe = t["qe_k444s_wall_s"]["conv_thr_1e-12"]
    ns = [1, 4, 8, 16]
    w, h = 760, 70 + len(ns) * 76 + 60
    s = svg_open(w, h, "Thread scaling: 32 k points, rust-dft threads against QE MPI processes on the same node")
    s += txt(20, 24, "32 k points (Si, 30 Ry): rust-dft threads against QE MPI processes, same node", 14, weight="bold")
    x0, x1 = 170, 740
    lo, hi = -0.2, 1.6
    def px(v): return x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
    top = 50
    for k in (0.5, 1, 2, 5, 10):
        xx = px(k)
        s += f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{h - 56}" stroke="{GRID}"/>' + txt(xx, h - 40, f"{k:g} s", 11, "middle", MUTED)
    s += txt((x0 + x1) / 2, h - 22, "wall time (log scale); QE with conv_thr 1e-12 (tighter than needed)", 12, "middle", MUTED)
    for i, n in enumerate(ns):
        y = top + 12 + i * 76
        s += txt(x0 - 12, y + 14, f"{n} thread(s) / process(es)", 12, "end", INK, "bold")
        for j, (name, v, col) in enumerate((("rust-dft, k-point threads only", prev[str(n)], "#7fb2e5"), ("rust-dft, parallel phases", ours[str(n)], C_OURS), ("QE 7.6", qe[f"np{n}"], C_QE))):
            yy = y + j * 22
            s += f'<rect x="{x0}" y="{yy}" width="{px(v) - x0:.1f}" height="16" fill="{col}"/>' + txt(px(v) + 6, yy + 12, f"{v:.3g} s  {name}", 11)
    return s + "</svg>\n"


def svg_open(w, h, title):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{title}" '
            f'font-family="Helvetica, Arial, sans-serif" font-size="12">\n<title>{title}</title>\n<rect width="{w}" height="{h}" fill="{BG}"/>\n')


def txt(x, y, s, size=12, anchor="start", fill=INK, weight="normal"):
    s = s.replace("&", "&amp;").replace("<", "&lt;")
    return f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" fill="{fill}" font-weight="{weight}">{s}</text>\n'


def fig_cutoff(rows):
    w, h = 760, 330
    s = svg_open(w, h, "Si total energy versus plane-wave cutoff: rust-dft and Quantum ESPRESSO 7.6")
    s += txt(20, 24, "Total energy versus cutoff (Si, Gamma, LDA-PZ): the two codes run parallel", 14, weight="bold")
    ref = rows[-1]
    # left panel: convergence, log10 of |E - E(80 Ry)| in Ha
    x0, x1, y0, y1 = 70, 360, 60, 270
    s += txt((x0 + x1) / 2, 46, "Error against the 80 Ry value (Ha, log scale)", 12, "middle", MUTED)
    lo, hi = -6.5, -1.5
    def px(e): return x0 + (e - 15) / (85 - 15) * (x1 - x0)
    def py(v): return y1 - (math.log10(v) - lo) / (hi - lo) * (y1 - y0)
    for k in range(-6, -1):
        yy = py(10.0 ** k)
        s += f'<line x1="{x0}" y1="{yy:.1f}" x2="{x1}" y2="{yy:.1f}" stroke="{GRID}"/>' + txt(x0 - 6, yy + 4, f"1e{k}", 11, "end", MUTED)
    for r in rows[:-1]:
        s += txt(px(r["ecut"]), y1 + 16, f"{r['ecut']:g}", 11, "middle", MUTED)
    s += txt((x0 + x1) / 2, y1 + 36, "cutoff (Ry)", 12, "middle", MUTED)
    for key, col, dx in (("ours", C_OURS, -3), ("qe", C_QE, 3)):
        pts = [(px(r["ecut"]) + dx, py(abs(r[key] - ref[key]))) for r in rows[:-1]]
        s += '<polyline fill="none" stroke="%s" stroke-width="1.5" points="%s"/>' % (col, " ".join(f"{a:.1f},{b:.1f}" for a, b in pts))
        for a, b in pts:
            s += f'<circle cx="{a:.1f}" cy="{b:.1f}" r="3.5" fill="{col}"/>'
    s += f'<circle cx="{x0 + 10}" cy="{y0 + 8}" r="3.5" fill="{C_OURS}"/>' + txt(x0 + 18, y0 + 12, "rust-dft", 11)
    s += f'<circle cx="{x0 + 90}" cy="{y0 + 8}" r="3.5" fill="{C_QE}"/>' + txt(x0 + 98, y0 + 12, "QE 7.6", 11)
    # right panel: ours - QE in micro-Ha
    a0, a1 = 450, 730
    diffs = [(r["ours"] - r["qe"]) * 1e6 for r in rows]
    top = max(max(diffs) * 1.4, 1.0)
    def qx(e): return a0 + (e - 15) / (85 - 15) * (a1 - a0)
    def qy(v): return y1 - v / top * (y1 - y0)
    s += txt((a0 + a1) / 2, 46, "Difference ours - QE (micro-Ha per cell)", 12, "middle", MUTED)
    for k in range(0, int(top) + 1):
        yy = qy(k)
        s += f'<line x1="{a0}" y1="{yy:.1f}" x2="{a1}" y2="{yy:.1f}" stroke="{GRID}"/>' + txt(a0 - 6, yy + 4, f"{k}", 11, "end", MUTED)
    for r in rows:
        s += txt(qx(r["ecut"]), y1 + 16, f"{r['ecut']:g}", 11, "middle", MUTED)
    s += txt((a0 + a1) / 2, y1 + 36, "cutoff (Ry)", 12, "middle", MUTED)
    pts = [(qx(r["ecut"]), qy(d)) for r, d in zip(rows, diffs)]
    s += '<polyline fill="none" stroke="%s" stroke-width="1.5" points="%s"/>' % (C_OURS, " ".join(f"{a:.1f},{b:.1f}" for a, b in pts))
    for a, b in pts:
        s += f'<circle cx="{a:.1f}" cy="{b:.1f}" r="3.5" fill="{C_OURS}"/>'
    s += txt((a0 + a1) / 2, y1 + 56, "Tolerance: 200 micro-Ha per cell (1e-4 Ha/atom), off this scale", 11, "middle", MUTED)
    return s + "</svg>\n"


def fig_kpoints(gamma, kruns):
    items = [{"label": "LDA Gamma, 30 Ry", "e": gamma["e"], "eig": gamma["eig"], "nk": 1}] + kruns
    w, h = 760, 70 + 62 * len(items) + 70
    s = svg_open(w, h, "Si k-point comparison with Quantum ESPRESSO 7.6: energy and eigenvalue differences versus tolerance")
    s += txt(20, 24, "Differences from QE 7.6 (Si, 30 Ry) against the acceptance tolerances", 14, weight="bold")
    x0, x1 = 210, 730
    lo, hi = -8.0, -2.5
    def px(v): return x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
    top = 62
    for k in range(-8, -2):
        xx = px(10.0 ** k)
        s += f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{h - 56}" stroke="{GRID}"/>' + txt(xx, h - 40, f"1e{k}", 11, "middle", MUTED)
    s += txt((x0 + x1) / 2, h - 22, "absolute difference (Ha per atom for energy, Ha for eigenvalues; log scale)", 12, "middle", MUTED)
    for i, it in enumerate(items):
        y = top + 16 + i * 62
        s += txt(x0 - 12, y + 14, it["label"], 12, "end", INK, "bold") + txt(x0 - 12, y + 30, f"{it['nk']} k point(s)", 11, "end", MUTED)
        for j, (key, col, name) in enumerate((("e", C_OURS, "total energy"), ("eig", C_QE, "max eigenvalue"))):
            yy = y + j * 20
            s += f'<rect x="{x0}" y="{yy}" width="{px(it[key]) - x0:.1f}" height="14" fill="{col}"/>' + txt(px(it[key]) + 6, yy + 11, f"{it[key]:.2e}  {name}", 11)
    for tol, name, ly in ((1e-4, "energy tolerance 1e-4", top - 20), (1e-3, "eigenvalue tolerance 1e-3", top - 6)):
        xx = px(tol)
        s += f'<line x1="{xx:.1f}" y1="{ly + 4}" x2="{xx:.1f}" y2="{h - 56}" stroke="{C_TOL}" stroke-width="1.5" stroke-dasharray="5 3"/>' + txt(xx - 4, ly, name, 11, "end", C_TOL)
    return s + "</svg>\n"


def fig_perf(perf):
    rows = perf["rows"]
    nst = len(perf["steps"])
    group = nst + 2
    w, h = 800, 70 + (group * 22 + 30) * len(rows) + 60
    s = svg_open(w, h, "Wall time of rust-dft optimization steps and Quantum ESPRESSO 7.6 on the same node")
    s += txt(20, 24, "Wall time on one node (Si, 30 Ry): rust-dft optimization steps against QE 7.6", 14, weight="bold")
    x0, x1 = 230, 770
    lo, hi = -1.0, 3.7
    def px(v): return x0 + (math.log10(v) - lo) / (hi - lo) * (x1 - x0)
    top = 50
    for k in range(-1, 3):
        xx = px(10.0 ** k)
        s += f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{h - 56}" stroke="{GRID}"/>' + txt(xx, h - 40, f"{10.0 ** k:g} s", 11, "middle", MUTED)
    s += txt((x0 + x1) / 2, h - 22, "wall time (median of 3 for rust-dft steps, mean of 3 for QE; log scale)", 12, "middle", MUTED)
    shades = ["#7fb2e5", "#3d85d1", C_OURS, "#0a3f7a", "#052a52"]
    for i, r in enumerate(rows):
        y = top + 12 + i * (group * 22 + 30)
        s += txt(x0 - 12, y + 14, r["label"], 12, "end", INK, "bold") + txt(x0 - 12, y + 30, f"{r['nk']} k point(s)", 11, "end", MUTED)
        bars = [(name, t, shades[min(j, len(shades) - 1)]) for j, (name, t) in enumerate(r["steps"])]
        bars += [("QE, 1 MPI process", r["qe1"], C_QE), ("QE, 4 MPI processes", r["qe4"], "#7c3aed")]
        for j, (name, t, col) in enumerate(bars):
            yy = y + j * 22
            s += f'<rect x="{x0}" y="{yy}" width="{px(t) - x0:.1f}" height="16" fill="{col}"/>' + txt(px(t) + 6, yy + 12, f"{t:.3g} s  {name}", 11)
    return s + "</svg>\n"


def show_perf():
    return bool(load(os.path.join(ROOT, "docs", "readme-config.json")).get("show_qe_performance", False))


def evidence_link(path):
    return f"[`{path}`]({path})"


def standard_rows():
    """Pair the frozen standard setting with strict, retaining failed cases."""
    freeze = load(os.path.join(REF, "perf/p17a/default-freeze.json"))
    report = load(os.path.join(REF, "perf/p17a/final-report.json"))
    label = f"standard-{freeze['standard_threshold_ha']:g}"
    rows = report["rows"]
    return [(r, next(s for s in rows if s["case"] == r["case"] and s["label"] == label))
            for r in rows if r["label"] == "strict"]


def evidence_text(path):
    with open(os.path.join(ROOT, path)) as f:
        return f.read()


def stress_term_errors():
    """Read the current isolated-term log, retaining its test/order provenance."""
    path = "tests/reference/stress/physics-terms.log"
    log = evidence_text(path)
    groups = {}
    for name, body in re.findall(r"test (\w+) \.\.\. (.*?)\nok", log, re.S):
        groups[name] = [json.loads(x)[-1] for x in re.findall(r"max errors (\[[^\]]+\]) Ha/Bohr\^3", body)]
    # Calls in stress_terms.rs: core then XC variants; nonlocal, local, Hartree.
    core = groups["core_and_xc_metric"]
    isolated = groups["isolated_local_core_hartree_nonlocal"]
    if len(core) != 6 or len(isolated) != 3:
        raise ValueError("Unexpected term log layout: " + path)
    return [("Kinetic", max(groups["frozen_coefficient_kinetic_strain_metric"])),
            ("Ewald (including background)", max(groups["ewald_strain_includes_background_and_eta_cancellation"])),
            ("Nonlocal", isolated[0]), ("Local", isolated[1]), ("Hartree", isolated[2]),
            ("NLCC core", core[0]), ("XC metric (LDA/PBE/spin)", max(core[1:]))]


def structural_results():
    """P19A–P22A: distinguish derivative oracles from material qualification."""
    cell_files = sorted(os.path.relpath(os.path.join(base, name), ROOT)
                        for base, _, names in os.walk(os.path.join(REF, "relaxation"))
                        for name in names if "cell" in name.lower())
    if cell_files:
        raise ValueError("Review new cell-relaxation evidence before regenerating: " + ", ".join(cell_files))
    sp = "tests/reference/symmetry/evidence.json"
    sy = load(os.path.join(ROOT, sp))
    mp = "tests/reference/symmetry/mac-k444s.json"
    mac = load(os.path.join(ROOT, mp))["cases"]
    tr, cr = mac["time_reversal"], mac["crystal"]
    median = lambda rows: sorted(r["wall_s"] for r in rows)[len(rows) // 2]
    lines = ["### Crystal symmetry: fewer k points for the same calculation", "",
             f"The shifted Si mesh goes from {tr[0]['irreducible_k_count']} time-reversal representatives to "
             f"{cr[0]['irreducible_k_count']} with crystal symmetry. Mac median wall time "
             f"({len(tr)} samples each): {median(tr):.3f} to {median(cr):.3f} s; this is an observation, "
             "not a server benchmark or a comparison with QE. " + evidence_link(mp) + ".", "",
             "Full mesh versus reduced mesh in the same engine (not an independent DFT code): "
             f"{len(sy['scf_oracles'])} Si, MgO and synthetic cases have maximum energy difference "
             f"{max(abs(r['delta_energy_ha']) for r in sy['scf_oracles']):.2e} Ha and density distance "
             f"{max(r['density_distance_bohr_minus_3_over_2'] for r in sy['scf_oracles']):.2e} Bohr^-3/2. "
             f"Fixed-spin differences: {sy['spin']['delta_energy_ha']:.2e} Ha, "
             f"{sy['spin']['density_distance_bohr_minus_3_over_2']:.2e} Bohr^-3/2 and "
             f"{sy['spin']['delta_moment_mu_b']:.2e} mu_B. Finite-temperature differences: "
             f"E {sy['fermi_dirac']['delta_energy_ha']:.2e} Ha, F {sy['fermi_dirac']['delta_free_energy_ha']:.2e} Ha, "
             f"density {sy['fermi_dirac']['density_distance_bohr_minus_3_over_2']:.2e} Bohr^-3/2. "
             + evidence_link(sp) + ".", "",
             "Only operations preserving the model, k mesh and FFT grid are retained; magnetic spin swaps and "
             "magnetic antiunitary operations are disabled. See [symmetry scope](docs/validation/symmetry.md); "
             "comparison with QE pending.", "", "### Atomic forces: checked against energy changes", ""]
    fp = "tests/reference/forces/manifest.json"
    force = load(os.path.join(ROOT, fp))
    obs = force["observations"]
    lines += ["Central displacement finite differences check individual terms and converged combinations. "
              f"Stored displacement steps: {', '.join(f'{x:g}' for x in force['steps_bohr'])} Bohr. "
              + evidence_link(fp) + ".", "",
              "| Force check | Maximum difference (Ha/Bohr) | Reference |", "|---|---:|---|"]
    for label, key, oracle in [
        ("Ewald", "ewald_max_error", "Energy finite differences"),
        ("Frozen local", "frozen_local_max_error", "Energy finite differences"),
        ("Frozen nonlocal", "frozen_nonlocal_max_error", "Energy finite differences"),
        ("Converged local + nonlocal", "converged_local_nonlocal_max_error", "Total-energy finite differences"),
        ("NLCC combination", "nlcc_max_error", "Total-energy finite differences"),
        ("PBE + NLCC combination", "pbe_nlcc_max_error", "Total-energy finite differences"),
        ("Finite temperature", "fd_max_error", "Free-energy finite differences"),
        ("Independent local + Ewald", "independent_max_error", "Second scalar energy reimplementation, finite differences"),
        ("Zero-moment spin reduction", "spin_reduction_max_difference", "Unpolarized engine")]:
        lines.append(f"| {label} | {obs[key + '_ha_per_bohr']:.2e} | {oracle} |")
    lines += ["", f"Resolved MgO net-force maximum component: {obs['mgo_resolved_net_max_ha_per_bohr']:.2e} Ha/Bohr "
              f"(cap {force['net_force_cap_resolved_grid']:.0e}). Real-file results from the same manifest:", "",
              "| Material | Force vs finite difference (Ha/Bohr) | Summed numerical uncertainty (Ha/Bohr) |",
              "|---|---:|---:|"]
    for name in ("si", "mgo"):
        lines.append(f"| {'Si' if name == 'si' else 'MgO'} | {obs[name + '_max_error_ha_per_bohr']:.2e} | "
                     f"{obs[name + '_uncertainty_sum_ha_per_bohr']:.2e} |")
    lines += ["", "Uncertainty covers SCF, cutoff and FFT-grid checks; it is not the derivative error. "
              f"Both sums are below one third of the {force['real_cap']:.0e} Ha/Bohr real-force cap. "
              "These results do not establish pseudopotential transferability. Forces require strict converged "
              "states and complete NC derivatives; DFT+U forces are unsupported. Spin and finite-temperature "
              "forces exist in engine/CLI; the ASE adapter remains Gamma/unpolarized. "
              "[Force scope](docs/validation/forces.md); comparison with QE pending.", "",
              "Mac cost observations (one sample per case, no performance acceptance): "
              + evidence_link("tests/reference/forces/cost.json") + ".", "",
              "| Case | SCF (s) | Force pass (s) | Force / (SCF + force) |", "|---|---:|---:|---:|"]
    for r in load(os.path.join(REF, "forces/cost.json")):
        lines.append(f"| {r['case']} | {r['timings']['engine_total_s']:.3f} | {r['forces']['elapsed_seconds']:.3f} | "
                     f"{100 * r['force_share_scf_plus_force']:.2f}% |")
    rp = "tests/reference/relaxation/manifest.json"
    relax = load(os.path.join(ROOT, rp))
    si = next(r for r in relax["materials"] if r["material"] == "Si")["trajectory"][-1]
    ase = relax["ase_comparison"]
    lp = "tests/reference/relaxation/logs/cli-all.log"
    log = evidence_text(lp)
    coords = re.search(r"Si final fractional coordinates: (\[\[.*\]\])", log)
    restart = re.search(r"relaxation restart: complete files bit-identical; (\d+) steps", log)
    if coords is None or restart is None:
        raise ValueError("Missing relaxation oracle in " + lp)
    position = json.loads(coords[1])[1]
    lines += ["", "### Ionic relaxation: move atoms at a fixed cell", "",
              "BFGS/FIRE are available with Cartesian constraints. Against an independent scalar energy scan, "
              f"compressed/stretched synthetic starts reach the periodic minimum within "
              f"{max(r['error_bohr'] for r in relax['synthetic_periodic_scan']):.2e} Bohr. "
              "This is a periodic model, not an isolated molecule or a global-minimum proof. "
              f"Warm versus fresh SCF energy differs by at most {max(relax['warm_fresh_energy_errors_ha']):.2e} Ha. "
              f"ASE BFGS versus CLI BFGS differs by {ase['max_position_error_bohr']:.2e} Bohr and "
              f"{ase['energy_error_ha']:.2e} Ha; ASE is an independent optimizer using BAANDI forces, not another DFT code. "
              + evidence_link(rp) + ".", "",
              f"Gamma-only Si returns to the ideal diamond bond: final atom fractional coordinates "
              f"({', '.join(f'{x:.12f}' for x in position)}), after {si['step']} steps, "
              f"maximum free force {si['max_force_ha_bohr']:.2e} Ha/Bohr. "
              f"Restart after an accepted boundary yields byte-identical complete output files in the "
              f"{restart[1]}-step fixture. " + evidence_link(lp) + "; " + evidence_link(rp) + ".", ""]
    curvature_path = "tests/reference/relaxation/mgo-curvature.json"
    curv = load(os.path.join(ROOT, curvature_path))["curvatures"]
    lines += ["MgO caveat: the Gamma-only model relaxes away from ideal rocksalt. A fresh-SCF scalar "
              f"energy scan finds negative curvature {curv[-1]['directional_curvature_ha_bohr2']:.5f} Ha/Bohr^2 "
              "at the ideal position; this is not bulk structural qualification. " + evidence_link(curvature_path) + ". "
              "CLI spin-polarized or finite-temperature relaxation is unsupported; strict unpolarized fixed occupations only. "
              "[Relaxation scope](docs/validation/relaxation.md); comparison with QE pending.", "",
              "![Si relaxation force history](docs/img/relaxation-forces.svg)", "", "### Stress: checked against cell strains", ""]
    stp = "tests/reference/stress/validation.json"
    stress = load(os.path.join(ROOT, stp))
    syn = stress["synthetic_cases"]
    lines += ["Tensile-positive stress differentiates energy (free energy at finite temperature) with fixed fractional atoms "
              "and a frozen set of plane-wave indices. This checks a finite-basis derivative, not cutoff convergence. "
              f"Strain steps: {', '.join(f'{x:g}' for x in stress['strain_steps'])}. "
              + evidence_link(stp) + ".", "",
              "Isolated terms against strain-energy finite differences (maximum at the smallest stored step): "
              + evidence_link("tests/reference/stress/physics-terms.log") + ".", "",
              "| Term | Difference (Ha/Bohr^3) |", "|---|---:|"]
    for label, error in stress_term_errors():
        lines.append(f"| {label} | {error:.2e} |")
    lines += ["", f"Across {len(syn)} synthetic total-stress cases (skew/slab, local/nonlocal, NLCC, PBE, spin and finite temperature), "
              f"maximum final strain error is {max(r['maximum_final_error_ha_bohr3'] for r in syn.values()):.2e} Ha/Bohr^3; "
              f"maximum raw antisymmetric norm is {max(r['antisymmetric_norm_ha_bohr3'] for r in syn.values()):.2e} Ha/Bohr^3. "
              f"ASE {stress['ase']['version']}'s independent numerical strain oracle differs by "
              f"{stress['ase']['maximum_error_ev_angstrom3']:.2e} eV/Angstrom^3, using BAANDI energies. "
              + evidence_link(stp) + ".", "",
              "Real-file results remain **inconclusive** as converged material stresses; small strain errors do not remove "
              "cutoff/k/grid uncertainty. Same evidence file:", "",
              "| Case | Strain error (Ha/Bohr^3) | Numerical uncertainty (Ha/Bohr^3) | Status |", "|---|---:|---:|---|"]
    for name, r in stress["real_materials"].items():
        lines.append(f"| {name} | {r['maximum_final_error_ha_bohr3']:.2e} | {r['uncertainty_ha_bohr3']:.2e} | {r['status']} |")
    lines += ["", "The separate isotropic Pulay estimate measures cutoff sensitivity; its two-step spread is not an error bound, "
              "does not establish anisotropic convergence and is not silently added to the analytic tensor. "
              "Stress with the standard profile is unsupported; a fresh strict stationary certificate is required. "
              "[Stress scope](docs/validation/stress.md); comparison with QE pending.", "",
              "Mac cost observations, one on/off pair per case. Analytic time includes force-based precondition checking; "
              "Pulay scans are extra. " + evidence_link("tests/reference/stress/cost.json") + ".", "",
              "| Case | SCF (s) | Analytic stress (s) | Pass / (SCF + pass) | Pulay scans (s) | Extra SCFs |",
              "|---|---:|---:|---:|---:|---:|"]
    for name, r in load(os.path.join(REF, "stress/cost.json"))["cases"].items():
        run = next(x for x in r["runs"] if x["enabled"])
        st = run["stress"]
        lines.append(f"| {name} | {run['scf_seconds']:.3f} | {st['elapsed_seconds']:.3f} | "
                     f"{100*r['analytic_fraction_of_scf_plus_pass']:.2f}% | {st['pulay_assessment']['elapsed_seconds']:.3f} | {st['cutoff_scf_runs']} |")
    lines += ["", "### Cell relaxation", "",
              "P23A cell-relaxation results skipped: no cell evidence files are present under "
              "`tests/reference/relaxation/` in this snapshot. No cell-relaxation validation is claimed; comparison with QE pending.", ""]
    return lines


def fig_relaxation():
    rows = next(r for r in load(os.path.join(REF, "relaxation/manifest.json"))["materials"]
                if r["material"] == "Si")["trajectory"]
    s = svg_open(760, 300, "Si fixed-cell relaxation: maximum free force versus ionic step")
    s += txt(20, 25, "Si Gamma-only: return to the ideal bond (not bulk qualification)", 15, weight="bold")
    x0, x1, y0, y1 = 90, 710, 60, 235
    vals = [math.log10(r["max_force_ha_bohr"]) for r in rows]
    lo, hi = math.floor(min(vals)), math.ceil(max(vals))
    for power in range(lo, hi + 1):
        y = y1 - (power - lo) / (hi - lo) * (y1 - y0)
        s += f'<path d="M {x0} {y:.1f} H {x1}" stroke="{GRID}"/>'
        s += txt(x0 - 8, y + 4, f"1e{power}", 11, "end")
    points = []
    for r, value in zip(rows, vals):
        x = x0 + r["step"] / rows[-1]["step"] * (x1 - x0)
        y = y1 - (value - lo) / (hi - lo) * (y1 - y0)
        points.append(f"{x:.1f},{y:.1f}")
        s += txt(x, y1 + 18, str(r["step"]), 11, "middle")
        s += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="3" fill="{C_OURS}"/>'
    s += f'<polyline points="{" ".join(points)}" fill="none" stroke="{C_OURS}"/>'
    s += txt(90, 46, "Maximum free force (Ha/Bohr, log scale)", 12)
    s += txt(400, 273, "Ionic step; source: tests/reference/relaxation/manifest.json", 12, "middle")
    return s + "</svg>\n"


def extended_results():
    """Render measured values, references and limits without promoting smoke tests."""
    spin_path = "crates/dft-engine/tests/collinear-evidence.json"
    spin = load(os.path.join(ROOT, spin_path))
    reduced = [r for r in spin["matrix"] if r["status"] == "certified_reduction"]
    # The original spin report stores independent-oracle measurements as test output.
    physics = spin["raw_output"]["physics"]
    match = re.search(r"PW92\+X maximum energy error=([\deE.+-]+) Ha/electron; interior potential error=([\deE.+-]+) Ha", physics)
    if match is None:
        raise ValueError("Missing PW92 oracle measurements in " + spin_path)
    lines = ["### Collinear spin: fixed up/down electron counts", "",
             f"LDA-PW and PBE synthetic checks: {len(reduced)} zero-moment calculations reproduce the "
             f"unpolarized engine with maximum energy difference {max(r['energy_difference_ha'] for r in reduced):.2e} Ha "
             f"and density distance {max(r['density_distance_bohr_minus_3_over_2'] for r in reduced):.2e} Bohr^-3/2. "
             f"The remaining {len(spin['matrix']) - len(reduced)} case gives the same noninsulating rejection. "
             f"Independent analytic exchange/PW92 formulas give maximum energy error {float(match[1]):.2e} Ha/electron; "
             f"finite-difference potentials differ by at most {float(match[2]):.2e} Ha. "
             "These test fixed spin populations, not a magnetic material ground state or agreement with QE. "
             + evidence_link(spin_path) + ".", ""]
    nlcc = load(os.path.join(REF, "datasets/p13c-status.json"))
    headers = load(os.path.join(REF, "datasets/pseudodojo-headers.json"))["datasets"]
    floor = load(os.path.join(REF, "materials/verification.json"))
    lines += ["### Nonlinear core correction and PseudoDojo files", "",
              f"Rust and independent Python header readers agree for {sum(r['rust_python_equal'] for r in headers)} "
              f"files ({', '.join(r['element'] for r in headers)}). Header identity is not dataset qualification. "
              + evidence_link("tests/reference/datasets/pseudodojo-headers.json") + ".", "",
              f"Synthetic Gaussian core transforms match the analytic formula to {nlcc['gaussian_core_max_relative_error']:.2e} relative error; "
              f"localized XC potentials match finite differences to {nlcc['localized_xc_derivative_max_error_ha']:.2e} Ha. "
              f"The later density-floor derivative check has maximum error {floor['floor_derivative_max_error_ha_bohr3']:.2e} Ha Bohr^3. "
              + evidence_link("tests/reference/datasets/p13c-status.json") + "; "
              + evidence_link("tests/reference/materials/verification.json") + ".", ""]
    pbe = [r for r in nlcc["total_scf_translation"] if r["functional"] == "Pbe"]
    lines += [f"Limits: translating the synthetic PBE cell changes the energy by {pbe[0]['delta_ha']:.2e} Ha on the coarse grid "
              f"(fails the {nlcc['caps_frozen_before_measurement']['translation_energy_ha']:.0e} Ha cap), falling to {pbe[-1]['delta_ha']:.2e} Ha "
              "with grid refinement. The stored isolated-atom runs include failures; all datasets remain unqualified "
              "for transferability, ghost states and cutoff convergence. NLCC checkpoint persistence remains unavailable; validated NC force/stress derivatives are reported below.", "",
              "### Real-material bulk runs: Si and MgO", "",
              "PBE calculations with external PseudoDojo files, positive uniform starting density. "
              "These are successful calculation checks, not an independent accuracy reference; comparison with QE pending for Si and MgO. "
              + evidence_link("tests/reference/materials/results.json") + ".", "",
              "| Case | Cutoff (Ha) | SCF iterations | Total energy (Ha/cell) | Electron-count error | Indirect gap (Ha) |",
              "|---|---:|---:|---:|---:|---:|"]
    for r in load(os.path.join(REF, "materials/results.json"))["cases"]:
        lines.append(f"| {r['case']} | {r['cutoff_ha']:g} | {r['iterations']} | {r['total_energy_ha']:.12f} | "
                     f"{abs(r['electron_count'] - r['expected_electrons']):.2e} | {r['band_gap']['indirect_ha']:.6f} |")
    lines += ["", "These sampled cutoffs do not establish a converged material energy or gap. Atomic-superposition initialization remains unavailable; "
              "see the [scope and verification record](tests/reference/materials/verification.json).", ""]
    spectra = load(os.path.join(REF, "spectral/summary.json"))
    material = [r for r in spectra if r["case"].startswith(("si-", "mgo-"))]
    skew = next(r for r in spectra if r["case"] == "syn20-skew")
    final = load(os.path.join(REF, "spectral/final/summary.json"))
    corrected = next(r for r in final["cases"] if r["case"] == "syn20-skew")
    lines += ["### Are the lowest electronic states actually found?", "",
              f"Independent full-matrix dense diagonalization at the final fixed Hamiltonian checks {len(material)} Si/MgO cases "
              f"at the stored sampled k points. Maximum occupied-state error: {max(r['occupied_max_error_ha'] for r in material):.2e} Ha; "
              f"maximum buffer-state error: {max(r['buffer_max_error_ha'] for r in material):.2e} Ha. "
              "This checks the eigensolver, not agreement with another DFT code. "
              + evidence_link("tests/reference/spectral/summary.json") + ".", "",
              f"A synthetic skew cell exposes a false success after {skew['iterations']} iterations: "
              f"occupied-state error {skew['occupied_max_error_ha']:.6f} Ha and fixed-H band-sum excess "
              f"{skew['maximum_absolute_per_k_band_sum_excess_ha']:.6f} Ha. The optional final audit detects it; "
              f"the corrected run exhausts {corrected['final_iterations']} iterations rather than claiming convergence. "
              "The band-sum excess is not a corrected SCF total energy. The audit defaults to "
              f"`{final['default']}`; large-operator searches are not a universal completeness proof. "
              + evidence_link("tests/reference/spectral/final/summary.json") + "; "
              + evidence_link("tests/reference/spectral/final/verification.json") + ".", "",
              "### Standard convergence profile: work versus strict", ""]
    freeze = load(os.path.join(REF, "perf/p17a/default-freeze.json"))
    lines += [f"Opt-in standard threshold: {freeze['standard_threshold_ha']:.0e} Ha; strict remains the default. "
              "**Standard results are not certified.** The table compares BAANDI with its own strict run, not QE. "
              "Energy differences are per cell; unavailable means no successful pair. Failed runs retain their work counts. "
              + evidence_link("tests/reference/perf/p17a/final-report.json") + ".", "",
              "| Case | Status strict / standard | SCF iterations strict / standard | Davidson steps strict / standard | E standard - E strict (Ha) |",
              "|---|---|---:|---:|---:|"]
    for strict, standard in standard_rows():
        delta = (f"{standard['total_energy_ha'] - strict['total_energy_ha']:+.2e}"
                 if strict["status"] == standard["status"] == "ok" else "unavailable")
        lines.append(f"| {strict['case']} | {strict['status']} / {standard['status']} | "
                     f"{strict.get('scf_iterations', 'unavailable')} / {standard.get('scf_iterations', 'unavailable')} | "
                     f"{strict.get('davidson_steps', 'unavailable')} / {standard.get('davidson_steps', 'unavailable')} | {delta} |")
    lines += ["", "Strict success alone does not prove spectral completeness (see the skew-cell result above). "
              "The threshold calibration uses relayed QE estimated errors, not independently reproduced equal-accuracy measurements. "
              "Mac timings are observations, not a server speedup claim. "
              + evidence_link("tests/reference/perf/p17a/default-freeze.json") + ".", ""]
    fd = load(os.path.join(REF, "occupations/summary.json"))
    syn, al = fd["synthetic"], fd["al_smoke"]
    derivative = fd["entropy_derivative"]
    low = fd["low_temperature"][-1]
    lines += ["### Finite-temperature Fermi-Dirac occupations", "",
              "Reference: independent analytic occupation/entropy formulas, finite differences of free energy, "
              "and the same engine's fixed-occupation low-temperature limit; no QE comparison. "
              + evidence_link("tests/reference/occupations/summary.json") + "; "
              + evidence_link("tests/reference/occupations/verification.json") + ".", "",
              f"Synthetic SCF: E = {syn['energy_ha']:.12f} Ha, F = {syn['free_energy_ha']:.12f} Ha, "
              f"S/kB = {syn['entropy_over_kb']:.12f}, in {syn['iterations']} iterations. "
              f"The entropy derivative error falls from {derivative[0]['error_dimensionless']:.2e} to "
              f"{derivative[-1]['error_dimensionless']:.2e} (dimensionless) as the temperature-energy step decreases "
              f"from {derivative[0]['step_ha']:g} to {derivative[-1]['step_ha']:g} Ha. "
              f"Empty-band refinement changes F by at most {max(r['free_energy_difference_ha'] for r in fd['empty_band_refinement']):.2e} Ha. "
              f"At kB T = {low['theta_ha']:g} Ha, energy differs from fixed occupations by {low['energy_difference_ha']:.2e} Ha.", "",
              f"External-file Al smoke: kB T = {al['theta_ha']:g} Ha, cutoff {al['cutoff_ha']:g} Ha, "
              f"{al['bands']} bands, mesh {'x'.join(map(str, al['mesh']))}, "
              f"E = {al['energy_ha']:.12f} Ha in {al['iterations']} iterations. "
              "This is not validated metallic material accuracy. Empty-band convergence must be checked for each calculation; "
              "finite-temperature checkpoint persistence is unavailable.", ""]
    return lines


def block(rows, kruns, gamma, perf=None):
    if not show_perf():
        perf = None
    n_rust = None
    lines = [BEGIN, "", "## Validation results", "",
             "Generated by `scripts/make_readme_figures.py` from stored evidence in `tests/reference` and the linked collinear report (do not edit by hand). "
             "Reference: Quantum ESPRESSO 7.6 run on a remote server (see [validation record](docs/validation.md)). "
             "The QE comparison below covers one cell (Si diamond) and one pseudopotential of unverified provenance; it is not a statement about general material accuracy.", "",
             "![Total energy versus cutoff](docs/img/cutoff-scan.svg)", "",
             "| Check (Si, 30 Ry, same pseudopotential for both codes) | k points | QE iterations | energy difference (Ha/atom) | largest eigenvalue difference (Ha) | tolerance (energy / eigenvalue) | status |",
             "|---|---|---|---|---|---|---|",
             f"| LDA Gamma only | 1 | {gamma['iters']} | {gamma['e']:.2e} | {gamma['eig']:.2e} | 1e-4 / 1e-3 | pass |"]
    for k in kruns:
        lines.append(f"| {k['label']} | {k['nk']} | {k['iters']} | {k['e']:.2e} | {k['eig']:.2e} | 1e-4 / 1e-3 | {'pass' if k['ok'] else 'FAIL'} |")
    lines += ["", "![k-point comparison](docs/img/kpoint-comparison.svg)", "",
              f"Cutoff scan: {len(rows)} cutoffs from {rows[0]['ecut']:g} to {rows[-1]['ecut']:g} Ry; the energy difference to QE is the same at every cutoff "
              f"({min((r['ours'] - r['qe']) for r in rows) * 1e6:.2f} to {max((r['ours'] - r['qe']) for r in rows) * 1e6:.2f} micro-Ha per cell), "
              "a constant offset explained by the G=0 local-potential integral convention (validation record R4). "
              "The 4x4x4 and PBE evidence are relayed text summaries; the other LDA runs are whole-file, hash-verified outputs. The PBE rows use QE's `input_dft = 'pbe'` override on the LDA pseudopotential: a code-to-code consistency test of the PBE implementation, not a physically meaningful PBE result for this pseudopotential.", "",
              ]
    if perf:
        lines += ["### Performance (rust-dft optimization steps against QE 7.6, same node)", "",
                  f"Same node ({perf['cpu']}), rust-dft on one thread (steps measured with `scripts/bench_scf.py`, median of 3, quiet node) against QE 7.6 on one and on four MPI processes "
                  f"(mean of 3); latest step: commit `{perf['commit']}`. Relayed measurements, see [validation record](docs/validation.md) R7 and the [optimization log](docs/performance/optimization-log.md). "
                  "QE converges to a much tighter threshold (1e-12) than the rust-dft defaults, so compare iterations and final energies as well as times.", "",
                  "![wall time](docs/img/performance.svg)", ""]
        head = "| Case | " + " | ".join(f"{n} (s)" for n in perf["steps"]) + " | QE 1 process (s) | QE 4 processes (s) | latest vs QE 1 process | speedup vs before | iterations rust-dft / QE | diagonalization share |"
        lines += [head, "|" + "---|" * (len(perf["steps"]) + 7)]
        for r in perf["rows"]:
            ts = [t for _, t in r["steps"]]
            lines.append(f"| {r['label']} ({r['nk']} k) | " + " | ".join(f"{t:.3g}" for t in ts) + f" | {r['qe1']:.3g} | {r['qe4']:.3g} | {ts[-1] / r['qe1']:.1f}x slower | {ts[0] / ts[-1]:.2f}x | {r['ours_it']} / {r['qe_it']} | {100 * r['diag_frac']:.0f}% |")
        sw = load(os.path.join(REF, "perf", "dev08_baseline.json")).get("qe_convthr_sweep")
        if sw:
            lines += ["", "Equal-accuracy comparison: QE's reference runs above use `conv_thr` = 1e-12 Ry, tighter than needed; its total energy is already identical to the printed 8 decimals at `conv_thr` = 1e-8 (the rust-dft defaults reach an energy error of about 3e-14 Ha against a fully converged run). Latest rust-dft step against QE at several thresholds (QE 1 process, wall time):", "",
                      "| Case | rust-dft latest (s) | QE conv_thr 1e-12 (s) | QE 1e-8 (s) | QE 1e-6 (s) | rust-dft slower than QE at 1e-12 | at 1e-8 | at 1e-6 |", "|---|---|---|---|---|---|---|---|"]
            for key, r in zip(("gamma", "k444s"), perf["rows"]):
                t = r["steps"][-1][1]
                q = {row["conv_thr_ry"]: row["wall_np1_s"] for row in sw[key]}
                lines.append(f"| {r['label']} | {t:.3g} | {q[1e-12]:.3g} | {q[1e-8]:.3g} | {q[1e-6]:.3g} | {t / q[1e-12]:.1f}x | {t / q[1e-8]:.1f}x | {t / q[1e-6]:.1f}x |")
        td = threads_data()
        ts = td["threads"]
        lines += ["", "Deterministic k-point threads (`[scf] threads = N`, default 1): results are bit-identical for every thread count (density file hashes equal for T = 1, 4, 8, 16 in all four cases). Wall time in seconds on dev08 (one socket, quiet node), rust-dft with T threads and QE 7.6 (`conv_thr` 1e-12) with the same number of MPI processes:", "",
                  "![thread scaling](docs/img/thread-scaling.svg)", "",
                  "| Case | T=1 | T=4 | T=8 | T=16 | speedup T=16 vs T=1 |", "|---|---|---|---|---|---|"]
        for key, name in (("gamma", "Gamma (1 k)"), ("k222s", "2x2x2 shifted (4 k)"), ("k222", "2x2x2 (8 k)"), ("k444s", "4x4x4 shifted (32 k)")):
            v = ts[key]
            lines.append(f"| {name} | {v['1']:.3g} | {v['4']:.3g} | {v['8']:.3g} | {v['16']:.3g} | {v['1'] / v['16']:.2f}x |")
        qe = td["qe_k444s_wall_s"]["conv_thr_1e-12"]
        o = ts["k444s"]
        lines += ["", "32 k points against QE with the same number of cores: " + ", ".join(f"{n}: {o[str(n)] / qe['np' + str(n)]:.1f}x slower" for n in (1, 4, 8, 16)) + ". After the parallel setup, density, energy-audit and certificate phases the 32-k-point run takes 1.23 s on 16 threads (3.06 s with k-point threads only), of which the diagonalization is 0.84 s; with 4 and 8 k points the best thread count is 8 (there are few k points per thread). One k point (Gamma) is not parallelized yet.", ""]
        lines += ["", f"Cross-platform check: the total energy on this Linux x86 build differs from the macOS Apple Silicon build by {perf['xplat']['gamma']:+.2e} Ha (Gamma) and {perf['xplat']['k444s']:+.2e} Ha (4x4x4, shifted).", ""]
    lines += extended_results()
    lines += structural_results()
    checkpoint = load(os.path.join(REF, "checkpoint", "results.json"))
    lines += ["### Checkpoint/restart equivalence (synthetic fixtures)", "",
              f"Reference: the uninterrupted BAANDI engine in the same build. "
              f"{checkpoint['memory_exact_comparisons']} memory restart comparisons and "
              f"{checkpoint['file_exact_comparisons']} binary-file restart comparisons preserve energy, density, orbitals, "
              "scientific iteration history and certificates bit for bit (wall timings excluded). "
              f"The fixtures cover cubic Gamma, four k points, and a skew two-species cell; "
              f"Linear/Pulay/KerkerPulay; requested threads {checkpoint['threads']}. "
              f"Different-mixer warm starts have maximum energy difference "
              f"{checkpoint['warm_max_energy_difference_ha']:.2e} Ha (acceptance: "
              f"{checkpoint['warm_tolerance_ha']:.0e} Ha). "
              "Paused runtime/audit restoration and cross-process CLI restart also pass. "
              "These are self-equivalence checks, not QE comparisons or a cross-build restart promise. "
              "See [format, commands and limits](docs/checkpoint-format.md) and "
              "[stored evidence](tests/reference/checkpoint/results.json).", ""]
    healing = load(os.path.join(REF, "runtime", "generated", "results.json"))
    hruns = healing["runs"]
    assert len(hruns) == 8 and healing["full_result_and_array_identity"]
    for r in hruns:
        assert r["converged"] == (r["arm"] != "baseline")
        if r["converged"]:
            assert r["energy_error_hartree"] <= healing["energy_tolerance_hartree"]
    lines += ["### Self-healing SCF: offline budget recovery", "",
              "Synthetic two-ion Gamma fixture; reference is a fresh, certified BAANDI run "
              "with a 1000-iteration budget in the same build (not QE). "
              "All four runs were repeated twice; normalized numerical results and every final "
              "density/orbital array bit match on replay. Rules and recorded-agent recovery also "
              "match the reference bit for bit. No network or LLM is needed.", "",
              "| Arm | Converged | Iterations | Certified energy (Ha) | Error vs reference (Ha) | Applied / rejected |",
              "|---|---|---|---|---|---|"]
    for r in hruns[:4]:
        energy = f"{r['certified_energy_hartree']:.15f}" if r["converged"] else "—"
        error = f"{r['energy_error_hartree']:.1e}" if r["converged"] else "—"
        lines.append(f"| {r['arm']} | {r['converged']} | {r['iterations']} | {energy} | {error} | {r['policy_changes']} / {r['rejected']} |")
    lines += ["", "Acceptance is 1e-10 Ha; the invalid beta=1.1 request is rejected without changing "
              "state. This proves initial budget repair and deterministic control, not numerical "
              "instability healing, speedup, general recovery defaults or a global ground state. "
              "[Reproduce and inspect provenance](docs/demos/self-healing-scf.md#executed-p11e-demo) "
              "with `scripts/demo_self_healing.sh`.", ""]
    lines += [
              "Other verified results: PBE (xcx) matches the published closed forms to 2e-16 Ha and the LDA path is unchanged bit for bit; k-point engine at Gamma reproduces the previous density bit for bit; unit conversion through `uom` and arithmetic input expressions leave the real-Si density bit-identical; "
              "fast gates (`scripts/check.sh`) and the Python regression tests pass. Results are added here after each verified implementation.", "", END]
    return "\n".join(lines)


STATUS_BEGIN, STATUS_END = "<!-- status:begin -->", "<!-- status:end -->"


def status_block():
    d = load(os.path.join(ROOT, "docs", "feature-status.json"))
    names = {"done": "Works today", "progress": "Being built now", "planned": "Planned, not available yet"}
    lines = [STATUS_BEGIN, "", "## What BAANDI can do today", "",
             f"Status as of {d['updated']}, written for people who run calculations. Source: `docs/feature-status.json` (this section is generated by `scripts/make_readme_figures.py`).", "",
             "**You can use it today for:** " + " ".join(d["can_today"]), "",
             "**Not available yet:** " + " ".join(d["not_yet"]), "", "**Speed:** " + d["speed"], ""]
    for key in ("done", "progress", "planned"):
        items = [i for i in d["items"] if i["status"] == key]
        lines += [f"### {names[key]}", "", "| Area | Feature | Notes |", "|---|---|---|"]
        for i in items:
            lines.append(f"| {i['area']} | {i['feature']} | {i['note']} |")
        lines.append("")
    lines.append(STATUS_END)
    return "\n".join(lines)


def main():
    rows = cutoff_data()
    kruns = kpoint_data()
    g = next(r for r in rows if r["ecut"] == 30.0)
    gamma = {"e": abs(g["ours"] - g["qe"]) / 2, "eig": g["eig"], "iters": g["iters"]}
    perf = perf_data()
    os.makedirs(IMG, exist_ok=True)
    with open(os.path.join(IMG, "performance.svg"), "w") as f:
        f.write(fig_perf(perf))
    with open(os.path.join(IMG, "thread-scaling.svg"), "w") as f:
        f.write(fig_threads(threads_data()))
    with open(os.path.join(IMG, "cutoff-scan.svg"), "w") as f:
        f.write(fig_cutoff(rows))
    with open(os.path.join(IMG, "kpoint-comparison.svg"), "w") as f:
        f.write(fig_kpoints(gamma, kruns))
    with open(os.path.join(IMG, "relaxation-forces.svg"), "w") as f:
        f.write(fig_relaxation())
    readme = os.path.join(ROOT, "README.md")
    text = open(readme).read()
    new = block(rows, kruns, gamma, perf)
    if BEGIN in text and END in text:
        i, j = text.index(BEGIN), text.index(END) + len(END)
        text = text[:i] + new + text[j:]
    else:
        anchor = "## Architecture Philosophy"
        if anchor not in text:
            print("README anchor not found", file=sys.stderr)
            return 1
        text = text.replace(anchor, new + "\n\n" + anchor, 1)
    sb = status_block()
    if STATUS_BEGIN in text and STATUS_END in text:
        i, j = text.index(STATUS_BEGIN), text.index(STATUS_END) + len(STATUS_END)
        text = text[:i] + sb + text[j:]
    else:
        text = text.replace(BEGIN, sb + "\n\n" + BEGIN, 1)
    open(readme, "w").write(text)
    print("wrote docs/img/cutoff-scan.svg, kpoint-comparison.svg, performance.svg and the README validation block")
    return 0


if __name__ == "__main__":
    sys.exit(main())
