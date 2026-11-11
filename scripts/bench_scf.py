#!/usr/bin/env python3 -I
"""Benchmark and equivalence check for the rust-dft SCF (standard library only).

Runs the release dft-cli on fixed Si cases, records wall time (median of N runs),
engine phase times, iterations and total energy, and compares them with a stored
baseline: an optimization must keep the total energy within 1e-12 Ha and, where
the change is meant to be numerically neutral, the density file bit-identical.
The pseudopotential is not in the repository: pass --upf (never commit it).

  python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /path/Si.pz-vbc.UPF \
      --cases gamma,k222 --repeats 3 --baseline docs/performance/baseline.json
  ... --record docs/performance/baseline.json     # write a new baseline (initial or after a decision)
Exit codes: 0 all cases within tolerance, 5 energy outside tolerance or density not
identical when --require-identical, 3 inconsistent inputs, 1 usage or I/O error.
"""
import argparse, hashlib, json, os, statistics, subprocess, sys, tempfile, time

TEMPLATE = """[system]
name = "Si diamond primitive, a = 5.43 Angstrom"
lattice_unit = "angstrom"
lattice_vectors = [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]

[[system.atoms]]
species = "Si"
fractional = [0.0, 0.0, 0.0]

[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.25]

[pseudopotentials]
Si = "{upf}"

[basis]
ecutwfc_ry = {ecut}
{kblock}

[xc]
functional = "lda-pz"

[scf]
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
"""
CASES = {
    "si8": ("8-atom conventional Gamma, 30 Ry", 30.0, "kpoints = [[0.0, 0.0, 0.0]]"),
    "gamma": ("Gamma only, 30 Ry", 30.0, "kpoints = [[0.0, 0.0, 0.0]]"),
    "k222s": ("2x2x2 shifted (4 k), 30 Ry", 30.0, "kgrid = [2, 2, 2]\nkshift = [0.5, 0.5, 0.5]"),
    "k222": ("2x2x2 (8 k), 30 Ry", 30.0, "kgrid = [2, 2, 2]\nkshift = [0.0, 0.0, 0.0]"),
    "k444s": ("4x4x4 shifted (32 k), 30 Ry", 30.0, "kgrid = [4, 4, 4]\nkshift = [0.5, 0.5, 0.5]"),
}
ENERGY_TOL_HA = 1e-12


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def case_input(upf, key, threads=None, extra_toml=()):
    _, ecut, kblock = CASES[key]
    text = TEMPLATE.format(upf=upf, ecut=ecut, kblock=kblock)
    if key == "si8":
        text = text.replace('Si diamond primitive', 'Si diamond conventional')
        text = text.replace('[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]',
                            '[[5.43, 0.0, 0.0], [0.0, 5.43, 0.0], [0.0, 0.0, 5.43]]')
        start = text.index('[[system.atoms]]')
        end = text.index('[pseudopotentials]')
        atoms = ""
        for site in [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0),
                     (.25, .25, .25), (.25, .75, .75), (.75, .25, .75), (.75, .75, .25)]:
            atoms += f'[[system.atoms]]\nspecies = "Si"\nfractional = {list(site)}\n\n'
        text = text[:start] + atoms + text[end:]
    if threads is not None:
        text = text.replace("[scf]", f"[scf]\nthreads = {threads}")
    return text + "\n".join(extra_toml) + ("\n" if extra_toml else "")


def run_case(binary, upf, key, repeats, threads=None, extra_toml=()):
    label, ecut, kblock = CASES[key]
    walls, phases, energies, hashes, engines, last = [], [], [], [], [], None
    with tempfile.TemporaryDirectory() as td:
        toml = os.path.join(td, "input.toml")
        text = case_input(upf, key, threads, extra_toml)
        open(toml, "w").write(text)
        for i in range(repeats):
            out = os.path.join(td, f"out{i}")
            t0 = time.perf_counter()
            p = subprocess.run([binary, "run", toml, "--output-dir", out], capture_output=True, text=True)
            walls.append(time.perf_counter() - t0)
            if p.returncode != 0:
                raise RuntimeError(f"{key}: dft-cli exit {p.returncode}\n{p.stderr[-400:]}")
            last = out
            report = json.load(open(os.path.join(out, "results.json")))
            phases.append(report["timings"]["phases"])
            engines.append(report["timings"]["engine_total_s"])
            energies.append(report["total_energy"]["ha"])
            hashes.append(sha256(os.path.join(out, "density.bin")))
        if len(set(energies)) != 1 or len(set(hashes)) != 1:
            raise RuntimeError(f"{key}: repeated runs were not identical")
        r = json.load(open(os.path.join(last, "results.json")))
        ph = r["timings"]["phases"]
        history = r["convergence"]["history"]
        certificates = [h["fixed_point_certificate"] for h in history if h.get("fixed_point_certificate")]
        # Switch manifest (scf.switches): every BAANDI_* variable set for this run
        # must be a known switch that the binary accepted from the environment.
        switches = r.get("scf", {}).get("switches", {})
        for name in sorted(k for k in os.environ if k.startswith("BAANDI_")):
            entry = switches.get(name)
            if switches and (entry is None or entry["source"] != "env"):
                raise RuntimeError(f"{key}: {name} is set but results.json reports {entry}")

        return {"operator_applications": sum(h["operator_applications"] for h in history),
                "davidson_iterations": sum(h["davidson_iterations"] for h in history),
                "certificate_applications": sum(h["operator_applications"] for h in certificates),
                "certificate_iterations": sum(h["davidson_iterations"] for h in certificates),
                "label": label, "wall_s": walls, "wall_median_s": statistics.median(walls),
                "engine_total_s": r["timings"]["engine_total_s"], "phases_s": ph,
                "engine_total_median_s": statistics.median(engines),
                "phase_medians_s": {key: statistics.median(p[key] for p in phases) for key in ph},
                "iterations": len(r["convergence"]["history"]), "total_energy_ha": r["total_energy"]["ha"],
                "density_sha256": sha256(os.path.join(last, "density.bin")),
                "switches_not_default": {n: e["value"] for n, e in switches.items() if e["source"] != "default"}}


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--upf", required=True)
    ap.add_argument("--cases", default="gamma,k222")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--threads", type=int, help="total thread budget (0 means 1); omitted preserves the input default")
    ap.add_argument("--extra-toml", action="append", default=[], help="append TOML lines (initially in [scf]); repeatable")
    ap.add_argument("--baseline")
    ap.add_argument("--record")
    ap.add_argument("--require-identical", action="store_true", help="fail unless the density is bit-identical to the baseline")
    ap.add_argument("--skip-legacy-si8-anchor", action="store_true",
                    help="skip the historical strict/complex Si8 absolute-energy anchor for opt-in profiles; baseline energy/density checks remain enabled")
    ap.add_argument("--note", default="")
    a = ap.parse_args(argv)
    if a.threads is not None and a.threads < 0:
        ap.error("--threads must be nonnegative")
    keys = a.cases.split(",")
    if any(k not in CASES for k in keys) or a.repeats < 1:
        print("usage error: unknown case or repeats < 1")
        return 1
    if not (os.path.isfile(a.binary) and os.path.isfile(a.upf)):
        print("I/O error: binary or UPF not found")
        return 1
    base = json.load(open(a.baseline))["cases"] if a.baseline and os.path.isfile(a.baseline) else {}
    result, status = {}, 0
    print(f"{'case':8} {'wall med (s)':>13} {'engine (s)':>11} {'diag (s)':>9} {'iters':>5} {'dE vs base (Ha)':>16} {'density':>10} {'speedup':>8}")
    for k in keys:
        try:
            r = run_case(a.binary, a.upf, k, a.repeats, a.threads, a.extra_toml)
        except RuntimeError as e:
            print(e)
            return 3
        result[k] = r
        if not a.skip_legacy_si8_anchor and k == "si8" and abs(r["total_energy_ha"] - (-31.3218036010316)) > 1e-9:
            status = 5
        b = base.get(k)
        de = abs(r["total_energy_ha"] - b["total_energy_ha"]) if b else None
        ident = (r["density_sha256"] == b["density_sha256"]) if b else None
        sp = b["wall_median_s"] / r["wall_median_s"] if b else None
        if b and (de > ENERGY_TOL_HA or (a.require_identical and not ident)):
            status = 5
        print(f"{k:8} {r['wall_median_s']:13.3f} {r['engine_total_s']:11.3f} {r['phases_s']['diagonalization_s']:9.3f} {r['iterations']:5d} "
              f"{('%.2e' % de) if de is not None else '-':>16} {('identical' if ident else 'differs') if b else '-':>10} {('%.2fx' % sp) if sp else '-':>8}")
    if a.record:
        meta = {"note": a.note, "binary": os.path.basename(a.binary), "machine_note": "record host and commit in the note",
                "energy_tolerance_ha": ENERGY_TOL_HA, "cases": result}
        json.dump(meta, open(a.record, "w"), indent=1)
        open(a.record, "a").write("\n")
        print("recorded", a.record)
    print("OK" if status == 0 else "FAIL: energy outside 1e-12 Ha or density not identical")
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
