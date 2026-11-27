#!/usr/bin/env python3
"""Convert existing BAANDI relax inputs (one per seed, same settings) into one case JSON.

  python3 -I import_toml.py NAME OUT.json SEED=PATH [SEED=PATH ...]

Only the flat layout written by the exp-tabpfn generators is accepted
(lattice_unit 'bohr', [[system.atoms]] species/fractional, absolute UPF paths,
[basis]/[xc]/[scf]/[relax] tables, optional fixed=[{atom=i},...]). Every seed
must have identical cell, species order, UPF files and settings; the UPF paths are
stored relative to the common pseudo root (the parent of the first UPF's
directory chain up to 'pseudo'). Stdlib only.
"""
import json
import os
import re
import sys

ATOM = re.compile(r"\[\[system\.atoms\]\]\s*species\s*=\s*'([A-Za-z]+)'\s*fractional\s*=\s*\[([^\]]*)\]")


def parse(path):
    s = open(path).read()
    assert re.search(r"lattice_unit\s*=\s*'bohr'", s), path
    lat = json.loads(re.search(r"lattice_vectors\s*=\s*(\[\[.*?\]\])", s, re.S).group(1))
    atoms = [(m.group(1), [float(x) for x in m.group(2).split(",")]) for m in ATOM.finditer(s)]
    pp = dict(re.findall(r"(?m)^([A-Z][a-z]?)\s*=\s*'(/[^']+)'", s.split("[pseudopotentials]")[1].split("[basis]")[0]))
    sections = "[basis]" + s.split("[basis]", 1)[1].split("[relax]")[0]
    relax = s.split("[relax]", 1)[1] if "[relax]" in s else ""
    num = lambda k: float(re.search(r"(?m)^%s\s*=\s*([0-9.eE+-]+)" % k, relax).group(1))
    fixed = [int(x) for x in re.findall(r"\{\s*atom\s*=\s*(\d+)\s*\}", relax)]
    return {"lattice_bohr": lat, "atoms": atoms, "pp": pp, "toml_sections": sections.strip() + "\n",
            "relax": {"force_tol_ha_bohr": num("force_tol_ha_bohr"), "energy_tol_ha": num("energy_tol_ha"), "max_displacement_bohr": num("max_displacement_bohr"),
                      "trust_radius_bohr": num("trust_radius_bohr")},
            "fixed": fixed}


def main():
    name, out, pairs = sys.argv[1], sys.argv[2], [a.split("=", 1) for a in sys.argv[3:]]
    parsed = {seed: parse(path) for seed, path in pairs}
    first = parsed[pairs[0][0]]
    for seed, p in parsed.items():
        for key in ("lattice_bohr", "pp", "toml_sections", "relax", "fixed"):
            assert p[key] == first[key], (seed, key)
        # Species order may differ between seeds (different random defect placement).
        assert len(p["atoms"]) == len(first["atoms"]), seed
    any_pp = next(iter(first["pp"].values()))
    root = any_pp[:any_pp.index("/pseudo/") + len("/pseudo/")]
    case = {
        "name": name,
        "source": {seed: os.path.abspath(path) for seed, path in pairs},
        "lattice_bohr": first["lattice_bohr"],
        "pseudopotentials": {s: os.path.relpath(p, root) for s, p in first["pp"].items()},
        "toml_sections": first["toml_sections"],
        "relax": first["relax"],
        "fixed": first["fixed"],
        "starts": {seed: p["atoms"] for seed, p in parsed.items()},
    }
    with open(out, "w") as f:
        json.dump(case, f, indent=1)
    print(out, len(first["atoms"]), "atoms", len(parsed), "starts")


if __name__ == "__main__":
    main()
