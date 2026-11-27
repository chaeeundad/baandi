"""Benchmark cases: a starting structure plus the fixed electronic settings.

A case file (JSON, ``cases/*.json``) holds the cell (rows = lattice vectors,
Bohr), the species/fractional coordinates of the unperturbed or pre-perturbed
structures, the pseudopotential file names (relative to a pseudo root), the
verbatim ``[basis]``/``[xc]``/``[scf]`` TOML text, the ionic step caps and the
fixed atoms. Seeds select either a pre-perturbed start (``starts``) or a
seeded perturbation of ``atoms`` (``perturbation``), never both. Stdlib only.
"""
import hashlib
import json
import math
import os
import random
import re

BOHR_A = 0.529177210903  # Angstrom per Bohr (CODATA 2018, same as dft-units)
HA_EV = 27.211386245988  # eV per Hartree (CODATA 2018)
HA_BOHR_TO_EV_A = HA_EV / BOHR_A

# SCF threshold (Ha) written into every arm's [scf] table (run_bench --scf-threshold-ha);
# None leaves the case text unchanged (the engine then resolves its own default,
# which for `baandi relax` follows the force tolerance, ADR 0057).
SCF_THRESHOLD_HA = None

# Primary target 0.02 eV/A and intermediate thresholds, as max atomic force NORM.
THRESHOLDS_EV_A = (0.1, 0.05, 0.02)
TARGET_EV_A = 0.02


def ev_a_to_ha_bohr(x):
    return x / HA_BOHR_TO_EV_A


def load(path):
    with open(path) as f:
        case = json.load(f)
    case["_path"] = os.path.abspath(path)
    return case


def start(case, seed):
    """-> list of (species, fractional) for this seed, and a provenance string."""
    if "starts" in case:
        key = str(seed)
        if key not in case["starts"]:
            raise ValueError("case %s has no pre-perturbed start for seed %s" % (case["name"], seed))
        return [tuple(a) for a in case["starts"][key]], "pre-perturbed start '%s'" % key
    p = case["perturbation"]
    assert p["kind"] == "uniform_per_axis_python_random"
    amp = float(p["amplitude_bohr"])
    lat = case["lattice_bohr"]
    lengths = [math.sqrt(sum(v * v for v in row)) for row in lat]
    # Same generator as exp-tabpfn/tools/relax_cases.py: random.Random(seed),
    # uniform(-amp, amp) Bohr along each lattice axis, in atom order.
    r = random.Random(seed)
    fixed = set(case.get("fixed", []))
    out = []
    for i, (s, f) in enumerate(case["atoms"]):
        d = [r.uniform(-amp, amp) / lengths[k] for k in range(3)]
        out.append((s, list(f) if i in fixed else [f[k] + d[k] for k in range(3)]))
    return out, "random.Random(%d) uniform +-%g Bohr per lattice axis" % (seed, amp)


def pseudo_paths(case, pseudo_root):
    return {s: os.path.join(pseudo_root, rel) for s, rel in case["pseudopotentials"].items()}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scf_text(case, threads):
    """[basis]/[xc]/[scf] text with the harness thread count and SCF threshold (identical for every method)."""
    text = re.sub(r"(?m)^threads\s*=.*\n", "", case["toml_sections"])
    if threads:
        text = text.replace("[scf]\n", "[scf]\nthreads=%d\n" % threads, 1)
    if SCF_THRESHOLD_HA is not None:
        text = re.sub(r"(?m)^standard_threshold_ha\s*=.*\n", "", text)
        text = text.replace("[scf]\n", "[scf]\nstandard_threshold_ha=%r\n" % SCF_THRESHOLD_HA, 1)
    return text


def system_toml(case, species_frac, pseudo_root, threads):
    t = "[system]\nname='%s'\nlattice_unit='bohr'\nlattice_vectors=%s\n" % (
        case["name"], json.dumps(case["lattice_bohr"]))
    for s, f in species_frac:
        t += "[[system.atoms]]\nspecies='%s'\nfractional=%s\n" % (s, json.dumps([float(x) for x in f]))
    t += "[pseudopotentials]\n" + "".join(
        "%s='%s'\n" % (s, p) for s, p in sorted(pseudo_paths(case, pseudo_root).items()))
    return t + scf_text(case, threads)


def max_norm(forces_ha_bohr, fixed):
    """Max atomic force NORM over free atoms (Ha/Bohr). Never a Cartesian component."""
    best = 0.0
    for i, f in enumerate(forces_ha_bohr):
        if i in fixed:
            continue
        best = max(best, math.sqrt(sum(x * x for x in f)))
    return best


def max_component(forces_ha_bohr, fixed):
    """Max |Cartesian component| over free atoms: reported only, never a stop rule here."""
    return max((abs(x) for i, f in enumerate(forces_ha_bohr) if i not in fixed for x in f), default=0.0)
