"""ASE calculator around ``baandi run`` ([output] forces=true): one COLD SCF per call.

``baandi run`` has no start-density input, so every evaluation starts from the
atomic superposition. Warm starts are therefore impossible on this path; runs
through this calculator are labelled ``start = "cold"`` and must only be compared
in wall time with other cold runs (spec section 6). Evaluation counts and SCF
iterations are recorded for every attempt, including failures.

ASE is imported here only: the engine never depends on Python or ASE.
"""
import json
import os
import shutil
import subprocess
import time

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from . import cases


class BaandiCLI(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, case, exe, pseudo_root, threads, workdir, recorder, main_atoms=None, keep=False, **kw):
        super().__init__(**kw)
        self.case, self.exe, self.pseudo_root, self.threads = case, exe, pseudo_root, threads
        self.workdir, self.recorder, self.keep = workdir, recorder, keep
        self.lat = np.array(case["lattice_bohr"], float)  # rows are lattice vectors
        self.main_atoms = main_atoms
        self.role = "candidate"
        self.memo, self.memo_hits = {}, 0
        os.makedirs(workdir, exist_ok=True)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        # Exact-geometry memo: ASE's single-state cache forgets the main structure
        # after an auxiliary probe on a copy (Exp mu estimation) and would rerun an
        # identical SCF. A re-query of an already evaluated geometry is not a new
        # DFT evaluation for any method; hits are counted separately.
        key = atoms.get_positions().tobytes() + atoms.cell.array.tobytes()
        if key in self.memo:
            self.results = dict(self.memo[key])
            self.memo_hits += 1
            return
        pos_bohr = atoms.get_positions() / cases.BOHR_A
        frac = pos_bohr @ np.linalg.inv(self.lat)
        species = atoms.get_chemical_symbols()
        text = cases.system_toml(self.case, list(zip(species, frac.tolist())), self.pseudo_root, self.threads)
        text += "[output]\nforces=true\n"
        n = len(self.recorder.evals)
        d = os.path.join(self.workdir, "eval%04d" % n)
        if os.path.exists(d):
            shutil.rmtree(d)
        os.makedirs(d)
        with open(os.path.join(d, "in.toml"), "w") as f:
            f.write(text)
        role = self.role
        if self.main_atoms is not None and atoms is not self.main_atoms:
            role = "auxiliary"  # e.g. ASE Exp mu estimation on a perturbed copy
        t0 = time.time()
        proc = subprocess.run(list(self.exe) + ["run", "in.toml", "--output-dir", "out"], cwd=d,
                              capture_output=True, text=True)
        wall = time.time() - t0
        res = None
        rp = os.path.join(d, "out", "results.json")
        if os.path.exists(rp):
            try:
                with open(rp) as f:
                    res = json.load(f)
            except ValueError:
                res = None
        ok = (res is not None and res.get("convergence", {}).get("converged")
              and res.get("forces", {}).get("status") == "available")
        rec = {"role": role, "wall_s": wall, "positions_bohr": pos_bohr.tolist(), "exit_code": proc.returncode}
        if not ok:
            rec["failed"] = True
            rec["reason"] = (proc.stderr or "")[-400:]
        else:
            rec.update({
                "energy_ha": res["total_energy"]["ha"],
                "forces_ha_bohr": res["forces"]["values"],
                "scf_iterations": res["convergence"]["iterations"],
                "h_vectors": res.get("work", {}).get("h_vectors"),
                "verified": res.get("verified"),
            })
        if not self.keep:
            shutil.rmtree(os.path.join(d, "out"), ignore_errors=True)
        if ok:
            f = np.array(rec["forces_ha_bohr"], float)
            self.results = {"energy": rec["energy_ha"] * cases.HA_EV,
                            "forces": f * cases.HA_BOHR_TO_EV_A}
            self.memo[key] = dict(self.results)
            if len(self.memo) > 8:
                self.memo.pop(next(iter(self.memo)))
        # The recorder may raise Converged / BudgetExhausted: results are already set.
        self.recorder.add(rec)
        if not ok:
            raise RuntimeError("BAANDI SCF/forces failed at evaluation %d (exit %d)" % (n + 1, proc.returncode))
