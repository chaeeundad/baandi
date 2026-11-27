#!/usr/bin/env python3
"""Phase 1 driver (PRELAX-PRECON): one (case, seed, method, start) relaxation.

Run with an EXTERNAL Python that has ASE (>= 3.23) and SciPy; the BAANDI engine
never depends on Python/ASE. Example (Mac smoke):

  PY=/Volumes/990-2TB/Develop/rust-dft-wt/venv-convert/bin/python
  $PY -I tests/reference/relax-precon/run_bench.py --case tests/reference/relax-precon/cases/si8-lda-g24.json \
      --seed 1 --method ase-precon-lbfgs --start cold --exe target/release/baandi \
      --pseudo-root /Volumes/990-2TB/Develop/rust-dft-wt/pseudo --threads 1 \
      --workdir /scratch/relax-precon/run --out results/x.json

Methods: baandi-bfgs, baandi-fire (``baandi relax``; --start cold|warm),
ase-bfgs, ase-lbfgs, ase-fire, ase-precon-lbfgs, a3-lbfgs (``baandi run`` per
evaluation through an ASE calculator; --start cold only, see relaxbench/clicalc.py).
All methods get the same structure, PPs, XC, cutoff, k points, SCF settings,
threads, fixed atoms, per-atom displacement cap (the case's max_displacement_bohr)
and the same external stop rule: the first optimizer candidate whose max atomic
force NORM over free atoms is below 0.02 eV/A (3.889e-4 Ha/Bohr). Every attempted
DFT evaluation is counted (failures, rejected line-search trials, auxiliary probes).
"""
import argparse
import datetime
import json
import os
import platform
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from relaxbench import cases  # noqa: E402
from relaxbench.recorder import BudgetExhausted, Converged, Recorder  # noqa: E402

METHODS = ("baandi-bfgs", "baandi-fire", "ase-bfgs", "ase-lbfgs", "ase-fire", "ase-precon-lbfgs", "a3-lbfgs")


def mark_accepted(evals, accepted_bohr, tol=1e-9):
    for e in evals:
        if e.get("failed"):
            e["accepted"] = False
            continue
        if e.get("role") == "auxiliary":
            e["accepted"] = None
            continue
        p = e["positions_bohr"]
        e["accepted"] = any(
            max(abs(a - b) for ra, rb in zip(p, q) for a, b in zip(ra, rb)) < tol for q in accepted_bohr)


def run_ase(args, case, start, rec, info):
    import numpy as np
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS, FIRE, LBFGS
    from relaxbench.clicalc import BaandiCLI

    lat_a = np.array(case["lattice_bohr"]) * cases.BOHR_A
    atoms = Atoms([s for s, _ in start], scaled_positions=[f for _, f in start], cell=lat_a, pbc=True)
    fixed = case.get("fixed", [])
    if fixed:
        atoms.set_constraint(FixAtoms(indices=fixed))
    calc = BaandiCLI(case, args.cmd, args.pseudo_root, args.threads, os.path.join(args.workdir, "evals"),
                     rec, main_atoms=atoms)
    atoms.calc = calc
    cap_a = case["relax"]["max_displacement_bohr"] * cases.BOHR_A
    accepted = []
    status = "stopped"
    try:
        if args.method == "a3-lbfgs":
            from relaxbench import a3
            free = np.ones(len(atoms), bool)
            free[fixed] = False
            info["method_parameters"] = {"memory": 10, "c1": 1e-4, "backtrack": 0.5, "max_trials": 25,
                                         "cap_mode": "per_atom", "cap_angstrom": cap_a,
                                         "metric": "fixed bond-spring model Hessian of the start geometry"}
            st, minfo = a3.run(atoms, free, cap_a, accepted)
            info["method_parameters"]["metric_info"] = minfo
            status = st
        else:
            if args.method == "ase-precon-lbfgs":
                from ase.optimize.precon import Exp, PreconLBFGS
                opt = PreconLBFGS(atoms, precon=Exp(A=3.0), use_armijo=True, maxstep=cap_a,
                                  logfile=os.path.join(args.workdir, "opt.log"))
                info["method_parameters"] = {"precon": "Exp(A=3.0) defaults (r_cut=2 r_NN, c_stab=0.1, mu estimated once)",
                                             "use_armijo": True, "c1": 0.23, "memory": 100, "maxstep_angstrom": cap_a}
            else:
                cls = {"ase-bfgs": BFGS, "ase-lbfgs": LBFGS, "ase-fire": FIRE}[args.method]
                opt = cls(atoms, maxstep=cap_a, logfile=os.path.join(args.workdir, "opt.log"))
                info["method_parameters"] = {"class": cls.__name__, "maxstep_angstrom": cap_a,
                                             "other": "ASE defaults"}
            for _ in opt.irun(fmax=1e-12, steps=10 ** 6):
                accepted.append(atoms.get_positions().copy())
            status = "optimizer_returned"
    except Converged:
        status = "converged"
    except BudgetExhausted:
        status = "budget_exhausted"
    except RuntimeError as e:
        status = "failed: %s" % str(e)[:200]
    info["memo_hits"] = calc.memo_hits
    return status, [(p / cases.BOHR_A).tolist() for p in accepted]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--method", required=True, choices=METHODS)
    ap.add_argument("--start", required=True, choices=("cold", "warm"))
    ap.add_argument("--exe", required=True)
    ap.add_argument("--pseudo-root", required=True)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--launcher", default="",
                    help="command prefix for every baandi process, identical for all methods "
                         "(e.g. 'taskset -c 3' or 'mpirun -np 8'); empty = direct")
    ap.add_argument("--max-evals", type=int, default=200)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--native-own-stop", action="store_true",
                    help="baandi-* only: keep the case's native stop rule (max component + energy, twice) and "
                         "do not stop at the external target (baseline audit); the target is still recorded")
    args = ap.parse_args()
    args.exe = os.path.abspath(args.exe)
    args.cmd = args.launcher.split() + [args.exe]
    if args.start == "warm" and not args.method.startswith("baandi-"):
        ap.error("ASE/A3 methods evaluate through `baandi run`, which has no warm start: use --start cold")

    case = cases.load(args.case)
    start, start_note = cases.start(case, args.seed)
    if args.native_own_stop and not args.method.startswith("baandi-"):
        ap.error("--native-own-stop applies to baandi-* methods only")
    rec = Recorder(case.get("fixed", []), args.max_evals, stop=not args.native_own_stop)
    info = {}
    os.makedirs(args.workdir, exist_ok=True)
    t0 = time.time()
    if args.method.startswith("baandi-"):
        from relaxbench import native
        status, wall, accepted = native.run(case, start, args.cmd, args.pseudo_root, args.threads,
                                            args.method.split("-")[1], args.start == "warm",
                                            args.workdir, rec, args.native_own_stop)
        info["method_parameters"] = {"relax": case["relax"], "native_own_stop": args.native_own_stop,
                                     "note": ("native stop rule as configured (baseline audit)" if args.native_own_stop
                                              else "native force_tol 1e-9: external evaluator stops the run")}
    else:
        status, accepted = run_ase(args, case, start, rec, info)
        wall = time.time() - t0
    mark_accepted(rec.evals, accepted)
    if rec.converged_at:
        # The evaluator stopped the optimizer at this candidate before it could accept it.
        rec.evals[rec.converged_at - 1]["accepted"] = True
    summary = rec.summary()
    summary["accepted_steps"] = sum(1 for e in rec.evals if e.get("accepted"))
    summary["rejected_evaluations"] = sum(1 for e in rec.evals if e.get("accepted") is False)
    summary["end_to_end_wall_s"] = wall
    summary["status"] = status
    try:
        import ase
        ase_version = ase.__version__
    except ImportError:
        ase_version = None
    result = {
        "schema": "relax-precon-run-1",
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "case": case["name"], "case_file": os.path.relpath(case["_path"]),
        "seed": args.seed, "start_structure": start_note,
        "method": args.method, "start": args.start, "native_own_stop": args.native_own_stop,
        "stop_rule": {"quantity": "max atomic force norm over free atoms", "target_ev_a": cases.TARGET_EV_A,
                      "target_ha_bohr": cases.ev_a_to_ha_bohr(cases.TARGET_EV_A),
                      "thresholds_ev_a": list(cases.THRESHOLDS_EV_A), "max_evals": args.max_evals},
        "threads": args.threads, "launcher": args.launcher,
        "executable": args.exe, "executable_sha256": cases.sha256_file(args.exe),
        "pseudopotentials_sha256": {s: cases.sha256_file(p) for s, p in cases.pseudo_paths(case, args.pseudo_root).items()},
        "host": {"platform": platform.platform(), "machine": platform.machine(), "node": platform.node(),
                 "python": platform.python_version(), "ase": ase_version},
        "info": info,
        "summary": summary,
        "evaluations": [{k: v for k, v in e.items() if k not in ("positions_bohr", "forces_ha_bohr")}
                        for e in rec.evals],
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=1)
    s = summary
    print("%s seed %d %s/%s: %s evals=%d (failed %d, rejected %d, aux %d) thr=%s scf=%d wall=%.1fs E=%s Fn=%s" % (
        case["name"], args.seed, args.method, args.start, status, s["attempted_evaluations"],
        s["failed_evaluations"], s["rejected_evaluations"], s["auxiliary_evaluations"],
        s["evaluations_to_threshold_ev_a"], s["total_scf_iterations"], wall,
        s["final_energy_ha"], s["final_max_force_norm_ha_bohr"]))


if __name__ == "__main__":
    main()
