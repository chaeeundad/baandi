"""BAANDI's own optimizers (``baandi relax``, method bfgs, fire or precon_lbfgs) under the common evaluator.

The native relax prints one line per accepted ionic step after it has published
results.json (atomic rename). The harness reads each line as it arrives, replays
the new frames and failure records into the common Recorder, and terminates the
process as soon as the external target (max atomic force norm) is met. The native
force/energy stop rule is disabled by a force tolerance below any reachable value
(1e-9 Ha/Bohr), while the energy tolerance keeps the case's value because it also
drives BFGS trust-radius halving: changing it would change the optimizer.
Per-evaluation wall comes from line arrival times (each includes output writing).

precon_lbfgs (ADR 0058) also evaluates line-search trials that it does not
accept; they are listed in ``relaxation.rejected_trials`` and published at once
(stdout line ``trial ...``). Every frame, failure and trial of that method carries
a global ``evaluation`` index, which orders the ledger; trials enter the recorder
with role ``trial`` (they are candidates for the common stop rule, as the A3 and
ASE line-search trials are).
Stdlib only.
"""
import json
import os
import signal
import subprocess
import time

from . import cases
from .recorder import BudgetExhausted, Converged


NATIVE_METHODS = {"bfgs": "bfgs", "fire": "fire", "precon": "precon_lbfgs"}


def relax_toml(case, species_frac, pseudo_root, threads, method, warm, max_steps, own_stop=False):
    r = case["relax"]
    ftol = r["force_tol_ha_bohr"] if own_stop else 1e-9
    t = cases.system_toml(case, species_frac, pseudo_root, threads)
    t += ("[relax]\nmethod='%s'\nforce_tol_ha_bohr=%r\nenergy_tol_ha=%r\nmax_steps=%d\n"
          "max_displacement_bohr=%r\ntrust_radius_bohr=%r\nwarm_start=%s\n" % (
              NATIVE_METHODS[method], ftol, r["energy_tol_ha"], max_steps, r["max_displacement_bohr"],
              r["trust_radius_bohr"], "true" if warm else "false"))
    if case.get("fixed"):
        t += "fixed=[%s]\n" % ",".join("{atom=%d}" % i for i in case["fixed"])
    return t


def ledger(report):
    """Attempt sequence (failures of step k precede frame k) from a relax results.json.

    precon_lbfgs: frames (optimizer.evaluation), failures and rejected trials
    (evaluation) are ordered by the global evaluation index; initial-structure
    failures (step 0, no index) come first."""
    rel = report["relaxation"]
    frames, failures = rel.get("frames", []), rel.get("failures", [])
    if "optimizer" in rel:
        keyed = [(fr["optimizer"]["evaluation"], "frame", fr) for fr in frames]
        keyed += [(x["evaluation"], "failure", x) for x in failures if "evaluation" in x]
        keyed += [(x["evaluation"], "trial", x) for x in rel.get("rejected_trials", [])]
        head = [("failure", x) for x in failures if "evaluation" not in x]
        return head + [(kind, x) for _, kind, x in sorted(keyed, key=lambda k: k[0])]
    seq = []
    for k, fr in enumerate(frames):
        seq += [("failure", x) for x in failures if x.get("step") == k]
        seq.append(("frame", fr))
    seq += [("failure", x) for x in failures if x.get("step", -1) >= len(frames)]
    return seq


def run(case, species_frac, exe, pseudo_root, threads, method, warm, workdir, recorder, own_stop=False):
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "in.toml"), "w") as f:
        f.write(relax_toml(case, species_frac, pseudo_root, threads, method, warm, recorder.max_evals, own_stop))
    out = os.path.join(workdir, "out")
    t0 = time.time()
    err = open(os.path.join(workdir, "stderr.txt"), "w")
    proc = subprocess.Popen(list(exe) + ["relax", "in.toml", "--output-dir", "out"], cwd=workdir,
                            stdout=subprocess.PIPE, stderr=err, text=True)
    state = {"seen": 0, "last_t": t0}
    accepted = []

    def replay():
        path = os.path.join(out, "results.json")
        if not os.path.exists(path):
            return
        now = time.time()
        with open(path) as f:
            report = json.load(f)
        seq = ledger(report)
        new = seq[state["seen"]:]
        state["seen"] = len(seq)
        share = (now - state["last_t"]) / max(len(new), 1)
        state["last_t"] = now
        for kind, x in new:
            if kind == "failure":
                recorder.add({"role": "candidate", "failed": True, "wall_s": share,
                              "positions_bohr": x.get("positions_bohr"), "reason": x.get("reason")})
                continue
            if kind == "trial":
                p = x["positions_bohr"]
                recorder.add({"role": "trial", "wall_s": share,
                              "positions_bohr": [p[i:i + 3] for i in range(0, len(p), 3)],
                              "energy_ha": x["energy_ha"], "forces_ha_bohr": x["forces_ha_bohr"],
                              "scf_iterations": x["scf_iterations"], "h_vectors": x.get("scf_h_vectors"),
                              "provenance": x.get("provenance")})
                continue
            pos = [x["positions"][i:i + 3] for i in range(0, len(x["positions"]), 3)]
            accepted.append(pos)
            recorder.add({"role": "candidate", "wall_s": share, "positions_bohr": pos,
                          "energy_ha": x["energy_ha"], "forces_ha_bohr": x["forces_ha_bohr"],
                          "scf_iterations": x["scf_iterations"], "h_vectors": x.get("scf_h_vectors"),
                          "provenance": x.get("provenance")})

    status = None
    try:
        for line in proc.stdout:
            if line.startswith("step ") or line.startswith("trial "):
                replay()
        proc.wait()
        replay()  # failures without a following accepted step
        status = "native_exit_%d" % proc.returncode
    except Converged:
        status = "converged"
    except BudgetExhausted:
        status = "budget_exhausted"
    finally:
        if proc.poll() is None:
            proc.send_signal(signal.SIGTERM)
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
        err.close()
    return status, time.time() - t0, accepted
