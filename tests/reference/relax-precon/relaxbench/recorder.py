"""Common external evaluator and per-evaluation ledger (stdlib only).

Every ATTEMPTED DFT evaluation is appended, including failed SCFs, rejected
line-search trials and auxiliary probes (for example ASE's mu estimation).
The stop rule is the same for every method: the first evaluation that is an
optimizer candidate (``role != "auxiliary"``) and whose max atomic force NORM
over free atoms is below the target. Intermediate thresholds are recorded as the
attempted-evaluation count at which they were first met.
"""
import time

from . import cases


class Converged(Exception):
    """Raised by the evaluator once the target is met (stops any optimizer)."""


class BudgetExhausted(Exception):
    pass


class Recorder:
    def __init__(self, fixed, max_evals, target_ev_a=cases.TARGET_EV_A,
                 thresholds_ev_a=cases.THRESHOLDS_EV_A, stop=True):
        self.stop = stop
        self.fixed = set(fixed)
        self.max_evals = max_evals
        self.target = cases.ev_a_to_ha_bohr(target_ev_a)
        self.thresholds = {str(t): cases.ev_a_to_ha_bohr(t) for t in thresholds_ev_a}
        self.first_below = {str(t): None for t in thresholds_ev_a}
        self.evals = []
        self.t0 = time.time()
        self.converged_at = None

    def add(self, rec):
        """rec: dict with role, failed, wall_s, and for successes energy_ha,
        forces_ha_bohr, positions_bohr, scf_iterations, h_vectors."""
        rec = dict(rec)
        rec["index"] = len(self.evals) + 1
        rec["elapsed_s"] = time.time() - self.t0
        if not rec.get("failed"):
            f = rec["forces_ha_bohr"]
            rec["max_force_norm_ha_bohr"] = cases.max_norm(f, self.fixed)
            rec["max_force_component_ha_bohr"] = cases.max_component(f, self.fixed)
        self.evals.append(rec)
        if not rec.get("failed") and rec.get("role") != "auxiliary":
            for key, tol in self.thresholds.items():
                if self.first_below[key] is None and rec["max_force_norm_ha_bohr"] < tol:
                    self.first_below[key] = rec["index"]
            if rec["max_force_norm_ha_bohr"] < self.target and self.converged_at is None:
                self.converged_at = rec["index"]
                if self.stop:
                    raise Converged()
        if len(self.evals) >= self.max_evals:
            raise BudgetExhausted()

    def summary(self):
        ok = [e for e in self.evals if not e.get("failed")]
        last = self.evals[self.converged_at - 1] if self.converged_at else (ok[-1] if ok else None)
        return {
            "attempted_evaluations": len(self.evals),
            "failed_evaluations": sum(1 for e in self.evals if e.get("failed")),
            "auxiliary_evaluations": sum(1 for e in self.evals if e.get("role") == "auxiliary"),
            "line_search_trials": sum(1 for e in self.evals if e.get("role") == "trial"),
            "evaluations_to_threshold_ev_a": self.first_below,
            "converged_at_evaluation": self.converged_at,
            "total_scf_iterations": sum(e.get("scf_iterations") or 0 for e in ok),
            "total_h_vectors": sum(e.get("h_vectors") or 0 for e in ok),
            "evaluation_wall_s": sum(e.get("wall_s") or 0.0 for e in self.evals),
            "final_energy_ha": last["energy_ha"] if last else None,
            "final_max_force_norm_ha_bohr": last["max_force_norm_ha_bohr"] if last else None,
            "final_positions_bohr": last["positions_bohr"] if last else None,
            "last_energy_ha": ok[-1]["energy_ha"] if ok else None,
            "last_max_force_norm_ha_bohr": ok[-1]["max_force_norm_ha_bohr"] if ok else None,
            "last_max_force_component_ha_bohr": ok[-1]["max_force_component_ha_bohr"] if ok else None,
        }
