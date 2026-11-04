#!/usr/bin/env python3
"""Offline synthetic SCF sweep. Exit zero means measured, not converged.

Every child has a 115 s wall cap. No engine defaults or convergence tolerances
are changed. Generated UPF files exist only in a temporary directory.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "tests/reference/stagnation"


def write(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", help="one model ID; runs all selected policies and recovery")
    p.add_argument("--policy", help="one policy ID")
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--iterations", type=int, default=200)
    a = p.parse_args()
    plan = json.loads((BASE / "inputs/sweep.json").read_text())
    models = [m for m in plan["models"] if not a.case or m["id"] == a.case]
    policies = [q for q in plan["policies"] if not a.policy or q["id"] == a.policy]
    if not models or not policies or a.iterations < 2:
        p.error("unknown case or policy")
    a.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, PATH="/opt/homebrew/opt/rustup/bin:" + os.environ["PATH"])
    subprocess.run(["cargo", "build", "--offline", "--locked", "--release", "--manifest-path",
                        str(BASE / "harness/Cargo.toml"), "--target-dir", str(ROOT / "target")],
                       cwd=ROOT, env=env, check=True)
    binary = ROOT / "target/release/stagnation-harness"
    spec = importlib.util.spec_from_file_location("synthetic", ROOT / "scripts/gen_synthetic_upf.py")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    sources = [ROOT / "Cargo.toml", ROOT / ".cargo/config.toml", ROOT / "Cargo.lock", ROOT / "rust-toolchain.toml", Path(__file__),
               BASE / "inputs/sweep.json", BASE / "harness/Cargo.lock",
               BASE / "harness/Cargo.toml", BASE / "harness/src/main.rs",
               ROOT / "scripts/gen_synthetic_upf.py"]
    sources += sorted((ROOT / "crates").glob("*/src/**/*.rs"))
    sources += sorted((ROOT / "crates").glob("*/Cargo.toml"))
    data = {"schema_version": 1, "scope": "nonphysical synthetic NC-format fixtures; same-engine reference, no QE",
            "units": {"energy": "Ha", "residual": "Bohr^-3/2", "lattice": "Bohr", "time": "s"},
            "energy_agreement_tolerance_hartree": 1e-10,
            "command": ["python3", "-I", "scripts/find_stagnation.py"] + sys.argv[1:],
            "platform": platform.platform(), "binary_sha256": digest(binary),
            "parent_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "compiler": subprocess.check_output(["rustc", "--version"], env=env, text=True).strip(),
            "sources": {str(s.relative_to(ROOT)): digest(s) for s in sources}, "runs": []}
    with tempfile.TemporaryDirectory(prefix="stagnation-") as tmp:
        upf = Path(tmp) / "synthetic.UPF"
        def execute(model, policy, arm):
            name = model["id"] + "--" + policy["id"] + "--" + arm
            config = {k: v for k, v in model.items() if k != "generator"}
            config.update({k: v for k, v in policy.items() if k != "id"})
            config.update(id=name, max_iterations=1000 if arm == "reference" else a.iterations, intervention=arm == "intervention")
            inp = a.output / (name + ".input.json")
            out = a.output / (name + ".json")
            write(inp, config)
            start = time.monotonic()
            try:
                process = subprocess.run([str(binary), str(inp.resolve()), str(upf), str(out.resolve())],
                                         cwd=ROOT, capture_output=True, text=True, timeout=115)
                result = json.loads(out.read_text()) if out.exists() else {"state": "ProcessError", "converged": False}
                result["process_exit_code"] = process.returncode
                result["stderr"] = process.stderr
            except subprocess.TimeoutExpired:
                result = {"state": "WallTimeout", "converged": False, "iterations": None,
                          "final_residual": None, "findings": [], "error": "115 s cap; no physical diagnosis inferred"}
            result.update(model=model["id"], policy=policy["id"], arm=arm,
                          input_sha256=digest(inp), upf_sha256=digest(upf),
                          generator=model["generator"], elapsed_seconds=time.monotonic()-start)
            write(out, result)
            data["runs"].append(result)
            write(a.output / "results.json", data)
            print(name, result["state"], result.get("iterations"), result.get("final_residual"), flush=True)
            return result
        for model in models:
            upf.write_text(generator.generate(**model["generator"]))
            reference = execute(model, {"id": "reference", "mixer": "pulay", "beta": 0.1, "solver": "davidson"}, "reference")
            for policy in policies:
                before = execute(model, policy, "baseline")
                if any(f["rule"] in ("Stagnation", "Oscillation") for f in before.get("findings", [])):
                    after = execute(model, policy, "intervention")
                    e = after.get("certified_energy_hartree")
                    ref = reference.get("certified_energy_hartree")
                    after["energy_error_hartree"] = abs(e-ref) if e is not None and ref is not None else None
                    after["recovery_validated"] = bool(after.get("intervention_applied") and e is not None and ref is not None and abs(e-ref) <= 1e-10)
                    write(a.output / (after["model"] + "--" + after["policy"] + "--intervention.json"), after)
                    write(a.output / "results.json", data)
    rows = ["# Measured sweep", "", "All energies are Ha; residuals are Bohr^-3/2. Full windows/audit in JSON.", "",
            "| Model | Policy | Arm | State | Iterations | Final residual | Certified energy | P11D triggered rules | Recovery ΔE | Sampled gap | Validated intervention |", "|---|---|---|---|---:|---:|---:|---|---:|---:|---|"]
    for r in data["runs"]:
        rules = ", ".join(sorted({f["rule"] for f in r.get("findings", [])})) or "none"
        rows.append("| " + " | ".join(str(x) for x in [r["model"], r["policy"], r["arm"], r["state"], r.get("iterations"), r.get("final_residual"), r.get("certified_energy_hartree"), rules, r.get("energy_error_hartree"), r.get("sampled_gap_hartree"), r.get("recovery_validated")]) + " |")
    (a.output / "table.md").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
