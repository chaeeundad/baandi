#!/usr/bin/env python3
"""Mac correctness only: no external clock, sampling, or performance summaries."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items()
                if k != "timings" and not k.endswith("_s")}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    if isinstance(value, float):
        return struct.pack("!d", value).hex()
    return value


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--upf", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    template = Path("examples/si_lda.toml").read_text().replace("/absolute/path/to/Si.UPF", str(args.upf.resolve()))
    cases = {
        "gamma": template,
        "skew": template.replace("[0.0, 2.715, 2.715]", "[0.08, 2.715, 2.715]"),
        "anisotropic": template.replace("[2.715, 2.715, 0.0]", "[2.79645, 2.79645, 0.0]"),
        "k2": template.replace("kpoints = [[0.0, 0.0, 0.0]]", 'kpoints = [[0.0, 0.0, 0.0, 0.5], [0.25, 0.0, 0.0, 0.5]]'),
    }
    evidence = []
    for name, text in cases.items():
        oracle = None
        for mode, binary in (("baseline", args.baseline), ("default", args.candidate), ("off", args.candidate), ("on", args.candidate), ("input", args.candidate)):
            input_path = args.output / f"{name}-{mode}.toml"
            text_arm = text + "\n[scf]\nthreads = 1\n"
            if mode == "input":
                text_arm += "\n[profile]\nphase_profile = true\n"
            input_path.write_text(text_arm)
            env = {k: v for k, v in os.environ.items() if not k.startswith("BAANDI_")}
            if mode not in ("input", "default"):
                env["BAANDI_PHASE_PROFILE"] = "on" if mode == "on" else "off"
            output = args.output / f"{name}-{mode}"
            command = [str(binary.resolve()), "run", str(input_path.resolve()), "--output-dir", str(output.resolve())]
            with (args.output / f"{name}-{mode}.log").open("w") as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
            if result.returncode:
                raise RuntimeError(f"{name}/{mode}: exit {result.returncode}; see log")
            doc = json.loads((output / "results.json").read_text())
            assert doc["convergence"]["converged"] and doc["verified"]
            fields = {k: scientific(doc[k]) for k in ("total_energy", "energies", "eigenvalues", "occupations", "verification")}
            fields["history"] = scientific(doc["convergence"]["history"])
            fields["density_bin"] = digest(output / "density.bin")
            fields["density_json"] = digest(output / "density.json")
            if oracle is None:
                oracle = fields
            assert fields == oracle, f"{name}/{mode}: scientific bits changed: {[k for k in fields if fields[k] != oracle[k]]}"
            assert ("launchio" in doc) == (mode in ("on", "input"))
            if mode in ("on", "input"):
                subprocess.run(["python3", "-I", "scripts/lane_launchio_parse.py", str(output), "--output", str(output / "parser-smoke.json")], check=True)
            evidence.append({"case": name, "mode": mode, "exit_code": result.returncode,
                             "energy_ha": doc["total_energy"]["ha"], "delta_energy_ha": 0,
                             "density_sha256": fields["density_bin"],
                             "iterations": doc["convergence"]["iterations"], "verified": doc["verified"],
                             "bit_identical": True})
    report = {"purpose": "correctness only; no Mac performance evidence", "baseline_sha256": digest(args.baseline),
              "candidate_sha256": digest(args.candidate), "upf_sha256": digest(args.upf), "runs": evidence}
    (args.output / "equivalence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(evidence)} executions; science and density bit-identical; 8 parser smoke checks")


if __name__ == "__main__":
    main()
