"""Mac diagnostic identity/byte observations; no speed-adoption claim.

Run from the repository root, with all scratch in --output. The external wall
starts before process launch and ends after process exit and stdout/stderr close.
No output is discarded. Baseline mode records the pre-edit reference executable.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def input_text(upf, rx, ry):
    a = 10.261212856717933
    base = [(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0),
            (.25, .25, .25), (.25, .75, .75), (.75, .25, .75), (.75, .75, .25)]
    text = f'''[system]
name = "diagnostic identity"
lattice_unit = "bohr"
lattice_vectors = [[{a*rx},0,0],[0,{a*ry},0],[0,0,{a}]]
'''
    i = 0
    for x in range(rx):
        for y in range(ry):
            for p in base:
                coords = [(p[0]+x)/rx + 1e-4*math.sin(i*7+1),
                          (p[1]+y)/ry + 1e-4*math.sin(i*11+3),
                          p[2] + 1e-4*math.sin(i*13+5)]
                text += f'[[system.atoms]]\nspecies="Si"\nfractional={coords}\n'
                i += 1
    return text + f'''[pseudopotentials]
Si = {json.dumps(str(upf))}
[basis]
ecutwfc_ry = 30.0
kpoints = [[0.0,0.0,0.0]]
[xc]
functional = "lda-pz"
[scf]
driver = "reference"
nbands = {2*i+4}
max_iterations = 200
mixing = "broyden"
mixing_beta = 0.7
mixing_history = 8
threads = 1
'''


def scientific(value):
    """Strip only clocks and explicitly optional diagnostic payloads."""
    if isinstance(value, list):
        return [scientific(v) for v in value]
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items()
                if k not in {"timings", "wall_s", "ritz_changes_ha", "inner_clocks"}
                and not k.endswith("_s")}
    return value


def science(doc):
    keys = ["total_energy", "energies", "eigenvalues", "occupations", "verified",
            "verification", "basis", "electrons", "bands", "single_path"]
    result = {k: scientific(doc[k]) for k in keys}
    result["history"] = scientific(doc["convergence"]["history"])
    result["reference_iterations"] = scientific(doc["reference_scf"]["iterations"])
    result["verification_calls"] = scientific(doc["reference_scf"]["run"]["verification_calls"])
    result["first_stop"] = doc["reference_scf"]["run"]["first_estimator_stop_iteration"]
    result["exit_reason"] = doc["reference_scf"]["exit_reason"]
    return result


def calls(doc):
    ref = doc["reference_scf"]
    return [c for r in ref["iterations"] for c in r["calls"]] + ref["run"]["verification_calls"]


def ritz_arrays(value):
    if isinstance(value, dict):
        for k, v in value.items():
            if k == "ritz_changes_ha":
                yield v
            else:
                yield from ritz_arrays(v)
    elif isinstance(value, list):
        for v in value:
            yield from ritz_arrays(v)


def measure(binary, inp, dest, mode, baseline=False):
    dest.mkdir(parents=True, exist_ok=False)
    env = {k: v for k, v in os.environ.items() if not k.startswith("BAANDI_")}
    # Same input bytes: this existing run-owned override is the only A/B change.
    env["BAANDI_DAVIDSON_STAGES"] = mode
    env["TMPDIR"] = str(dest)
    cmd = [str(binary), "run", str(inp), "--output-dir", str(dest)]
    load = os.getloadavg()
    start = time.perf_counter()
    with (dest / "stdout").open("wb") as stdout, (dest / "stderr").open("wb") as stderr:
        run = subprocess.run(cmd, env=env, stdout=stdout, stderr=stderr)
    wall = time.perf_counter() - start
    assert run.returncode == 0, (cmd, run.returncode, dest / "stderr")
    doc = json.loads((dest / "results.json").read_text())
    assert doc["run"]["completed"] and doc["verified"] and doc["convergence"]["converged"]
    assert doc["standard_threshold_ha"] == 5e-9
    arrays = list(ritz_arrays(doc))
    if not baseline:
        assert doc["reference_scf"]["diagnostics"]["level"] == ("full" if mode == "on" else "lightweight")
        for c in calls(doc):
            trace = c["ritz_changes_ha"]
            if mode == "on":
                assert len(trace) == c["davidson_iterations"]
                assert all(len(row) == doc["bands"] for row in trace)
                assert c["inner_clocks"]["total_s"] > 0
                assert any(v > 0 for v in c["inner_clocks"]["stage_s"].values())
            else:
                assert trace == [] and c["inner_clocks"] is None
        assert all(bool(a) == (mode == "on") for a in arrays)
    files = {p.name: p.stat().st_size for p in dest.iterdir() if p.is_file()}
    summary = dict(command=cmd, exit_code=run.returncode, mode=mode,
                   external_wall_s=wall, load_before=load, input_sha256=sha(inp),
                   density_sha256=sha(dest / "density.bin"), energy_ha=doc["total_energy"]["ha"],
                   steps=doc["convergence"]["iterations"],
                   verification=scientific(doc["verification"]),
                   ritz_f64_payload_bytes=8*sum(len(row) for c in calls(doc) for row in c["ritz_changes_ha"]),
                   ritz_json_value_bytes_all_copies=sum(len(json.dumps(a, separators=(",", ":"))) for a in arrays),
                   output_bytes=sum(files.values()), files=files)
    return summary, science(doc)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--upf", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--baseline", action="store_true")
    p.add_argument("--compare-baseline", type=Path)
    args = p.parse_args()
    binary, upf, output = (x.resolve() for x in (args.binary, args.upf, args.output))
    output.mkdir(parents=True, exist_ok=True)
    evidence = dict(platform=platform.platform(), binary_sha256=sha(binary), upf_sha256=sha(upf),
                    scope="Mac serial correctness and byte observations, not dev08 speed evidence", cases={})
    for name, rx, ry in [("si8", 1, 1), ("si16", 2, 1), ("si32", 2, 2)]:
        inp = output / f"{name}.toml"
        inp.write_text(input_text(upf, rx, ry))
        samples = []
        previous = None
        # Crossed pairs, all output and every verification attempt included.
        modes = ["off"] if args.baseline else ["off", "on", "on", "off"]
        for i, mode in enumerate(modes):
            dest = output / f"{name}-{i}-{mode}"
            sample, current = measure(binary, inp, dest, mode, args.baseline)
            if previous is not None:
                assert current == previous, f"scientific mismatch {name} sample {i}"
                assert sample["density_sha256"] == samples[0]["density_sha256"]
            if args.compare_baseline:
                old = args.compare_baseline / f"{name}-0-off"
                old_doc = json.loads((old / "results.json").read_text())
                assert current == science(old_doc)
                assert sample["density_sha256"] == sha(old / "density.bin")
                if mode == "on":
                    new_doc = json.loads((dest / "results.json").read_text())
                    assert list(ritz_arrays(new_doc)) == list(ritz_arrays(old_doc))
            previous = current
            samples.append(sample)
            print(name, mode, f"{sample['external_wall_s']:.6f} s", sample["output_bytes"], "bytes", flush=True)
        evidence["cases"][name] = dict(samples=samples, scientific_identity=True,
                                        baseline_identity=bool(args.compare_baseline),
                                        baseline_full_trace_identity=bool(args.compare_baseline and not args.baseline))
        (output / "evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")


if __name__ == "__main__":
    main()
