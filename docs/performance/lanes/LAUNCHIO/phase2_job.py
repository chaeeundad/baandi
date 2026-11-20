#!/usr/bin/env python3
"""Master-only dev08 A/A floor and crossed A/B; one case/model per fresh job."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def paired_stats(pairs):
    savings = [a - b for a, b in pairs]
    mean = statistics.mean(savings)
    sd = statistics.stdev(savings) if len(savings) > 1 else None
    se = sd / math.sqrt(len(savings)) if sd is not None else None
    return {"pairs": len(pairs), "mean_saving_s": mean, "sd_s": sd, "se_s": se,
            "abs_mean_over_se": abs(mean) / se if se else None,
            "median_a_s": statistics.median(a for a, _ in pairs),
            "median_b_s": statistics.median(b for _, b in pairs), "savings_s": savings}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", choices=["si16", "si32", "si64"], required=True)
    p.add_argument("--ranks", type=int, choices=[1, 8], required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if platform.system() != "Linux" or platform.node().split(".")[0] != "dev08":
        p.error("remote timing is restricted to the master on dev08")
    if args.ranks == 1 and args.case != "si64":
        p.error("this Phase 2 request covers only si64 one core")
    import tomllib
    source, input_path = args.binary.resolve(), args.input.resolve()
    config = tomllib.loads(input_path.read_text())
    if len(config["system"]["atoms"]) != int(args.case[2:]) or config.get("scf", {}).get("threads", 1) != 1:
        p.error("case size mismatch or threads != 1")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = root / "dft-cli"
    shutil.copy2(source, binary)
    expected_hash = sha(binary) # independent master oracle, identical for A and B
    (root / "input.original.toml").write_bytes(input_path.read_bytes())
    env = os.environ.copy()
    env.update(MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    spec = importlib.util.spec_from_file_location("identity", Path(__file__).with_name("verify_phase2.py"))
    identity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(identity)
    launcher = ["taskset", "-c", "0"] if args.ranks == 1 else [
        "mpirun", "-np", "8", "--cpu-set", "0-7", "--map-by", "core", "--bind-to", "core", "--report-bindings"]
    manifest = {"case": args.case, "ranks": args.ranks, "source_binary": str(source),
                "binary": str(binary), "binary_sha256": expected_hash, "binary_bytes": binary.stat().st_size,
                "input": str(input_path), "input_sha256": sha(input_path),
                "upf_sha256": {name: sha(input_path.parent / path) for name, path in config["pseudopotentials"].items()},
                "environment": {k: v for k, v in env.items() if k.startswith(("BAANDI_", "MKL_", "OMP_"))},
                "commands": [], "runs": []}
    reference = {}

    def save():
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")

    def run(label, arm, profile=False, cold=False):
        deadline = time.monotonic() + 600
        while os.getloadavg()[0] >= 0.30:
            if time.monotonic() > deadline:
                raise RuntimeError("load gate did not clear; reschedule, preserve this failed job")
            time.sleep(10)
        out = root / label
        env.update(BAANDI_BINARY_HASH_CACHE="on" if arm == "B" else "off",
                   BAANDI_PHASE_PROFILE="on" if profile else "off", BAANDI_STAGE_TIMERS="on" if profile else "off")
        command = launcher + [str(binary), "run", str(input_path), "--output-dir", str(out)]
        manifest["commands"].append({"label": label, "arm": arm, "profile": profile, "argv": command})
        save()
        before = list(os.getloadavg())
        with (root / (label + ".stdout")).open("w") as stdout, (root / (label + ".stderr")).open("w") as stderr:
            start = time.monotonic()
            proc = subprocess.run(command, env=env, stdout=stdout, stderr=stderr)
            wall = time.monotonic() - start
        entry = {"label": label, "arm": arm, "profile": profile, "wall_s": wall,
                 "exit_code": proc.returncode, "load_before": before, "load_after": list(os.getloadavg())}
        manifest["runs"].append(entry)
        save()
        if proc.returncode:
            raise RuntimeError(f"{label}: process exit {proc.returncode}")
        paths = [out / "results.json", *sorted(out.glob("rank*/results.json"))]
        if len(paths) != args.ranks:
            raise RuntimeError(f"{label}: missing rank outputs")
        seen = set()
        for path in paths:
            doc = json.loads(path.read_text())
            rank = doc["run"]["rank"]
            seen.add(rank)
            assert doc["run"]["completed"] and doc["convergence"]["converged"] and doc["verified"]
            assert doc["run"]["binary"]["sha256"] == expected_hash
            if args.ranks > 1:
                assert doc["parallel"]["model"] != "mpi_replicated"
            c = doc["run"]["binary"].get("hash_cache")
            if arm == "B":
                assert c["status"] in (("miss", "hit") if cold else ("hit",))
                assert c["hashed_bytes"] == (binary.stat().st_size if c["status"] == "miss" else 0)
            else:
                assert c is None
            density_hash = sha(path.with_name("density.bin"))
            metadata_hash = sha(path.with_name("density.json"))
            assert density_hash == json.loads(path.with_name("density.json").read_text())["sha256"]
            # Full normalized JSON for uninstrumented acceptance. The diagnostic
            # runs also verify all science against A, while their nested profiles
            # are read separately rather than compared as deterministic data.
            if profile:
                compare = {k: doc[k] for k in ("total_energy", "energies", "eigenvalues", "occupations", "verification")}
                compare["history"] = doc["convergence"]["history"]
            else:
                compare = doc
            normalized = identity.normalized(compare)
            key = (rank, profile)
            state = (normalized, density_hash, metadata_hash)
            if key in reference:
                assert state == reference[key], f"{label}/rank{rank}: output identity failure"
            else:
                reference[key] = state
        assert seen == set(range(args.ranks))
        return wall

    # New private binary path => sidecar absent. Cold cost is retained separately.
    run("warmup-A", "A")
    run("cold-B", "B", cold=True)
    run("warmup-B", "B")
    aa = [(run(f"aa-{i}-A0", "A"), run(f"aa-{i}-A1", "A")) for i in range(6)]
    floor = paired_stats(aa)
    manifest["aa_floor"] = floor
    save()
    if floor["abs_mean_over_se"] is not None and floor["abs_mean_over_se"] >= 3:
        raise RuntimeError("resolved A/A drift; do not start A/B, reserve a new window")
    ab = []
    for i in range(12):
        order = ("A", "B") if i % 2 == 0 else ("B", "A")
        times = {arm: run(f"ab-{i}-{arm}", arm) for arm in order}
        ab.append((times["A"], times["B"]))
    result = paired_stats(ab)
    manifest["ab"] = result
    manifest["resolved_positive"] = (result["mean_saving_s"] > 0 and result["abs_mean_over_se"] is not None
                                      and result["abs_mean_over_se"] >= 3 and result["median_b_s"] < result["median_a_s"])
    save()
    # Not pooled with acceptance: 3 instrumented crossed pairs for mechanism.
    for i in range(3):
        for arm in (("A", "B") if i % 2 == 0 else ("B", "A")):
            run(f"diagnostic-{i}-{arm}", arm, profile=True)
    manifest["final_binary_sha256"] = sha(binary)
    assert manifest["final_binary_sha256"] == expected_hash
    assert sha(input_path) == manifest["input_sha256"]
    for name, path in config["pseudopotentials"].items():
        assert sha(input_path.parent / path) == manifest["upf_sha256"][name]
    save()
    print(json.dumps({"aa_floor": floor, "ab": result, "resolved_positive": manifest["resolved_positive"]}, indent=2))


if __name__ == "__main__":
    main()
