#!/usr/bin/env python3
"""Correctness only; compare every JSON field after explicit volatile-field exclusions."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized(value, path=()):
    """Keep all non-clock fields, including work counts and full scientific history."""
    if isinstance(value, dict):
        out = {}
        for k, v in value.items():
            full = path + (k,)
            if k == "timings" or k.endswith("_s"):
                continue
            if full in {("convergence", "work", "cutoff_ladder", "coarse_wall_fine_step_equivalents"), ("timestamp",), ("input_sha256",),
                        ("run", "run_id"), ("run", "switch_agreement_hash"),
                        ("run", "binary", "path"), ("run", "binary", "sha256"),
                        ("run", "binary", "hash_cache"),
                        ("scf", "switches", "BAANDI_BINARY_HASH_CACHE"),
                        ("execution", "non_default_switches", "BAANDI_BINARY_HASH_CACHE")}:
                continue
            out[k] = normalized(v, full)
        return out
    if isinstance(value, list):
        return [normalized(v, path + (str(i),)) for i, v in enumerate(value)]
    if isinstance(value, float):
        return {"f64_bits": struct.pack("!d", value).hex()}
    return value


def differences(a, b, prefix=""):
    if isinstance(a, dict) and isinstance(b, dict):
        return [item for k in a.keys() | b.keys() for item in (
            [prefix + "/" + k] if k not in a or k not in b else differences(a[k], b[k], prefix + "/" + k))]
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        return [item for i, (x, y) in enumerate(zip(a, b)) for item in differences(x, y, prefix + f"/{i}")]
    return [] if a == b else [prefix]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--upf", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    candidate = root / "candidate"
    shutil.copy2(args.candidate, candidate)
    cache = Path(str(candidate) + ".baandi-sha256.json")
    template = Path("examples/si_lda.toml").read_text().replace("/absolute/path/to/Si.UPF", str(args.upf.resolve()))
    cases = {
        "gamma": template,
        "skew": template.replace("[0.0, 2.715, 2.715]", "[0.08, 2.715, 2.715]"),
        "anisotropic": template.replace("[2.715, 2.715, 0.0]", "[2.79645, 2.79645, 0.0]"),
        "k2": template.replace("kpoints = [[0.0, 0.0, 0.0]]", "kpoints = [[0.0, 0.0, 0.0, 0.5], [0.25, 0.0, 0.0, 0.5]]"),
    }
    records = []
    for case, text in cases.items():
        if cache.exists():
            cache.unlink() # only this harness's sidecar next to its private binary
        oracle = None
        for mode in ("baseline", "default", "off", "cold", "warm", "malformed", "input", "unavailable"):
            binary = args.baseline.resolve() if mode == "baseline" else candidate
            config = text + "\n[scf]\nthreads = 1\n"
            if mode == "input":
                config += "\n[performance]\nbinary_hash_cache = true\n"
            input_path = root / f"{case}-{mode}.toml"
            input_path.write_text(config)
            env = {k: v for k, v in os.environ.items() if not k.startswith("BAANDI_")}
            env.update(BAANDI_PHASE_PROFILE="off", BAANDI_STAGE_TIMERS="off")
            if mode not in ("baseline", "default", "input"):
                env["BAANDI_BINARY_HASH_CACHE"] = "off" if mode == "off" else "on"
            if mode == "malformed":
                cache.write_text('{"broken":')
            if mode == "unavailable":
                cache.unlink()
                cache.mkdir() # publication failure must not remove provenance
            output = root / f"{case}-{mode}"
            command = [str(binary), "run", str(input_path), "--output-dir", str(output)]
            with (root / f"{case}-{mode}.log").open("w") as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError(f"{case}/{mode}: exit {result.returncode}")
            doc = json.loads((output / "results.json").read_text())
            assert doc["verified"] and doc["convergence"]["converged"]
            assert doc["run"]["binary"]["sha256"] == sha(binary), "provenance differs from independent hashlib"
            norm = normalized(doc)
            normalized_bytes = json.dumps(norm, sort_keys=True, separators=(",", ":")).encode()
            files = {name: sha(output / name) for name in ("density.bin", "density.json")}
            if oracle is None:
                oracle = norm, normalized_bytes, files, {name: (output / name).read_bytes() for name in files}
            assert normalized_bytes == oracle[1], f"{case}/{mode}: JSON changes {differences(oracle[0], norm)[:20]}"
            assert files == oracle[2], f"{case}/{mode}: density file hashes changed"
            assert all((output / name).read_bytes() == oracle[3][name] for name in files), f"{case}/{mode}: density file bytes changed"
            cache_record = doc["run"]["binary"].get("hash_cache")
            if mode in ("baseline", "default", "off"):
                assert not cache.exists() and cache_record is None
            else:
                expected = "hit" if mode in ("warm", "input") else "miss"
                assert cache_record["status"] == expected
                assert cache_record["hashed_bytes"] == (0 if expected == "hit" else candidate.stat().st_size)
                if mode == "unavailable":
                    assert cache.is_dir() and not cache_record["published"]
                    cache.rmdir()
            records.append({"case": case, "mode": mode, "exit_code": result.returncode,
                            "energy_ha": doc["total_energy"]["ha"], "delta_energy_ha": 0,
                            "density_sha256": files["density.bin"],
                            "normalized_results_sha256": hashlib.sha256(normalized_bytes).hexdigest(),
                            "hash_cache": cache_record, "provenance_verified": True})
    report = {"scope": "Mac correctness only; full normalized results and raw density bytes; no performance measurements",
              "baseline_sha256": sha(args.baseline), "candidate_sha256": sha(candidate),
              "upf_sha256": sha(args.upf), "runs": records}
    (root / "equivalence.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {len(records)} executions; all density bytes/normalized JSON identical; independent binary digests match")


if __name__ == "__main__":
    main()
