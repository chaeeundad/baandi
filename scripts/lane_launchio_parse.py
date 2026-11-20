#!/usr/bin/env python3
"""Validate rank-local Phase 0 output and emit LAUNCHIO tables (stdlib only)."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def number(value, label, nonnegative=True):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{label}: expected finite number")
    if nonnegative and value < 0:
        raise ValueError(f"{label}: negative value")
    return value


def load(path):
    def unique(pairs):
        out = {}
        for key, value in pairs:
            if key in out:
                raise ValueError(f"{path}: duplicate JSON key {key}")
            out[key] = value
        return out
    return json.loads(path.read_text(), object_pairs_hook=unique)


def parse_run(root):
    root = Path(root)
    if root.is_dir():
        root /= "results.json"
    paths = [root, *sorted(root.parent.glob("rank*/results.json"))]
    rows = []
    identities = set()
    expected = None
    for path in paths:
        doc = load(path)
        if doc["run"]["completed"] is not True or doc["convergence"]["converged"] is not True:
            raise ValueError(f"{path}: incomplete or unconverged run")
        if doc.get("verified") is not True:
            raise ValueError(f"{path}: fixed-point verification required")
        launch, lane = doc["launch"], doc["launchio"]
        if launch["schema"] != "launch-1" or lane["schema"] != "launchio-1":
            raise ValueError(f"{path}: unsupported schema")
        rank, ranks = launch["rank"], launch["ranks"]
        if type(rank) is not int or type(ranks) is not int or ranks < 1 or not 0 <= rank < ranks:
            raise ValueError(f"{path}: invalid rank")
        if rank != lane["rank"] or (expected is not None and expected != ranks):
            raise ValueError(f"{path}: rank mismatch")
        expected = ranks
        if ranks > 1 and (launch["mpi_backend"] != "mpi" or doc["parallel"]["model"] == "mpi_replicated"):
            raise ValueError(f"{path}: not a distributed MPI run")
        identities.add((doc["run"]["run_id"], doc["input_sha256"], doc["run"]["binary"]["sha256"]))
        density = path.with_name("density.bin")
        density_hash = hashlib.sha256(density.read_bytes()).hexdigest()
        if density_hash != load(path.with_name("density.json"))["sha256"]:
            raise ValueError(f"{path}: density hash mismatch")
        sizes = lane["sizes"]
        for key, value in sizes.items():
            if type(value) is not int or value < 0:
                raise ValueError(f"{path}: invalid size {key}")
        if density.stat().st_size != sizes["density_bytes"] or sizes["density_bytes"] != 8 * sizes["grid_points"]:
            raise ValueError(f"{path}: density byte count mismatch")
        seconds = {key: number(value, key) for key, value in launch["seconds"].items()}
        for key in ("scf_driver", "input_read", "input_parse", "density_encode", "density_sha256", "history_json", "results_serialize", "results_write", "mpi_init", "mpi_finalize"):
            if key not in seconds:
                raise ValueError(f"{path}: missing timer {key}")
        for key in ("scf_driver", "input_read", "input_parse", "results_serialize", "results_write"):
            if seconds[key] == 0:
                raise ValueError(f"{path}: unusable zero timer {key}")
        main = number(launch["main_wall_s"], "main_wall_s")
        covered = number(launch["covered_s"], "covered_s")
        closure = sum(number(v, k, False) for k, v in launch["closure"].items())
        outside = number(launch["outside_covered_s"], "outside_covered_s", False)
        if not math.isclose(main - covered, outside, abs_tol=1e-9) or not math.isclose(closure, outside, abs_tol=1e-9):
            raise ValueError(f"{path}: launch closure mismatch")
        buckets = {}
        for b in lane["buckets"]:
            label = b["bucket"]
            if label in buckets:
                raise ValueError(f"{path}: duplicate bucket {label}")
            times = {key: number(value, key) for key, value in b["seconds"].items()}
            maps = b["counts"]["maps"]
            if type(maps) is not int or maps < 0:
                raise ValueError(f"{path}: invalid map count")
            buckets[label] = {"seconds": times, "maps": maps,
                              "wall_s_per_map": times["wall"] / maps if maps else None}
        if set(buckets) != {"init", "coarse", "fine", "fine_redo", "verify_pass", "verify_fail", "materialize", "between"}:
            raise ValueError(f"{path}: missing or unexpected phase bucket")
        if not math.isclose(sum(b["seconds"]["wall"] for b in buckets.values()), covered, abs_tol=1e-9):
            raise ValueError(f"{path}: local bucket closure mismatch")
        for record in lane["upf_records"]:
            for key in ("bytes", "projectors", "read_s", "parse_s", "hash_s"):
                number(record[key], key)
        rows.append({"rank": rank, "file": str(path), "seconds": seconds,
                     "main_wall_s": main, "covered_s": covered,
                     "outside_covered_s": outside, "main_minus_scf_s": main - seconds["scf_driver"],
                     "sizes": sizes, "upf_records": lane["upf_records"], "buckets": buckets,
                     "results_bytes_including_launch_append": path.stat().st_size,
                     "energy_ha": number(doc["total_energy"]["ha"], "energy", False),
                     "density_sha256": density_hash})
    if sorted(row["rank"] for row in rows) != list(range(expected)) or rows[0]["rank"] != 0:
        raise ValueError(f"{root}: missing/duplicate ranks or root is not rank zero")
    if len(identities) != 1:
        raise ValueError(f"{root}: mixed run IDs, inputs or binaries")
    if len({(float(row["energy_ha"]).hex(), row["density_sha256"]) for row in rows}) != 1:
        raise ValueError(f"{root}: ranks disagree on delivered science")
    for row in rows:
        resource_path = root.parent / f"resource-rank-{row['rank']}.json"
        row["resources"] = None
        if resource_path.exists():
            resource = load(resource_path)
            if resource["exit_code"] != 0:
                raise ValueError(f"{resource_path}: failed rank")
            for key, value in resource.items():
                number(value, key)
            row["resources"] = resource
    keys = set.intersection(*(set(row["seconds"]) for row in rows))
    ranges = {key: {"min_s": min(r["seconds"][key] for r in rows),
                    "max_s": max(r["seconds"][key] for r in rows)} for key in sorted(keys)}
    external = None
    wall_path = root.parent / "launch-wall.json"
    if wall_path.exists():
        wall = load(wall_path)
        if wall["exit_code"] != 0:
            raise ValueError(f"{wall_path}: failed external launch")
        external = {"wall_s": number(wall["wall_s"], "external wall"),
                    "wall_minus_max_main_s": wall["wall_s"] - max(r["main_wall_s"] for r in rows),
                    "wall_minus_max_scf_s": wall["wall_s"] - max(r["seconds"]["scf_driver"] for r in rows)}
    return {"root": str(root), "ranks": expected, "identity": list(next(iter(identities))),
            "per_rank": rows, "timer_ranges": ranges, "external": external}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, nargs="+", help="one root results.json or directory per launch")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-external", action="store_true", help="require dev08 launch and every rank's resource sidecars")
    args = parser.parse_args()
    try:
        report = {"schema": "launchio-tables-1", "runs": [parse_run(p) for p in args.results],
                  "note": "Seconds are nested diagnostic clocks; never sum timer maxima across ranks. Per-map values are bucket means. External minus max local duration is a residual, not an isolated MPI/loader measurement. No speedup or adoption verdict."}
        if args.require_external and any(run["external"] is None or any(row["resources"] is None for row in run["per_rank"]) for run in report["runs"]):
            raise ValueError("missing launch-wall.json or resource-rank-N.json")
    except (KeyError, TypeError, ValueError, OSError) as exc:
        parser.exit(2, f"invalid LAUNCHIO evidence: {exc}\n")
    data = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
