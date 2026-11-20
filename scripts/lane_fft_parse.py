#!/usr/bin/env python3
"""Parse LANE-FFT Phase-0 results; never infer a speedup from a profile.

Input: results.json with phase_profile.fft_lane.per_rank[].buckets[].
Each stage has seconds, calls, points, logical_bytes and fft_flops_estimate.
Times are inclusive of instrumentation overhead. Bytes are logical traffic,
not a DRAM performance counter. Do not sum rank times into a wall time.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

FIELDS = ["source", "model", "ranks", "rank", "bucket", "maps", "stage",
          "status", "calls", "seconds", "us_per_call", "ms_per_map",
          "points", "logical_bytes", "fft_flops_estimate"]


def rows(path):
    data = json.loads(Path(path).read_text())
    profile = data.get("phase_profile", {})
    lane = profile.get("fft_lane", {})
    if lane.get("schema") != "lane-fft-1" or not lane.get("enabled"):
        raise ValueError(f"{path}: FFT lane instrumentation unavailable or disabled")
    if data.get("threading", {}).get("threads", 1) > 1:
        raise ValueError(f"{path}: lane profile supports one thread per rank only")
    if lane.get("stage_timers") not in ("on", "off"):
        raise ValueError(f"{path}: unknown stage timer mode")
    ranks = profile.get("ranks")
    parts = lane.get("per_rank", [])
    if not isinstance(ranks, int) or ranks < 1 or len(parts) != ranks:
        raise ValueError(f"{path}: incomplete per-rank FFT records")
    if sorted(p.get("rank") for p in parts) != list(range(ranks)):
        raise ValueError(f"{path}: duplicate or missing ranks")
    timed = lane.get("stage_timers") == "on"
    for part in parts:
        for bucket in part["buckets"]:
            maps = bucket["maps"]
            for stage in bucket["stages"]:
                calls = stage["calls"]
                sec = stage["seconds"]
                values = [maps, calls, sec, stage["points"], stage["logical_bytes"], stage["fft_flops_estimate"]]
                if any(not isinstance(v, (int, float)) or not 0 <= v < float("inf") for v in values):
                    raise ValueError(f"{path}: invalid nonnegative finite metric")
                if not timed and sec != 0:
                    raise ValueError(f"{path}: nonzero timer in count-only mode")
                yield dict(zip(FIELDS, [str(path), data.get("parallel", {}).get("model", "single-process"),
                    ranks, part["rank"], bucket["bucket"], maps, stage["stage"],
                    "unusable-zero-calls" if not calls else ("timed" if timed else "count-only"),
                    calls, sec if timed else None,
                    sec * 1e6 / calls if timed and calls else None,
                    sec * 1e3 / maps if timed and maps else None,
                    stage["points"], stage["logical_bytes"], stage["fft_flops_estimate"]]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+")
    parser.add_argument("--json", action="store_true", help="emit records as JSON instead of CSV")
    args = parser.parse_args()
    try:
        records = [r for p in args.results for r in rows(p)]
        if args.json:
            json.dump(records, sys.stdout, indent=2, allow_nan=False)
            print()
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(records)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"lane_fft_parse: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
