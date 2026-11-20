#!/usr/bin/env python3
"""Parse LANE-INIT Phase 0 results, preserving rank and invocation boundaries.

No timings across jobs or ranks are added. Existing phase-profile min/max
are aggregated envelopes, not this file's local times. Lane counters are
local, cumulative per process/run; ns/call is an inclusive mean, not a trace.
"""
import argparse
import json
from pathlib import Path

PREFIX = "init_lane."


def parse(path):
    report = json.loads(path.read_text())
    run = report.get("run", {})
    counters = run.get("counters", {})
    lane = {k[len(PREFIX):]: v for k, v in counters.items() if k.startswith(PREFIX)}
    if not lane:
        raise ValueError(f"{path}: missing init_lane counters (wrong binary/schema)")
    if any(not isinstance(v, int) or isinstance(v, bool) or v < 0 for v in lane.values()):
        raise ValueError(f"{path}: lane counters must be nonnegative integers")
    calls = {k[:-6]: v for k, v in lane.items() if k.endswith(".calls")}
    entered = any(calls.values())
    timed = any(v for k, v in lane.items() if k.endswith(".ns"))
    if timed and not entered:
        raise ValueError(f"{path}: timers without calls")
    spans = []
    for name, count in calls.items():
        ns = lane.get(name + ".ns")
        if ns is None:
            raise ValueError(f"{path}: missing timer for {name}")
        spans.append(dict(span=name, calls=count,
                          seconds=ns / 1e9 if timed else None,
                          seconds_per_call=ns / 1e9 / count if timed and count else None))
    profile = report.get("phase_profile") or report.get("parallel", {}).get("phases") or {}
    envelope = {}
    for bucket in profile.get("buckets", []):
        if bucket.get("bucket") == "init":
            envelope = {bound: {k: v for k, v in bucket.get(bound, {}).items()
                               if k.startswith(("init.", "atomic.")) or k in ("wall", "atomic_start")}
                        for bound in ("min", "max")}
    session_envelope = {}
    for bucket in profile.get("buckets", []):
        if bucket.get("bucket") == "between":
            session_envelope = {bound: bucket.get(bound, {}).get("between.session_setup")
                                for bound in ("min", "max")}
    work = {k: v for k, v in lane.items() if not k.endswith((".calls", ".ns"))}
    # Logical payload floors, not measured traffic: no claim about cache reads,
    # write allocation, allocator metadata, FFT plan storage or effective GB/s.
    payload = {
        "factor_write_bytes": work.get("factors.bytes"),
        "projector_write_bytes": work.get("projectors.materialized_bytes"),
        "one_real_grid_bytes_per_problem_sum": 8 * work.get("problem.grid_points", 0),
        "one_complex_grid_bytes_per_problem_sum": 16 * work.get("problem.grid_points", 0),
    }
    launch = report.get("launch", {})
    return dict(file=str(path), run_id=run.get("run_id"), rank=run.get("rank"),
                ranks=run.get("ranks"), source_commit=run.get("source_commit"),
                binary_sha256=run.get("binary", {}).get("sha256"),
                input_sha256=report.get("input_sha256"),
                state="timed" if timed else "counts_only" if entered else "disabled_or_not_reached",
                spans=spans, work=work, payload_floors=payload,
                phase_profile_envelope=envelope,
                session_setup_envelope=session_envelope,
                atomic_start=report.get("atomic_start"), effective_start=report.get("effective_start"),
                final_shapes={k: report.get(k) for k in ("basis", "bands", "threading")},
                launch_seconds=launch.get("seconds", {}),
                existing_counters={k: v for k, v in counters.items()
                                   if k.startswith(("owned.", "phase_validate.")) or k in ("ewald_hits", "ewald_computed", "problems")},
                correctness=dict(energy_ha=report.get("total_energy", {}).get("ha"),
                                 converged=report.get("convergence", {}).get("converged"),
                                 verified=report.get("verified")))


def markdown(data):
    out = ["# INIT lane profile", "", data["label"], "",
           "Inclusive local timers; nested spans must not be added. A dash means unavailable.",
           "Counts and sizes are cumulative over all constructors (including ladder stages).", ""]
    for r in data["runs"]:
        out += [f"## {r['file']}", "", f"Rank {r['rank']}/{r['ranks']}; {r['state']}; run {r['run_id']}.", "",
                "| span | calls | local seconds | seconds/call |", "| --- | ---: | ---: | ---: |"]
        for s in r["spans"]:
            fmt = lambda v: "-" if v is None else f"{v:.9g}"
            out.append(f"| {s['span']} | {s['calls']} | {fmt(s['seconds'])} | {fmt(s['seconds_per_call'])} |")
        out += ["", "| work / size | count |", "| --- | ---: |"]
        out += [f"| {k} | {v} |" for k, v in r["work"].items()]
        out += ["", "| logical payload floor | bytes |", "| --- | ---: |"]
        out += [f"| {k} | {v} |" for k, v in r["payload_floors"].items()]
        out += ["", "Existing atomic/startup phase envelopes (already reduced min/max, seconds):", "",
                "| component | min | max |", "| --- | ---: | ---: |"]
        env = r["phase_profile_envelope"]
        for k in sorted(set(env.get("min", {})) | set(env.get("max", {}))):
            out.append(f"| {k} | {env.get('min', {}).get(k, '-')} | {env.get('max', {}).get(k, '-')} |")
        out.append("")
        out += [f"Session setup envelope (all sessions, not only first): {r['session_setup_envelope']}",
                f"Effective starts: {json.dumps(r['effective_start'], sort_keys=True)}", ""]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("paths", type=Path, nargs="+")
    ap.add_argument("--json-out", type=Path, required=True)
    ap.add_argument("--markdown-out", type=Path, required=True)
    ap.add_argument("--label", default="dev08 Phase 0; diagnostic instrumentation, not adoption timing.")
    args = ap.parse_args()
    files = sorted({p.resolve() for arg in args.paths
                    for p in (arg.rglob("results.json") if arg.is_dir() else [arg])})
    if not files:
        ap.error("no results.json files found")
    try:
        data = dict(schema="lane-init-profile-1", label=args.label, runs=[parse(p) for p in files])
    except (ValueError, KeyError, TypeError) as e:
        ap.error(str(e))
    args.json_out.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    args.markdown_out.write_text(markdown(data) + "\n")
    print(f"Parsed {len(files)} result files; rank boundaries preserved.")


if __name__ == "__main__":
    main()
