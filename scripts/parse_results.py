#!/usr/bin/env python3
"""Validate a dft-cli results.json before its numbers are used (RUN-REPORT,
2026-11-18). Standard library only; run with `python3 -I`.

    python3 -I scripts/parse_results.py RESULTS.json [--manifest M.json |
        --binary PATH | --expect-sha HEX] [--json]

Exit 0 and one `ACCEPT` line when the file is usable, exit 2 and one `REJECT:`
line per problem otherwise (exit 1 on usage errors). Rejected:

* truncated or interrupted results: not valid JSON, or no completion marker
  (`run.completed == true`, written last by `dft-cli run`);
* stale binary: the binary SHA-256 recorded in `run.binary.sha256` differs from
  the expected one (`--manifest` with `binary_sha256`/`sha256`, `--binary` to
  hash a file, or `--expect-sha`), or is missing when one is expected;
* an active measured path whose clocks are missing or zero: the phase profile
  switch on but no `phase_profile` block; with stage timers on, a profile
  bucket that ran maps but has zero wall, a phase profile whose per-rank
  covered time is zero, Davidson stage timers on with every stage second zero
  after Davidson iterations, MPI collectives with calls but zero seconds;
* impossible rows: a profile entry whose minimum over ranks exceeds its
  maximum, a negative `solve.other` / `other_*` entry, a negative per-rank
  closure residual beyond clock granularity (covered time larger than the
  rank's own wall), a closure whose covered minimum exceeds its maximum.
"""
import argparse
import hashlib
import json
import sys

EPS = 1e-9


def _switch(doc, name):
    sw = ((doc.get("scf") or {}).get("switches") or {}).get(name) or {}
    return sw.get("value")


def _residual_tolerance(wall):
    # Clock granularity and float sums: 1 ms or 0.1 % of the wall.
    return max(1e-3, 1e-3 * (wall or 0.0))


def check_profile(profile, stage_timers_on):
    problems = []
    if profile.get("stage_timers") == "off":
        stage_timers_on = False
    for bucket in profile.get("buckets", []):
        label = bucket.get("bucket", "?")
        mx, mn = bucket.get("max", {}), bucket.get("min", {})
        for key, hi in mx.items():
            lo = mn.get(key)
            if lo is not None and hi is not None and lo > hi + EPS:
                problems.append(f"impossible row: phase_profile {label}.{key} min {lo} > max {hi}")
            if (key == "solve.other" or key.startswith("other_")) and hi is not None and hi < -EPS:
                problems.append(f"impossible row: phase_profile {label}.{key} negative ({hi})")
        if stage_timers_on and (mx.get("maps") or 0) > 0 and (mx.get("wall") or 0.0) <= 0.0:
            problems.append(f"clocks missing: phase_profile bucket {label} ran {mx.get('maps')} maps with zero wall")
    per_rank = profile.get("per_rank")
    if per_rank is not None:
        for r in per_rank:
            wall, residual = r.get("wall_s"), r.get("residual_s")
            if residual is not None and residual < -_residual_tolerance(wall):
                problems.append(
                    f"impossible row: rank {r.get('rank')} covered {r.get('covered_s')} s exceeds its wall {wall} s (residual {residual})"
                )
        if stage_timers_on and per_rank and all((r.get("covered_s") or 0.0) <= 0.0 for r in per_rank):
            problems.append("clocks missing: phase profile active but every rank's covered time is zero")
    closure = profile.get("closure")
    if closure:
        if closure.get("covered_min_s", 0.0) > closure.get("covered_max_s", 0.0) + EPS:
            problems.append("impossible row: closure covered_min > covered_max")
        rmin, rmax = closure.get("residual_min_s"), closure.get("residual_max_s")
        if rmin is not None and rmax is not None and rmin > rmax + EPS:
            problems.append("impossible row: closure residual_min > residual_max")
    return problems


def check(doc, expected_sha=None):
    """Problems of a parsed results document (empty list: accepted)."""
    problems = []
    if not isinstance(doc, dict):
        return ["truncated: top level is not an object"]
    run = doc.get("run")
    if not isinstance(run, dict) or run.get("completed") is not True:
        problems.append("truncated or interrupted: no completion marker (run.completed)")
        run = run if isinstance(run, dict) else {}
    recorded = (run.get("binary") or {}).get("sha256")
    if expected_sha is not None:
        if recorded is None:
            problems.append("stale binary: results carry no binary SHA-256 to compare with the manifest")
        elif recorded.lower() != expected_sha.lower():
            problems.append(f"stale binary: results from {recorded}, manifest expects {expected_sha}")
    stage_timers_on = _switch(doc, "BAANDI_STAGE_TIMERS") != "off"
    profile = doc.get("phase_profile")
    if _switch(doc, "BAANDI_PHASE_PROFILE") == "on" and profile is None:
        problems.append("clocks missing: BAANDI_PHASE_PROFILE is on but results have no phase_profile")
    if profile is not None:
        problems += check_profile(profile, stage_timers_on)
    work = (doc.get("convergence") or {}).get("work") or {}
    stages = work.get("davidson_stages")
    iterations = work.get("davidson_iterations_recorded") or 0
    if _switch(doc, "BAANDI_DAVIDSON_STAGES") == "on" and stage_timers_on and iterations > 0:
        stage_s = (stages or {}).get("stage_s")
        if not isinstance(stage_s, dict):
            problems.append("clocks missing: BAANDI_DAVIDSON_STAGES is on but davidson_stages.stage_s is absent")
        elif all((v or 0.0) <= 0.0 for v in stage_s.values()):
            problems.append(f"clocks missing: BAANDI_DAVIDSON_STAGES is on, {iterations} Davidson iterations, every stage second zero")
    parallel = doc.get("parallel") or {}
    if parallel.get("stage_timers") == "on":
        for c in parallel.get("collectives") or []:
            if (c.get("calls") or 0) > 0 and (c.get("seconds") or 0.0) <= 0.0:
                problems.append(f"clocks missing: collective {c.get('op')} has {c.get('calls')} calls and zero seconds")
    return problems


def expected_sha_from(args):
    if args.expect_sha:
        return args.expect_sha
    if args.binary:
        with open(args.binary, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()
    if args.manifest:
        with open(args.manifest) as f:
            m = json.load(f)
        for key in ("binary_sha256", "sha256"):
            if isinstance(m.get(key), str):
                return m[key]
        raise SystemExit(f"error: {args.manifest} has no binary_sha256 or sha256")
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--manifest", help="provenance manifest with binary_sha256")
    g.add_argument("--binary", help="binary whose SHA-256 the results must carry")
    g.add_argument("--expect-sha", help="expected binary SHA-256")
    ap.add_argument("--json", action="store_true", help="print a JSON verdict")
    args = ap.parse_args(argv)
    expected = expected_sha_from(args)
    try:
        with open(args.results, "rb") as f:
            raw = f.read()
        doc = json.loads(raw)
    except (OSError, ValueError) as e:
        problems, doc = [f"truncated: not valid JSON ({e.__class__.__name__}: {e})"], None
    else:
        problems = check(doc, expected)
    run = (doc or {}).get("run") or {} if isinstance(doc, dict) else {}
    if args.json:
        print(json.dumps({"results": args.results, "accepted": not problems, "problems": problems,
                          "run_id": run.get("run_id"), "exit_reason": run.get("exit_reason")}))
    elif problems:
        for p in problems:
            print(f"REJECT: {p}")
    else:
        closure = ((doc.get("phase_profile") or {}).get("closure")) or {}
        extra = ""
        if closure:
            extra = (f" covered_max={closure.get('covered_max_s'):.6g}s imbalance={closure.get('imbalance'):.3g}"
                     f" residual_max={closure.get('residual_max_s')}")
        print(f"ACCEPT {args.results} run_id={run.get('run_id')} exit_reason={run.get('exit_reason')}{extra}")
    return 2 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
