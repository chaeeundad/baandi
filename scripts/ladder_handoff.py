#!/usr/bin/env python3
"""Offline handoff/fraction experiments; outputs stay under worktree target/scratch.

Build: cargo test --release -p dft-engine --lib --no-run
Run with the resulting dft_engine test executable (--binary). Fraction candidates
use explicit run-owned switches; defaults stay unchanged. Captures and accurate
references are deliberately outside timing claims.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import statistics
import struct
import shutil
from collections import defaultdict
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("fast_screen", ROOT / "scripts/fast_screen.py")
fast = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fast)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixture(case, seed, amp, pseudo):
    # si64 is shipped as si64-gamma.toml too; covered by ladder_2x2_tests.py.
    aliases = {f"si{n}": f"si{n}-gamma" for n in [8, 16, 32, 64]}
    pp = []
    if case.startswith("si"):
        c = fast.load_case(ROOT / "benchmarks/cases" / (aliases.get(case, case) + ".toml"))
        rows, atoms, cutoff = c["lattice"], c["atoms"], c["cutoff_ry"] / 2
        labels = ["Si"] * len(atoms)
        pp = [pseudo / "Si.pz-vbc.UPF"]
        extra = [f"pp Si {pp[0]}"]
    elif case == "mgo8":
        rows = [[7.956, 0, 0], [0, 7.956, 0], [0, 0, 7.956]]
        atoms, labels = [], []
        for s in [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]:
            atoms.extend([s, [(s[0] + .5) % 1, s[1], s[2]]])
            labels.extend(["Mg", "O"])
        pp = [pseudo / "pseudodojo/ex" / (s + ".upf") for s in ["Mg", "O"]]
        extra = [f"pp {s} {p}" for s, p in zip(["Mg", "O"], pp)] + ["pbe"]
        cutoff = 35
    elif case in ["syn-skew2", "syn-slab"]:
        rows = [[5, 0, 0], [.5, 4, 0], [.2, .4, 6]] if case == "syn-skew2" else [[4, 0, 0], [0, 4, 0], [0, 0, 10]]
        atoms = [[.13, .21, .32], [.6, .7, .8]] if case == "syn-skew2" else [[.1, .2, 0], [.6, .7, .2]]
        labels, cutoff = ["Xx", "Yy"], 4
        extra = ["synthetic Xx 2 1 nonlocal", "synthetic Yy 2 1.3 local", "bands 4"]
    else:
        raise ValueError(f"unknown case: {case}")
    atoms = fast.displaced(atoms, seed, amp)
    lines = ["row " + " ".join(map(repr, r)) for r in rows]
    lines += [f"atom {s} " + " ".join(map(repr, a)) for s, a in zip(labels, atoms)]
    lines += [f"cutoff {cutoff}", f"seed {seed}"] + extra
    return "\n".join(lines) + "\n", {str(p): sha(p) for p in pp}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--cases", default="si8,si16,si32,si-skew,mgo8,syn-skew2,syn-slab")
    p.add_argument("--seeds", default="601,602,603")
    p.add_argument("--amp", type=float, default=1e-3)
    p.add_argument("--arms", default=",".join(f"0.5:{s}:default" for s in ["1e-4", "1e-2", "3e-3", "1e-3", "3e-4", "3e-5"]))
    p.add_argument("--pseudo", type=Path, default=ROOT.parent / "pseudo")
    p.add_argument("--capture", action="store_true")
    p.add_argument("--trace-pulay", action="store_true")
    p.add_argument("--timeout", type=int, default=900)
    args = p.parse_args()
    scratch = (ROOT / "target/scratch").resolve()
    out = args.out.resolve()
    if not out.is_relative_to(scratch):
        p.error("--out must be inside this worktree's target/scratch")
    if out.exists() and any(out.iterdir()):
        p.error("--out must be empty; use a fresh directory to prevent stale evidence")
    out.mkdir(parents=True, exist_ok=True)
    source_binary = args.binary.resolve()
    # A later build must not replace the executable between geometry jobs.
    binary = out / "experiment"
    shutil.copy2(source_binary, binary)
    # Do not inherit algorithm overrides from the shell; the experiment is the
    # current default stack plus explicitly listed run-owned profiling controls.
    env = {k: v for k, v in os.environ.items() if not k.startswith("BAANDI_")}
    env.update(TMPDIR_990=str(scratch), TMPDIR=str(scratch), HANDOFF_ARMS=args.arms,
               HANDOFF_TRACE="on" if args.trace_pulay else "off",
               HANDOFF_CAPTURE="on" if args.capture else "off")
    provenance = {"binary": str(binary), "binary_sha256": sha(binary),
                  "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  "engine_source_sha256": sha(ROOT / "crates/dft-engine/src/cutoff_ladder.rs"),
                  "script_sha256": sha(__file__), "synthetic_source_sha256": sha(ROOT / "crates/dft-engine/tests/common/mod.rs"), "platform": platform.platform(),
                  "arms": args.arms, "capture": args.capture, "trace_pulay": args.trace_pulay, "amplitude_fractional": args.amp,
                  "load_start": os.getloadavg(), "threads": 1}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    all_rows = []
    failures = []
    for case in args.cases.split(","):
        for seed in map(int, args.seeds.split(",")):
            job = out / f"{case}-{seed}"
            job.mkdir(exist_ok=True)
            text, hashes = fixture(case, seed, args.amp, args.pseudo.resolve())
            inp = job / "fixture.txt"
            inp.write_text(text)
            metadata = {"case": case, "seed": seed, "input_sha256": sha(inp), "pp_sha256": hashes}
            (job / "input.json").write_text(json.dumps(metadata, indent=2) + "\n")
            job_env = dict(env, HANDOFF_INPUT=str(inp), HANDOFF_OUTPUT=str(job))
            if args.capture:
                job_env["BAANDI_SNAPSHOT_DIR"] = str(job / "snapshots")
            cmd = [str(binary), "scf::cutoff_ladder::handoff_experiment::measure", "--exact", "--ignored", "--nocapture", "--test-threads=1"]
            print(f"RUN {case} seed={seed}", flush=True)
            with (job / "run.log").open("w") as log:
                try:
                    rc = subprocess.run(cmd, cwd=ROOT, env=job_env, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout).returncode
                except subprocess.TimeoutExpired:
                    rc = 124
            measurements = sorted(job.glob("arm-*/measurement.json"))
            if rc == 0 and len(measurements) != len(args.arms.split(",")):
                rc = 2
            (job / "exit.json").write_text(json.dumps({"command": cmd, "exit_code": rc}) + "\n")
            for path in measurements:
                row = json.loads(path.read_text())
                row.update(metadata)
                ratio = math.prod(row["coarse_grid"]) / math.prod(row["fine_grid"])
                # Atomic starts are already in real H counters. The excess
                # wave transforms and kept forward-only columns are delta-V
                # updates, not additional full H applications.
                assert row["coarse_apply_h"] == row["coarse_h"] + row["coarse_start_h"]
                assert row["fine_apply_h"] == row["fine_history_h"]
                assert row["verification_apply_h"] == row["verification_h"]
                delta_calls = 0
                for stage in ["coarse", "fine", "verification"]:
                    full_delta = row[f"{stage}_wave_h"] - row[f"{stage}_apply_h"]
                    delta = full_delta + row[f"{stage}_kept_h"]
                    assert full_delta >= 0 and delta % row["bands"] == 0
                    row[f"{stage}_delta_full_h"] = full_delta
                    row[f"{stage}_delta_calls"] = delta // row["bands"]
                    delta_calls += delta // row["bands"]
                assert delta_calls == row["delta_apply_calls"]
                row["total_h"] = sum(row[f"{s}_apply_h"] for s in ["coarse", "fine", "verification"])
                row["total_calls"] = sum(row[f"{s}_apply_calls"] for s in ["coarse", "fine", "verification"])
                row["total_operator_calls"] = row["total_calls"] + delta_calls
                row["total_h_eq"] = ratio * row["coarse_apply_h"] + row["fine_apply_h"] + row["verification_apply_h"]
                row["with_probe_h_eq"] = row["total_h_eq"] + row["probe_h"]
                row["fft_column_eq"] = sum((ratio if stage == "coarse" else 1) *
                    (row[f"{stage}_wave_h"] + .5 * row[f"{stage}_kept_h"])
                    for stage in ["coarse", "fine", "verification"])
                all_rows.append(row)
            if rc:
                failures.append({"case": case, "seed": seed, "exit_code": rc})
                print(f"FAIL {case} seed={seed} exit={rc}; see {job / 'run.log'}", flush=True)
            else:
                print(f"DONE {case} seed={seed}", flush=True)
            (out / "rows.json").write_text(json.dumps(all_rows, indent=2) + "\n")
    (out / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    if all_rows:
        keys = sorted(set().union(*(row.keys() for row in all_rows)))
        with (out / "rows.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_rows)
    print(f"Rows: {len(all_rows)}, failed jobs: {len(failures)}", flush=True)
    return int(bool(failures))


def validate_captures(root):
    """Independent shape/header reader for the snapshot-1 binary layout."""
    count = 0
    for path in root.glob("*/*/snapshots/*.json"):
        m = json.loads(path.read_text())
        assert m["schema"] == "baandi-snapshot-1" and m["ranks"] == 1
        density = path.parent / m["density_file"]
        assert density.stat().st_size == 8 * math.prod(m["grid"])
        if m["start"] != "none":
            rank = path.parent / m["rank_files"].replace("<r>","0")
            with rank.open("rb") as f:
                header = f.read(48)
            assert header[:8] == b"BAANDSN1"
            nr,nc,hx,_,layout = struct.unpack("<5Q",header[8:])
            assert (nr,nc,layout)==(m["plane_waves"],m["bands"],2)
            assert bool(hx)==m["has_hx"]
            assert rank.stat().st_size==48+8*nr*nc*(1+hx)
        assert m["result"] is not None
        count += 1
    for path in root.glob("*/arm-*/measurement.json"):
        m = json.loads(path.read_text())
        for name,size in {
            "coarse-orbitals.bin":16*m["coarse_plane_waves"]*m["bands"],
            "handoff-orbitals.bin":16*m["fine_plane_waves"]*m["bands"],
            "coarse-density.bin":8*math.prod(m["coarse_grid"]),
            "handoff-density.bin":8*math.prod(m["fine_grid"]),
        }.items():
            assert (path.parent/name).stat().st_size == size
    return count


def analyze(argv):
    """Reproduce the registered seed split and emit reviewable tables."""
    parser = argparse.ArgumentParser(description=analyze.__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    out = args.out.resolve()
    if not out.is_relative_to((ROOT / "target/scratch").resolve()):
        parser.error("--out must be in target/scratch")
    out.mkdir(parents=True, exist_ok=True)
    rows = json.loads(args.rows.read_text())
    capture_count = validate_captures(args.rows.parent)
    groups = defaultdict(list)
    for row in rows:
        groups[(row["case"], row["seed"])].append(row)
    cases = sorted({r["case"] for r in rows})
    for rs in groups.values():
        longest = max(rs,key=lambda r:r["coarse_maps"])
        for row in rs:
            for key in ["coarse_tolerances","coarse_standard_errors"]:
                assert row[key] == longest[key][:row["coarse_maps"]]
    for case in cases:
        data = [r for r in rows if r["case"]==case]
        assert len({tuple(r["fine_grid"]) for r in data}) == 1
        assert len({r["fine_plane_waves"] for r in data}) == 1
        assert {r["seed"] for r in data} == {601,602,603}
        assert all(len(groups[(case,seed)])==6 for seed in [601,602,603])
    qs = [.001, .01, .1, 1.]
    stops = sorted({r["stop"] for r in rows})
    def baseline(rs):
        return next(r for r in rs if r["stop"] == 1e-4)
    def policy(rs, q):
        base = baseline(rs)
        candidates = {r["coarse_maps"]: r for r in rs
                      if r["stop"] in [1e-3, 3e-4] and r["coarse_maps"] < base["coarse_maps"]}
        probes = 0
        for _, row in sorted(candidates.items()):
            probes += 1
            if row["coarse_standard_errors"][-1] <= q * row["indicator_high_ha"][1]:
                return row, probes
        return base, probes
    def charge(row, probes):
        return row["fft_column_eq"] + probes * row["probe_h"]
    training = {}
    for q in qs:
        ratios = []
        for case in cases:
            rs = groups[(case, 601)]
            row, n = policy(rs, q)
            ratios.append(charge(row,n) / baseline(rs)["fft_column_eq"] if row["error_ha"] <= 5e-9 else float("inf"))
        training[q] = statistics.mean(ratios)
    selected = min(qs, key=lambda q: (training[q], q))
    fixed_training = {}
    for stop in stops:
        ratios = []
        for case in cases:
            rs = groups[(case, 601)]
            row = next(r for r in rs if r["stop"] == stop)
            ratios.append(row["fft_column_eq"] / baseline(rs)["fft_column_eq"] if row["error_ha"] <= 5e-9 else float("inf"))
        fixed_training[stop] = statistics.mean(ratios)
    fixed = min(stops, key=lambda stop: (fixed_training[stop], stop))
    verdicts = []
    for case in cases:
        measured = []
        for seed in [602,603]:
            rs = groups[(case,seed)]
            row, probes = policy(rs, selected)
            base = baseline(rs)
            fixed_row = next(r for r in rs if r["stop"] == fixed)
            measured.append({"seed":seed,"stop":row["stop"],"probes":probes,
                "baseline_cost":base["fft_column_eq"],"policy_cost":charge(row,probes),
                "fixed_cost":fixed_row["fft_column_eq"],"policy_error_ha":row["error_ha"],
                "ratio":charge(row,probes)/base["fft_column_eq"]})
        ratio = statistics.median(r["ratio"] for r in measured)
        verdicts.append({"case":case,"median_paired_ratio":ratio,
                         "win":ratio<1 and all(r["policy_error_ha"]<=5e-9 for r in measured),"rows":measured})
    calibration = []
    def corr(xs, ys):
        mx, my = statistics.mean(xs), statistics.mean(ys)
        covariance = sum((x-mx)*(y-my) for x,y in zip(xs,ys))
        scale = math.sqrt(sum((x-mx)**2 for x in xs)*sum((y-my)**2 for y in ys))
        return covariance/scale if scale else None
    for case in cases:
        unique = {(r["seed"],r["coarse_maps"]):r for r in rows if r["case"]==case}
        data = list(unique.values())
        high = [r["indicator_high_ha"][1] for r in data]
        low = [r["indicator_low_ha"][1] for r in data]
        fine = [sum(r[f"{stage}_wave_h"]+.5*r[f"{stage}_kept_h"]
                    for stage in ["fine","verification"]) for r in data]
        calibration.append({"case":case,"unique_states":len(data),
            "high_min_ha":min(high),"high_max_ha":max(high),
            "high_relative_span":(max(high)-min(high))/statistics.median(high),
            "remaining_work_min":min(fine),"remaining_work_max":max(fine),
            "high_vs_work_pearson":corr(high,fine),"low_vs_work_pearson":corr(low,fine)})
    analysis = {"rows_sha256":sha(args.rows),"validated_snapshot_manifests":capture_count,"q_training":training,"selected_q":selected,
                "fixed_training":fixed_training,"selected_fixed_stop":fixed,
                "held_out":verdicts,"winning_systems":sum(v["win"] for v in verdicts),
                "max_error_ha":max(r["error_ha"] for r in rows),"calibration":calibration}
    # Infinite training scores mean an accuracy failure, not valid JSON numbers.
    for scores in [analysis["q_training"],analysis["fixed_training"]]:
        for k,v in scores.items():
            if not math.isfinite(v): scores[k] = "accuracy failure"
    (out / "analysis.json").write_text(json.dumps(analysis,indent=2,allow_nan=False)+"\n")
    lines = ["### Complete candidate costs (per geometry)", "",
        "Hc includes atomic start; Hf includes redo; Hv includes every verification attempt. C is full H block calls; D is delta-potential block calls. W is the grid-weighted FFT-column proxy including delta updates, before probes. eta is shift=1 Ha. All final energies passed verification.", "",
        "| System | Seed | Stop | Coarse/fine/verify maps | Hc/Hf/Hv | Total H | C/D | W | eta (Ha) | Error (Ha) |",
        "| --- | ---: | ---: | --- | --- | ---: | --- | ---: | ---: | ---: |"]
    for r in sorted(rows,key=lambda r:(r["case"],r["seed"],-r["stop"])):
        lines.append(f"| {r['case']} | {r['seed']} | {r['stop']:.0e} | {r['coarse_maps']}/{r['fine_maps']}/{r['verification_maps']} | {r['coarse_apply_h']}/{r['fine_apply_h']}/{r['verification_apply_h']} | {r['total_h']} | {r['total_calls']}/{r['delta_apply_calls']} | {r['fft_column_eq']:.2f} | {r['indicator_high_ha'][1]:.6g} | {r['error_ha']:.2e} |")
    lines += ["", "### Held-out policy comparison", "",
        f"Training seed 601 selects q={selected:g}, fixed stop={fixed:g}; evaluation seeds 602/603. Probe columns are included in the policy cost.", "",
        "| System | Baseline W | Policy W | Fixed W | Paired policy/baseline | Win |",
        "| --- | ---: | ---: | ---: | ---: | --- |"]
    for v in verdicts:
        rs=v["rows"]
        med=lambda k:statistics.median(r[k] for r in rs)
        lines.append(f"| {v['case']} | {med('baseline_cost'):.2f} | {med('policy_cost'):.2f} | {med('fixed_cost'):.2f} | {v['median_paired_ratio']:.5f} | {v['win']} |")
    (out / "tables.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(analysis,indent=2,allow_nan=False))
    return 0


def analyze_config(argv):
    parser = argparse.ArgumentParser(description="Summarize coarse-configuration experiments")
    parser.add_argument("--rows",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args(argv)
    out=args.out.resolve()
    if not out.is_relative_to((ROOT/"target/scratch").resolve()):
        parser.error("--out must be in target/scratch")
    out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(args.rows.read_text())
    assert len(rows)==105
    validate_captures(args.rows.parent)
    fine_work=lambda r:sum(r[f"{s}_wave_h"]+.5*r[f"{s}_kept_h"] for s in ["fine","verification"])
    summaries=[]
    timing=[]
    for case in sorted({r["case"] for r in rows}):
        data=[r for r in rows if r["case"]==case]
        base={r["seed"]:r for r in data if r["arm"]=="0.5:1e-4:default"}
        assert set(base)=={601,602,603}
        assert len({tuple(r["fine_grid"]) for r in data})==1
        for arm in sorted({r["arm"] for r in data}):
            rs=[r for r in data if r["arm"]==arm]
            med=lambda key:statistics.median(r[key] for r in rs)
            summaries.append({"case":case,"arm":arm,"grid":rs[0]["coarse_grid"],
                "coarse_maps":med("coarse_maps"),"fine_maps":med("fine_maps"),"verify_maps":med("verification_maps"),
                "coarse_h":med("coarse_apply_h"),"fine_h":med("fine_apply_h"),"verify_h":med("verification_apply_h"),
                "total_h":med("total_h"),"total_calls":med("total_calls"),"delta_calls":med("delta_apply_calls"),
                "cost":med("fft_column_eq"),
                "paired_cost_ratio":statistics.median(r["fft_column_eq"]/base[r["seed"]]["fft_column_eq"] for r in rs),
                "paired_fine_ratio":statistics.median(fine_work(r)/fine_work(base[r["seed"]]) for r in rs),
                "seeds_with_more_fine_work":sum(fine_work(r)>fine_work(base[r["seed"]]) for r in rs),
                "observed_wall_ratio":statistics.median(r["wall_s"]/base[r["seed"]]["wall_s"] for r in rs),
                "max_error_ha":max(r["error_ha"] for r in rs)})
        rs=list(base.values())
        phases={key:statistics.median(r[key] for r in rs) for key in [
            "coarse_solve_s","coarse_density_s","coarse_potential_s","coarse_energy_s","coarse_mixing_s",
            "coarse_session_s","transfer_s","probe_apply_s","probe_s","probe_fields_s","probe_build_s","wall_s"]}
        timing.append({"case":case,"coarse_maps":[r["coarse_maps"] for r in rs],
            "last_tolerances_ha":[r["coarse_tolerances"][-1] for r in rs],
            "probe_h":rs[0]["probe_h"],**phases})
    result={"rows_sha256":sha(args.rows),"max_error_ha":max(r["error_ha"] for r in rows),
        "loose_controls":sum(r["loose_driver_identity"] for r in rows),"summaries":summaries,"timings":timing}
    (out/"config-analysis.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    lines=["### Coarse configuration screen (medians over seeds 601–603)","",
        "W includes all fine work and verification; ratios are medians of paired ratios. Fine ratio includes verification and delta updates. Wall ratios are loaded-Mac observations, not acceptance evidence. No capture I/O is in these walls.","",
        "| System | Arm (fraction:stop:mode) | Coarse grid | Coarse/fine/verify maps | Hc/Hf/Hv | Total H | C/D | W | W/base | Fine/base | Wall/base (O) |",
        "| --- | --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: |"]
    for r in summaries:
        lines.append(f"| {r['case']} | {r['arm']} | {'x'.join(map(str,r['grid']))} | {r['coarse_maps']}/{r['fine_maps']}/{r['verify_maps']} | {r['coarse_h']}/{r['fine_h']}/{r['verify_h']} | {r['total_h']} | {r['total_calls']}/{r['delta_calls']} | {r['cost']:.2f} | {r['paired_cost_ratio']:.4f} | {r['paired_fine_ratio']:.4f} | {r['observed_wall_ratio']:.3f} |")
    lines += ["", "### Coarse phase and probe observations (seconds)","",
        "Baseline, capture off, median of three different geometries. Probe total includes independent problem construction and all three shifts plus low-G diagnostics; H-only is the single occupied-block apply. Neither is a measured production-controller wall.","",
        "| System | Maps (601/602/603) | Solve | Density | Potential | Energy | Mixing | Session | Handoff | Probe H columns | Probe H-only | Probe total |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in timing:
        lines.append(f"| {r['case']} | {'/'.join(map(str,r['coarse_maps']))} | {r['coarse_solve_s']:.6f} | {r['coarse_density_s']:.6f} | {r['coarse_potential_s']:.6f} | {r['coarse_energy_s']:.6f} | {r['coarse_mixing_s']:.6f} | {r['coarse_session_s']:.6f} | {r['transfer_s']:.6f} | {r['probe_h']} | {r['probe_apply_s']:.6f} | {r['probe_s']:.6f} |")
    lines += ["", "### Coarse tolerance schedules (Ha)", "", "All baseline schedules start at the standard Ritz cap 0.005 Ha; later entries are the actual requested tolerances, not a fitted schedule.", ""]
    for r in rows:
        if r["arm"]=="0.5:1e-4:default":
            lines.append(f"- {r['case']} seed {r['seed']}: " + ", ".join(f"{t:.6g}" for t in r["coarse_tolerances"]))
    (out/"config-tables.md").write_text("\n".join(lines)+"\n")
    print(f"Validated {len(rows)} rows; {result['loose_controls']} exact loose-driver controls; max error {result['max_error_ha']:.6e} Ha")
    return 0


if __name__ == "__main__":
    command=sys.argv[1:2]
    if command==["analyze"]:
        sys.exit(analyze(sys.argv[2:]))
    if command==["analyze-config"]:
        sys.exit(analyze_config(sys.argv[2:]))
    sys.exit(main())
