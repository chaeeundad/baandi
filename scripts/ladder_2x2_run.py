#!/usr/bin/env python3
"""Captured-state 2x2 experiment (one core, release library test binary).

No production switch. All raw evidence is restricted to target/scratch. Each
geometry has its own strict reference, three unsplit natural controls, two exact
handoff states, four crossed continuations, a 0.45 cross-check, and an A/A repeat.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("handoff", ROOT / "scripts/ladder_handoff.py")
handoff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(handoff)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cases", default="si32,si64")
    parser.add_argument("--seeds", default="601,602,603")
    parser.add_argument("--pairs", type=int, default=1)
    parser.add_argument("--amp", type=float, default=1e-3)
    parser.add_argument("--pseudo", type=Path, default=ROOT.parent / "pseudo")
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    scratch = (ROOT / "target/scratch").resolve()
    out = args.out.resolve()
    if not out.is_relative_to(scratch) or args.pairs < 1:
        parser.error("--out must be in target/scratch and --pairs must be positive")
    if out.exists() and any(out.iterdir()):
        parser.error("--out must be empty")
    out.mkdir(parents=True, exist_ok=True)
    binary = out / "experiment"
    shutil.copy2(args.binary.resolve(), binary)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("BAANDI_", "HANDOFF_"))}
    env.update(TMPDIR_990=str(scratch), TMPDIR=str(scratch), HANDOFF_PAIRS=str(args.pairs),
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    provenance = dict(command=sys.argv, binary_sha256=handoff.sha(binary),
        source_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        source_sha256=handoff.sha(ROOT / "crates/dft-engine/src/cutoff_ladder.rs"),
        runner_sha256=handoff.sha(__file__), platform=platform.platform(),
        load_start=os.getloadavg(), threads=1, pairs=args.pairs,
        amplitude_fractional=args.amp)
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    failures = []
    for case in args.cases.split(","):
        for seed in map(int, args.seeds.split(",")):
            job = out / f"{case}-{seed}"
            job.mkdir()
            fixture, hashes = handoff.fixture(case, seed, args.amp, args.pseudo.resolve())
            inp = job / "fixture.txt"
            inp.write_text(fixture)
            (job / "input.json").write_text(json.dumps(dict(case=case, seed=seed,
                input_sha256=handoff.sha(inp), pp_sha256=hashes), indent=2) + "\n")
            cmd = [str(binary), "scf::cutoff_ladder::handoff_experiment::measure_2x2",
                   "--exact", "--ignored", "--nocapture", "--test-threads=1"]
            print(f"RUN {case}/{seed}", flush=True)
            with (job / "run.log").open("w") as log:
                try:
                    rc = subprocess.run(cmd, cwd=ROOT,
                        env=dict(env, HANDOFF_INPUT=str(inp), HANDOFF_OUTPUT=str(job)),
                        stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout).returncode
                except subprocess.TimeoutExpired:
                    rc = 124
            if rc == 0 and len(list(job.glob("pair-*/measurement.json"))) != 6 * args.pairs:
                rc = 2
            (job / "exit.json").write_text(json.dumps(dict(command=cmd, exit_code=rc)) + "\n")
            print(f"{'FAIL' if rc else 'DONE'} {case}/{seed}, exit={rc}", flush=True)
            if rc:
                failures.append(dict(case=case, seed=seed, exit_code=rc))
    (out / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
