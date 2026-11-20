#!/usr/bin/env python3
"""Correctness-only instrumentation check. Deliberately records no timings."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k not in ("timings", "wall_s")}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def float_bits(value):
    import struct
    if isinstance(value, float):
        return struct.pack('!d', value).hex()
    if isinstance(value, dict):
        return {k: float_bits(v) for k, v in value.items()}
    if isinstance(value, list):
        return [float_bits(v) for v in value]
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--baseline", type=Path)
    ap.add_argument("--allow-nonconverged", action="store_true")
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--json-out", type=Path, required=True)
    ap.add_argument("inputs", type=Path, nargs="+")
    a = ap.parse_args()
    records = []
    for src in a.inputs:
        reference = None
        counts_reference = None
        modes = (["baseline"] if a.baseline else []) + ["off", "counts", "timed", "phase_off"]
        for mode in modes:
            binary = a.baseline if mode == "baseline" else a.binary
            out = (a.work / src.stem / mode).resolve()
            env = dict(os.environ, BAANDI_GLUE_PROFILE="off" if mode in ("off", "baseline") else "on",
                       BAANDI_FFT_PROFILE="on", BAANDI_PHASE_PROFILE="off" if mode == "phase_off" else "on", BAANDI_STAGE_TIMERS="off" if mode == "counts" else "on")
            cmd = [str(binary.resolve()), "run", str(src.resolve()), "--output-dir", str(out)]
            p = subprocess.run(cmd, env=env, capture_output=True, text=True)
            out.mkdir(parents=True, exist_ok=True)
            (out / "stdout.log").write_text(p.stdout)
            (out / "stderr.log").write_text(p.stderr)
            assert p.returncode in ((0, 3) if a.allow_nonconverged else (0,)), (cmd, p.returncode, p.stderr[-1000:])
            r = json.loads((out / "results.json").read_text())
            signature = dict(exit_code=p.returncode, converged=r["convergence"]["converged"], energy=r["total_energy"], density=digest(out / "density.bin"),
                             history=scientific(r["convergence"]["history"]),
                             eigenvalues=r.get("eigenvalues"), kpoints=r.get("kpoints"),
                             energies=r.get("energies"), work=r.get("work"), occupations=r.get("occupations"), verified=r.get("verified"),
                             verification=scientific(r.get("verification")))
            signature = float_bits(signature)
            if reference is None:
                reference = signature
            if not a.allow_nonconverged:
                assert r["verified"] and r["convergence"]["converged"]
            assert signature == reference, f"{src} {mode}: scientific output differs"
            lane = {k: v for k, v in r["run"]["counters"].items() if k.startswith("glue_lane.")}
            work = {k: v for k, v in lane.items() if not k.endswith(".ns")}
            timers = [v for k, v in lane.items() if k.endswith(".ns")]
            if mode in ("off", "phase_off"):
                assert lane and not any(lane.values()), "disabled instrumentation recorded work"
            elif mode in ("counts", "timed"):
                assert work["glue_lane.serial.mix.calls"] > 0
                if counts_reference is None:
                    counts_reference = work
                assert counts_reference == work, "timers changed work accounting"
                assert any(timers) == (mode == "timed"), "clock gating failed"
            records.append(dict(input=str(src), input_sha256=digest(src), mode=mode,
                                command=cmd, exit_code=p.returncode, energy_ha=r["total_energy"]["ha"],
                                delta_energy_ha=0.0, density_sha256=signature["density"],
                                species=r.get("species"), history_sha256=hashlib.sha256(json.dumps(signature["history"], sort_keys=True).encode()).hexdigest(),
                                scientific_bits_equal=True, converged=r["convergence"]["converged"],
                                verified=r.get("verified"), work=work))
            print(f"{src.stem} {mode}: exit {p.returncode}; scientific bits equal", flush=True)
    a.json_out.write_text(json.dumps(dict(binary_sha256=digest(a.binary),
                                         baseline_sha256=digest(a.baseline) if a.baseline else None,
                                         purpose="Mac correctness only; no performance conclusions",
                                         records=records), indent=2) + "\n")


if __name__ == "__main__":
    main()
