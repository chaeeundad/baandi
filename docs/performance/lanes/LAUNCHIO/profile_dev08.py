#!/usr/bin/env python3
"""Phase 1 launcher. Master only, exclusive dev08 window, one case per job."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", required=True, choices=["si16", "si32", "si64"])
    p.add_argument("--ranks", type=int, choices=[1, 8], required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if platform.system() != "Linux" or platform.node().split(".")[0] != "dev08":
        p.error("profiling is restricted to dev08 through the master session")
    if args.case == "si16" and args.ranks == 1:
        p.error("si16 one-core is outside this request")
    # Do not silently change an inherited algorithm switch or a production input.
    env = os.environ.copy()
    env.update(MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary, input_path = args.binary.resolve(), args.input.resolve()
    import tomllib
    config = tomllib.loads(input_path.read_text())
    if config.get("scf", {}).get("threads", 1) != 1:
        p.error("one thread per rank required")
    if len(config["system"]["atoms"]) != int(args.case[2:]):
        p.error("case name and atom count differ")
    pseudo = {name: {"path": str((input_path.parent / path).resolve()),
                     "sha256": sha(input_path.parent / path)}
              for name, path in config["pseudopotentials"].items()}
    manifest = {"case": args.case, "ranks": args.ranks, "binary": str(binary),
                "binary_sha256": sha(binary), "input": str(input_path), "input_sha256": sha(input_path),
                "pseudopotentials": pseudo,
                "environment": {k: v for k, v in env.items() if k.startswith(("BAANDI_", "MKL_", "OMP_"))},
                "commands": [], "scope": "diagnostic instrumentation overhead, no optimization A/B"}
    (root / "input.original.toml").write_bytes(input_path.read_bytes())
    wrapper = root / "rank.sh"
    wrapper.write_text('''#!/usr/bin/env bash
set -euo pipefail
rank=${OMPI_COMM_WORLD_RANK:-0}
exec /usr/bin/time -f '{"exit_code":%x,"user_s":%U,"system_s":%S,"max_rss_kib":%M,"minor_faults":%R,"major_faults":%F,"voluntary_context_switches":%w,"involuntary_context_switches":%c,"filesystem_inputs":%I,"filesystem_outputs":%O}' -o "$LAUNCHIO_RESOURCE_DIR/resource-rank-$rank.json" "$@"
''')
    launcher = ["taskset", "-c", "0"]
    if args.ranks == 8:
        launcher = ["mpirun", "-np", "8", "--cpu-set", "0-7", "--map-by", "core", "--bind-to", "core", "--report-bindings"]
    reference = None
    # Pair zero warms both modes; five retained crossed pairs follow.
    for rep in range(6):
        order = ("off", "on") if rep % 2 == 0 else ("on", "off")
        for mode in order:
            # The preceding run also raises the load average; let it decay.
            deadline = time.monotonic() + 600
            while os.getloadavg()[0] >= 0.30:
                if time.monotonic() >= deadline:
                    raise RuntimeError("load >= 0.30 for 10 minutes: reschedule the exclusive window")
                time.sleep(10)
            out = root / f"{mode}-{rep}"
            out.mkdir()
            env.update(BAANDI_PHASE_PROFILE=mode, BAANDI_STAGE_TIMERS=mode,
                       LAUNCHIO_RESOURCE_DIR=str(out))
            command = launcher + ["bash", str(wrapper), str(binary), "run", str(input_path), "--output-dir", str(out)]
            manifest["commands"].append(command)
            (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
            load_before = list(os.getloadavg())
            with (out / "stdout.log").open("w") as stdout, (out / "stderr.log").open("w") as stderr:
                start = time.monotonic()
                result = subprocess.run(command, env=env, stdout=stdout, stderr=stderr, check=False)
                wall = time.monotonic() - start
            (out / "launch-wall.json").write_text(json.dumps({"wall_s": wall, "exit_code": result.returncode,
                                                              "load_before": load_before, "load_after": list(os.getloadavg())}) + "\n")
            if result.returncode:
                raise RuntimeError(f"failed launch {out}: exit {result.returncode}")
            doc = json.loads((out / "results.json").read_text())
            if not doc["verified"] or not doc["convergence"]["converged"]:
                raise RuntimeError(f"unverified output {out}")
            science = (doc["total_energy"]["ha"].hex(), sha(out / "density.bin"), doc["convergence"]["iterations"])
            if reference is None:
                reference = science
            if science != reference:
                raise RuntimeError(f"instrumentation changed science {out}")
    subprocess.run(["python3", "-I", "scripts/lane_launchio_parse.py",
                    *[str(root / f"on-{i}" / "results.json") for i in range(1, 6)],
                    "--require-external", "--output", str(root / "tables.json")], check=True)


if __name__ == "__main__":
    main()
