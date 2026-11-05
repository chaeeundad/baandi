#!/usr/bin/env python3
"""Recheck the negative ideal-position curvature of the tested Gamma MgO model.

The stored center input is the recipe. UPFs remain external; no QE data involved.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--nc-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("mgo-curvature.json"))
    args = parser.parse_args()
    prior = json.loads(args.output.read_text())
    # Both original and refined evidence retain the center input as sample 1.
    template = prior["samples"][1]["input"]
    for symbol in ["Mg", "O"]:
        template = re.sub(rf"^{symbol}=.*$", f"{symbol}={json.dumps(str((args.nc_dir / (symbol + '.upf')).resolve()))}", template, flags=re.M)
    samples = []
    curvatures = []
    with tempfile.TemporaryDirectory() as tmp:
        folder = Path(tmp)
        for h in [0.005, 0.0025, 0.00125]:
            values = []
            for shift in [-h, 0., h]:
                f = 0.5 + shift / 7.96
                text = template.replace("fractional=[0.5,0.5,0.5]", f"fractional=[{f},{f},{f}]")
                inp = folder / "input.toml"
                inp.write_text(text)
                out = folder / f"out-{h}-{shift}"
                command = [str(args.binary.resolve()), "run", str(inp), "--output-dir", str(out)]
                run = subprocess.run(command, capture_output=True, text=True)
                if run.returncode:
                    raise RuntimeError(run.stderr)
                data = json.loads((out / "results.json").read_text())
                energy = data["total_energy"]["ha"]
                values.append(energy)
                samples.append({"cartesian_shift_each_axis_bohr": shift, "energy_ha": energy, "input": text})
            curvature = (values[0] - 2 * values[1] + values[2]) / (3 * h * h)
            curvatures.append({"step_each_axis_bohr": h, "directional_curvature_ha_bohr2": curvature})
            print(f"h={h} Bohr: curvature={curvature:.17e} Ha/Bohr^2")
            assert curvature < 0, "negative curvature of this finite model no longer reproduced"
    result = {
        "scope": "Gamma-only MgO primitive fixed-cell model; not bulk structural qualification",
        "samples": samples,
        "curvatures": curvatures,
        "command": f"python3 -I tests/reference/relaxation/check_mgo_curvature.py --binary {args.binary} --nc-dir {args.nc_dir}",
        "pp_sha256": {s: hashlib.sha256((args.nc_dir / (s + '.upf')).read_bytes()).hexdigest() for s in ["Mg", "O"]},
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
