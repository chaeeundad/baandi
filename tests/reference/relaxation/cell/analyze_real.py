"""Summarize actual P23A material trajectories; never substitutes for a run."""
import argparse
import json
from pathlib import Path
import shutil


def summarize(run_dir, output):
    records = json.loads((run_dir / "summary.json").read_text())
    for row in records:
        result_path = run_dir / row["case"] / "results.json"
        if not result_path.exists():
            continue
        relaxation = json.loads(result_path.read_text())["relaxation"]
        frames = relaxation.get("frames", [])
        if not frames or "lattice_bohr" not in row:
            continue
        a = row["lattice_bohr"]
        p = row["frozen_final"]["pressure_ha_bohr3"]
        slope = None
        for frame in reversed(frames[:-1]):
            old_a = 2 * frame["cell_columns_bohr"][0][1]
            if abs(a - old_a) > 1e-5:
                tensor = frame["stress_ha_bohr3"]
                old_p = -sum(tensor[i][i] for i in range(3)) / 3
                slope = (p - old_p) / (a - old_a)
                break
        row["local_dp_da_ha_bohr4"] = slope
        tolerance = relaxation["settings"]["stress_tol_ha_bohr3"]
        # For isotropic D=d*I, the driver stops on |pressure|/d.
        physical_pressure_tolerance = tolerance * a / row["start_lattice_bohr"]
        row["physical_pressure_tolerance_ha_bohr3"] = physical_pressure_tolerance
        row["estimated_lattice_tolerance_bohr"] = (
            physical_pressure_tolerance / abs(slope) if slope else None
        )
    start_pairs = []
    cutoff_pairs = []
    for material, cutoffs in [("Si", (30, 40)), ("MgO", (50, 60))]:
        selected = [r for r in records if r["material"] == material and "lattice_bohr" in r]
        for cutoff in cutoffs:
            pair = sorted([r for r in selected if r["cutoff_ha"] == cutoff], key=lambda r: r["start_factor"])
            if len(pair) != 2:
                continue
            delta = pair[1]["lattice_bohr"] - pair[0]["lattice_bohr"]
            bounds = [r.get("estimated_lattice_tolerance_bohr") for r in pair]
            bound = sum(bounds) if all(x is not None for x in bounds) else None
            start_pairs.append(dict(material=material, cutoff_ha=cutoff,
                                    expanded_minus_compressed_bohr=delta,
                                    estimated_combined_optimizer_tolerance_bohr=bound,
                                    within_estimated_optimizer_tolerance=abs(delta) <= bound if bound else None))
        for factor in (0.98, 1.02):
            pair = sorted([r for r in selected if r["start_factor"] == factor], key=lambda r: r["cutoff_ha"])
            if len(pair) == 2:
                cutoff_pairs.append(dict(material=material, start_factor=factor,
                                         higher_minus_lower_cutoff_bohr=pair[1]["lattice_bohr"]-pair[0]["lattice_bohr"]))
    report = dict(cases=records, start_comparisons=start_pairs, cutoff_comparisons=cutoff_pairs,
                  tolerance_method="First-order conversion of the unchanged stress stopping tolerance using each trajectory's last resolved dP/da; an estimated lattice tolerance, not an independent bound or a changed acceptance threshold.",
                  scope="Finite frozen-Miller minima and rebuilt nominal-cutoff diagnostics; no experiment or external-code comparison.")
    output.mkdir(parents=True, exist_ok=True)
    (output / "real-validation.json").write_text(json.dumps(report, indent=2)+"\n")
    lines = ["# P23A actual material runs", "", "Pressure units: Ha/Bohr^3; lattice: Bohr; time: seconds.", "",
             "| Case | Lattice | Frozen P | Rebuilt P | Rebuilt - frozen P | Time | Exit |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in records:
        if "lattice_bohr" not in row:
            lines.append(f"| {row['case']} | unavailable | | | | {row['elapsed_s']:.3f} | {row['exit_code']} |")
        else:
            frozen = row["frozen_final"]["pressure_ha_bohr3"]
            rebuilt = row["rebuilt_final"]["pressure_ha_bohr3"]
            lines.append(f"| {row['case']} | {row['lattice_bohr']:.12f} | {frozen:.9e} | {rebuilt:.9e} | {rebuilt-frozen:.9e} | {row['elapsed_s']:.3f} | {row['exit_code']} |")
        for suffix in ("toml", "stdout", "stderr"):
            source = run_dir / (row["case"]+"."+suffix)
            if source.exists():
                shutil.copyfile(source, output/source.name)
    lines += ["", "## Start differences", "", "Optimizer lattice tolerances below are local slope estimates from the original stress tolerance; they are not an independently specified lattice tolerance.", ""]
    for row in start_pairs:
        lines.append(f"- {row['material']} {row['cutoff_ha']} Ha: expanded minus compressed = {row['expanded_minus_compressed_bohr']:.9e} Bohr; combined estimated tolerance = {row['estimated_combined_optimizer_tolerance_bohr']:.9e} Bohr; within estimate = {row['within_estimated_optimizer_tolerance']}.")
    lines += ["", "## Cutoff dependence", ""]
    for row in cutoff_pairs:
        lines.append(f"- {row['material']} start factor {row['start_factor']}: higher minus lower cutoff lattice = {row['higher_minus_lower_cutoff_bohr']:.9e} Bohr.")
    lines += ["", report["scope"], ""]
    (output/"real-validation.md").write_text("\n".join(lines))
    print("summarized", len(records), "cases")
    # Frozen starts use different Miller sets. Preserve the diagnostic numbers,
    # but assert start agreement only in rebuild_protocol.py.


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summarize(args.run_dir, args.output)
