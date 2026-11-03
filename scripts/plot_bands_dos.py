"""Render band CSV and optional Gaussian DOS CSV using only the standard library."""
import argparse
import csv
import html
import json
import math
from pathlib import Path


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def render(bands_path, dos_path=None):
    path = Path(bands_path)
    rows = read_csv(path)
    if not rows:
        raise ValueError("bands CSV is empty")
    names = sorted((k for k in rows[0] if k.startswith("band_")),
                   key=lambda k: int(k[5:]))
    if not names:
        raise ValueError("bands CSV has no band columns")
    xs = [float(r["distance_inv_bohr"]) for r in rows]
    bands = [[float(r[k]) for r in rows] for k in names]
    dos = read_csv(dos_path) if dos_path else []
    metadata = path.with_suffix(".json")
    meta = json.loads(metadata.read_text()) if metadata.exists() else {}
    breaks = meta.get("path_break_indices", [i for i in range(1, len(xs)) if xs[i] == xs[i-1]])
    if any(not isinstance(i, int) or i <= 0 or i >= len(rows) for i in breaks) or breaks != sorted(set(breaks)):
        raise ValueError("invalid path break indices")
    if dos_path:
        dos_metadata = Path(dos_path).with_suffix(".json")
        dos_meta = json.loads(dos_metadata.read_text()) if dos_metadata.exists() else {}
        if meta and dos_meta:
            if meta.get("input_sha256") != dos_meta.get("input_sha256"):
                raise ValueError("bands and DOS must come from the same input/density setup")
            for key in ["scf_density_sha256", "pseudopotentials", "lattice_columns_bohr", "functional"]:
                if meta.get(key) != dos_meta.get(key):
                    raise ValueError("bands and DOS have incompatible " + key)
            # Both tables use their sampled VBM; align the DOS to the band's
            # absolute reference. Conversion constant is recorded by the CLI.
            shift = (dos_meta["energy_reference_hartree"] - meta["energy_reference_hartree"]) * meta["hartree_to_ev"]
            for row in dos:
                row["energy_eV"] = float(row["energy_eV"]) + shift
    energies = [v for band in bands for v in band] + [float(r["energy_eV"]) for r in dos]
    values = xs + energies + [float(r["dos_states_per_eV_per_cell"]) for r in dos]
    if not all(math.isfinite(v) for v in values):
        raise ValueError("CSV contains nonfinite numbers")
    if any(b < a for a, b in zip(xs, xs[1:])) or xs[-1] <= xs[0]:
        raise ValueError("band distances must increase cumulatively")
    low, high = min(energies + [0.0]), max(energies + [0.0])
    pad = max((high - low) * .05, .1)
    low, high = low - pad, high + pad
    x = lambda v: 80 + (v - xs[0]) / (xs[-1] - xs[0]) * 640
    y = lambda v: 520 - (v - low) / (high - low) * 460
    elements = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="600" viewBox="0 0 1000 600">',
                '<rect width="1000" height="600" fill="white"/>']
    def text(px, py, value):
        elements.append(f'<text x="{px}" y="{py}" font-family="sans-serif" font-size="14" text-anchor="middle">{html.escape(str(value))}</text>')
    def line(x1, y1, x2, y2, extra=""):
        elements.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#888" {extra}/>')
    def poly(points, color):
        data = " ".join(f"{a:.6f},{b:.6f}" for a, b in points)
        elements.append(f'<polyline points="{data}" fill="none" stroke="{color}" stroke-width="1.5"/>')
    title = meta.get("title", "Band structure")
    text(500, 30, title)
    for i in range(7):
        e = low + (high - low) * i / 6
        text(42, y(e) + 5, f"{e:.2f}")
    line(80, 60, 80, 520)
    line(80, 520, 720, 520)
    line(80, y(0), 950 if dos else 720, y(0), 'stroke-dasharray="6 4"')
    labels_path = path.with_name("bands_labels.csv")
    if labels_path.exists():
        grouped = {}
        for label in read_csv(labels_path):
            grouped.setdefault(float(label["distance"]), []).append(label["label"])
        for distance, labels in grouped.items():
            px = x(distance)
            line(px, 60, px, 520)
            text(px, 545, "|".join(labels))
    boundaries = [0] + breaks + [len(rows)]
    for band in bands:
        for start, end in zip(boundaries, boundaries[1:]):
            poly(zip(map(x, xs[start:end]), map(y, band[start:end])), "#2458a6")
    text(400, 580, "Reciprocal path (1/Bohr)")
    text(160, 48, "Energy relative to VBM (eV)")
    if dos:
        densities = [float(r["dos_states_per_eV_per_cell"]) for r in dos]
        if min(densities) < 0 or max(densities) <= 0:
            raise ValueError("DOS must be nonnegative with positive maximum")
        poly(((780 + 170 * d / max(densities), y(float(r["energy_eV"])))
              for d, r in zip(densities, dos)), "#b44623")
        line(780, 60, 780, 520)
        text(865, 555, "DOS (states/eV/cell)")
        text(865, 580, f"0 to {max(densities):.3g}")
    return "\n".join(elements + ["</svg>"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bands")
    parser.add_argument("dos", nargs="?")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    inputs = [Path(args.bands), Path(args.bands).with_suffix(".json"), Path(args.bands).with_name("bands_labels.csv")]
    if args.dos:
        inputs.extend([Path(args.dos), Path(args.dos).with_suffix(".json")])
    if any(output.resolve() == p.resolve() or (output.exists() and p.exists() and output.samefile(p)) for p in inputs):
        parser.error("output must not overwrite an input")
    output.write_text(render(args.bands, args.dos), encoding="utf-8")


if __name__ == "__main__":
    main()
