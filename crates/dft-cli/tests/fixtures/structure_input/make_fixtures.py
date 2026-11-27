#!/usr/bin/env python3
"""Write the SI1 structure-reader fixtures and the independent ASE reading of each.

Test-only and offline: needs ASE (validated with 3.26.0, for example the
workstation venv ``rust-dft-wt/venv-convert``). The BAANDI engine never runs
this. Every structure file below is either written by ASE's own writers or
hand-written to cover layouts ASE does not emit (negative volume scale, three
scale factors, Cartesian selective dynamics, extra columns, several frames);
``reference.json`` holds what ASE's *reader* makes of each file, which
``tests/all/structure_input.rs`` compares with BAANDI's native Rust readers.

Run from this directory:
    /path/to/venv/bin/python -I make_fixtures.py
"""
import json
import math
from pathlib import Path

import ase
import ase.io
from ase import Atoms
from ase.constraints import FixAtoms, FixCartesian, FixScaled

HERE = Path(__file__).resolve().parent


def reciprocal_axis(cell, i):
    """Cartesian axis parallel to the reciprocal direction a_j x a_k, or None."""
    a, b = cell[(i + 1) % 3], cell[(i + 2) % 3]
    r = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    norm = math.sqrt(sum(x * x for x in r))
    for c in range(3):
        if all(o == c or abs(r[o]) <= 1e-12 * norm for o in range(3)):
            return c
    return None


def fixed_mask(atoms):
    """Per-atom Cartesian freeze mask from ASE's constraint objects."""
    masks = [None] * len(atoms)
    cell = atoms.cell[:].tolist()

    def merge(i, m):
        old = masks[i] or [False] * 3
        new = [bool(x or y) for x, y in zip(old, m)]
        masks[i] = new if any(new) else None

    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            for i in c.index:
                merge(int(i), [True] * 3)
        elif isinstance(c, FixCartesian):
            for i in c.index:
                merge(int(i), [bool(x) for x in c.mask])
        elif isinstance(c, FixScaled):
            for i in c.index:
                frozen = [bool(x) for x in c.mask]
                if all(frozen):
                    merge(int(i), [True] * 3)
                    continue
                m = [False] * 3
                for d in range(3):
                    if frozen[d]:
                        axis = reciprocal_axis(cell, d)
                        assert axis is not None, "fixture uses a non-axis lattice constraint"
                        m[axis] = True
                merge(int(i), m)
        else:
            raise SystemExit(f"unexpected constraint {c!r}")
    return masks


def reference(name, fmt, index=-1):
    atoms = ase.io.read(str(HERE / name), format=fmt, index=index)
    moments = None
    if atoms.has("initial_magmoms"):
        moments = [float(x) for x in atoms.get_initial_magnetic_moments()]
    return {
        "file": name,
        "format": "poscar" if fmt == "vasp" else fmt,
        "frame": index,
        "lattice_rows_angstrom": atoms.cell[:].tolist(),
        "species": atoms.get_chemical_symbols(),
        "fractional": atoms.get_scaled_positions(wrap=False).tolist(),
        "fixed": fixed_mask(atoms),
        "magmoms": moments,
        "pbc": [bool(x) for x in atoms.pbc],
    }


def write(name, text):
    (HERE / name).write_text(text)


def main():
    cases = []
    # 1-2: ASE writers, fcc primitive Si (Direct and Cartesian).
    si = Atoms(
        "Si2",
        cell=[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]],
        scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]],
        pbc=True,
    )
    ase.io.write(str(HERE / "si_direct.POSCAR"), si, format="vasp", direct=True)
    cases.append(reference("si_direct.POSCAR", "vasp"))
    ase.io.write(str(HERE / "si_cartesian.vasp"), si, format="vasp", direct=False)
    cases.append(reference("si_cartesian.vasp", "vasp"))

    # 3: orthorhombic slab, selective dynamics from FixAtoms and a z-only FixScaled.
    slab = Atoms(
        "Al2O2",
        cell=[[4.1, 0.0, 0.0], [0.0, 3.9, 0.0], [0.0, 0.0, 16.0]],
        scaled_positions=[[0.0, 0.0, 0.1], [0.5, 0.5, 0.1], [0.5, 0.0, 0.21], [0.0, 0.5, 0.2]],
        pbc=True,
    )
    slab.set_constraint([FixAtoms(indices=[0, 1]), FixScaled(3, mask=[False, False, True])])
    ase.io.write(str(HERE / "slab_selective.POSCAR"), slab, format="vasp", direct=True)
    cases.append(reference("slab_selective.POSCAR", "vasp"))

    # 4: hand-written: triclinic cell, negative scale = target volume, two species.
    write(
        "triclinic_volume.POSCAR",
        "LiF triclinic, scale = -volume\n"
        "  -71.25\n"
        "  4.02  0.10  0.05\n"
        "  1.30  3.95 -0.07\n"
        " -0.20  0.60  4.11\n"
        "  Li  F\n"
        "  2  2\n"
        "direct\n"
        "  0.0  0.0  0.0 Li1\n"
        "  0.5  0.5  0.0 Li2\n"
        "  0.5  0.0  0.5 F1\n"
        "  0.03 0.49 0.52 F2\n",
    )
    cases.append(reference("triclinic_volume.POSCAR", "vasp"))

    # 5: hand-written: three scale factors, Cartesian, selective dynamics.
    write(
        "scale3_cartesian.POSCAR",
        "three scale factors\n"
        "  1.0 1.1 0.9\n"
        "  5.0 0.0 0.0\n"
        "  0.0 5.0 0.0\n"
        "  0.0 0.0 5.0\n"
        "  H O\n"
        "  2 1\n"
        "Selective dynamics\n"
        "Cartesian\n"
        "  0.0 0.0 0.0  F F F\n"
        "  0.7 0.6 0.1  T T T\n"
        "  1.1 1.3 2.2  T F T\n",
    )
    cases.append(reference("scale3_cartesian.POSCAR", "vasp"))

    # 6: ASE extxyz writer: triclinic, FixAtoms (move_mask L:1), initial moments.
    tri = Atoms(
        "FeNiO2",
        cell=[[4.2, 0.0, 0.0], [0.4, 4.1, 0.0], [0.3, -0.2, 4.4]],
        scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.71]],
        pbc=True,
    )
    tri.set_initial_magnetic_moments([2.5, -1.25, 0.0, 0.0])
    tri.set_constraint(FixAtoms(indices=[2]))
    ase.io.write(str(HERE / "magnetic.extxyz"), tri, format="extxyz")
    cases.append(reference("magnetic.extxyz", "extxyz"))

    # 7: FixCartesian -> move_mask L:3 (ASE's writer, read back by ASE).
    fc = si.copy()
    fc.set_constraint(FixCartesian(1, mask=[True, False, True]))
    ase.io.write(str(HERE / "fixcartesian.extxyz"), fc, format="extxyz")
    cases.append(reference("fixcartesian.extxyz", "extxyz"))

    # 8: three frames (a short trajectory); default = last frame, also frame 0.
    frames = []
    for k in range(3):
        f = si.copy()
        f.positions[1] += [0.01 * k, -0.02 * k, 0.005 * k]
        frames.append(f)
    ase.io.write(str(HERE / "trajectory.extxyz"), frames, format="extxyz")
    cases.append(reference("trajectory.extxyz", "extxyz", -1))
    cases.append(reference("trajectory.extxyz", "extxyz", 0))

    # 9: hand-written: extra columns (forces, integer tags), quoted info values with spaces.
    write(
        "extra_columns.extxyz",
        "3\n"
        'config_type="a b c" energy=-12.5 Lattice="6.0 0.0 0.0 0.0 6.5 0.0 0.0 0.0 7.0" '
        "Properties=species:S:1:pos:R:3:forces:R:3:tags:I:1 pbc=\"T T T\"\n"
        "O 0.0 0.0 0.0  0.1 0.2 0.3 1\n"
        "H 0.95 0.1 0.0  -0.1 0.0 0.0 0\n"
        "H -0.2 0.9 0.05  0.0 -0.1 0.0 0\n",
    )
    cases.append(reference("extra_columns.extxyz", "extxyz"))

    (HERE / "reference.json").write_text(
        json.dumps(
            {
                "reader": "ase",
                "reader_version": ase.__version__,
                "note": "ASE's reading of each file; fixed = Cartesian freeze mask (true = frozen), FixScaled mapped to the parallel Cartesian axis",
                "cases": cases,
            },
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
