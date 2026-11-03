#!/usr/bin/env python3
"""Independent standard-library validator for scalar Bohr Gaussian cube files.

Only the single scalar, positive-axis-count variant emitted by rust-dft is
accepted. Negative atom counts (orbital-ID blocks), Angstrom axes and multi-value
voxels are deliberately rejected. Values are z-fastest. Integral uses |det voxel|.
"""
import argparse
import json
import math
from pathlib import Path


def determinant(rows):
    a, b, c = rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def read_cube(path):
    with open(path, encoding="utf-8") as stream:
        comments = [stream.readline().rstrip("\n") for _ in range(2)]
        def row(length):
            fields = stream.readline().split()
            if len(fields) != length:
                raise ValueError(f"expected {length} header fields")
            return fields
        header = row(4)
        count = int(header[0])
        if count < 0:
            raise ValueError("orbital cube unsupported")
        origin = list(map(float, header[1:]))
        dims, voxels = [], []
        for _ in range(3):
            fields = row(4)
            dims.append(int(fields[0]))
            voxels.append(list(map(float, fields[1:])))
        if min(dims) <= 0:
            raise ValueError("positive Bohr axis counts required")
        atoms = []
        for _ in range(count):
            fields = row(5)
            z = int(fields[0])
            if not 1 <= z <= 118:
                raise ValueError("invalid atomic number")
            atoms.append(dict(atomic_number=z, charge=float(fields[1]),
                              position_bohr=list(map(float, fields[2:]))))
        values = [float(x) for line in stream for x in line.split()]
    if len(values) != math.prod(dims):
        raise ValueError("scalar sample count mismatch")
    numbers = origin + [x for row in voxels for x in row] + values
    for atom in atoms:
        numbers += [atom['charge']] + atom['position_bohr']
    if not all(map(math.isfinite, numbers)):
        raise ValueError("nonfinite number")
    voxel_volume = abs(determinant(voxels))
    if not math.isfinite(voxel_volume) or voxel_volume == 0:
        raise ValueError("degenerate or unrepresentable voxel")
    integral = math.fsum(values)*voxel_volume
    if not math.isfinite(integral):
        raise ValueError("unrepresentable integral")
    return dict(comments=comments, dims=dims, origin_bohr=origin,
                lattice_vectors_bohr=[[x*n for x in row] for n, row in zip(dims, voxels)],
                atoms=atoms, values=values, integral=integral,
                minimum=min(values), maximum=max(values))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cube', type=Path)
    parser.add_argument('--values', action='store_true', help='include every sample')
    args = parser.parse_args()
    try:
        result = read_cube(args.cube)
        if not args.values:
            result.pop('values')
        print(json.dumps(result, allow_nan=False))
    except (ValueError, OSError, OverflowError) as error:
        parser.exit(1, f'{error}\n')
