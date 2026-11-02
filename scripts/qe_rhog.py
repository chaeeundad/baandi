"""Strict QE records and independent, standard-library density Fourier analysis."""
import cmath
import hashlib
import itertools
import json
import math
from pathlib import Path
import struct


class FormatError(ValueError):
    """Inconsistent scientific input or unsupported format."""


def require(condition, message):
    if not condition:
        raise FormatError(message)


def finite(value):
    require(type(value) in (int, float) and math.isfinite(value), 'expected finite number')
    return value


def miller(value):
    require(isinstance(value, (list, tuple)) and len(value) == 3
            and all(type(x) is int for x in value), 'invalid Miller triplet')
    return tuple(value)


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (ValueError, UnicodeError) as exc:
        raise FormatError(f'invalid JSON: {exc}') from exc


def read_qe(path):
    """Return (Miller -> complex coefficient, metadata); never infer a format."""
    with Path(path).open('rb') as stream:
        def record(size, number):
            marker = stream.read(4)
            require(len(marker) == 4, f'record {number}: missing leading marker')
            length, = struct.unpack('<i', marker)
            require(length == size, f'record {number}: wrong size {length}, expected {size}')
            # Check remaining size before allocating from an untrusted header.
            position = stream.tell()
            stream.seek(0, 2)
            remaining = stream.tell() - position
            stream.seek(position)
            require(remaining >= size + 4, f'record {number}: truncated payload or marker')
            payload = stream.read(size)
            require(stream.read(4) == marker, f'record {number}: marker mismatch')
            return payload
        gamma, count, spin = struct.unpack('<iii', record(12, 1))
        require(gamma == 0, 'unsupported gamma_only: expected 0')
        require(spin == 1, 'unsupported nspin: expected 1')
        require(count > 0, 'ngm_g must be positive')
        vectors = struct.unpack('<9d', record(72, 2))
        for x in vectors:
            finite(x)
        indices = list(struct.iter_unpack('<iii', record(12 * count, 3)))
        require(len(set(indices)) == count, 'duplicate Miller triplet')
        values = list(struct.iter_unpack('<dd', record(16 * count * spin, 4)))
        require(not stream.read(1), 'trailing bytes after record 4')
    coefficients = {m: complex(finite(re), finite(im)) for m, (re, im) in zip(indices, values)}
    return coefficients, dict(gamma_only=gamma, ngm_g=count, nspin=spin,
                              b_vectors_as_stored=vectors)


def write_qe(path, coefficients, vectors=(1., 0., 0., 0., 1., 0., 0., 0., 1.)):
    """Write a synthetic fixture, exclusively to a fresh file."""
    records = [struct.pack('<iii', 0, len(coefficients), 1), struct.pack('<9d', *vectors),
               b''.join(struct.pack('<iii', *m) for m in coefficients),
               b''.join(struct.pack('<dd', z.real, z.imag) for z in coefficients.values())]
    with Path(path).open('xb') as stream:
        for payload in records:
            marker = struct.pack('<i', len(payload))
            stream.write(marker + payload + marker)


GRID_RULE = 'Each axis: -(n//2) <= m <= (n-1)//2; negative even Nyquist; unique m mod n.'


def validate_grid(grid, indices):
    require(isinstance(grid, (list, tuple)) and len(grid) == 3
            and all(type(n) is int and n > 0 for n in grid), 'invalid fft_grid')
    seen = set()
    for item in indices:
        m = miller(item)
        require(all(-(n // 2) <= x <= (n - 1) // 2 for x, n in zip(m, grid)),
                f'Miller {m} unrepresentable without aliasing: {GRID_RULE}')
        wrapped = tuple(x % n for x, n in zip(m, grid))
        require(wrapped not in seen, f'duplicate/aliased Miller {m}')
        seen.add(wrapped)


def dft(density, grid, requested):
    """Separable 3D DFT, exp(-2*pi*i*m.r), normalized by total grid size."""
    requested = list(requested)
    validate_grid(grid, requested)
    require(len(density) == math.prod(grid), 'density size differs from grid')
    n1, n2, n3 = grid
    data = [complex(x) for x in density]
    for axis, n in enumerate(grid):
        stride = math.prod(grid[axis + 1:])
        block = stride * n
        phases = [[cmath.exp(-2j * math.pi * k * j / n) / n for j in range(n)]
                  for k in range(n)]
        out = [0j] * len(data)
        for start in range(0, len(data), block):
            for offset in range(stride):
                lane = [data[start + offset + j * stride] for j in range(n)]
                for k in range(n):
                    out[start + offset + k * stride] = sum(z * w for z, w in zip(lane, phases[k]))
        data = out
    return {m: data[((m[0] % n1) * n2 + m[1] % n2) * n3 + m[2] % n3]
            for m in requested}


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def dot(a, b):
    return math.fsum(x*y for x, y in zip(a, b))


def load_ours(directory):
    """Require density metadata plus results.json for independent density-sphere membership."""
    directory = Path(directory)
    meta = read_json(directory / 'density.json')
    require(meta['schema_version'] == 1, 'unsupported density schema_version')
    for key, expected in [('dtype', 'f64'), ('endianness', 'little'),
                          ('ordering', 'C: (i*n2+j)*n3+k; k fastest'),
                          ('units', 'electrons per Bohr^3')]:
        require(meta[key] == expected, f'unsupported density {key}')
    grid = meta['fft_grid']
    validate_grid(grid, [])
    raw = (directory / 'density.bin').read_bytes()
    require(hashlib.sha256(raw).hexdigest() == meta['sha256'], 'density.bin SHA-256 mismatch')
    require(len(raw) == 8 * math.prod(grid), 'density.bin wrong size')
    density = [finite(x[0]) for x in struct.iter_unpack('<d', raw)]
    volume = finite(meta['volume_bohr3'])
    require(volume > 0, 'volume must be positive')
    integral = math.fsum(density) * volume / len(density)
    expected = finite(meta['normalization']['expected_electrons'])
    require(abs(integral - expected) <= 1e-10, 'density integral differs from expected electrons')
    require(abs(integral - finite(meta['normalization']['integral_electrons'])) <= 1e-10,
            'recorded integral differs from density')
    lattice = meta['lattice_columns_bohr']
    require(len(lattice) == 3 and all(len(row) == 3 for row in lattice), 'invalid lattice')
    for row in lattice:
        for x in row:
            finite(x)
    columns = list(zip(*lattice))
    determinant = dot(columns[0], cross(columns[1], columns[2]))
    require(abs(abs(determinant) - volume) <= 1e-10 * volume, 'lattice/volume mismatch')
    reciprocal = [[2 * math.pi * x / determinant for x in cross(columns[j], columns[k])]
                  for j, k in [(1, 2), (2, 0), (0, 1)]]
    results = read_json(directory / 'results.json')
    require(results['basis']['fft_grid'] == grid, 'results/density grid mismatch')
    require(results['cell']['lattice_columns_bohr'] == lattice, 'results/density lattice mismatch')
    cutoff = finite(results['cutoff']['ecutwfc_ha'])
    require(cutoff > 0, 'cutoff must be positive')
    limit = 8 * cutoff * (1 + 1e-12)
    bounds = [math.ceil(math.sqrt(dot(a, a) * limit) / (2 * math.pi)) for a in columns]
    indices = []
    for m in itertools.product(*(range(-b, b + 1) for b in bounds)):
        g = [sum(m[j] * reciprocal[j][i] for j in range(3)) for i in range(3)]
        if dot(g, g) <= limit:
            indices.append(m)
    require(len(indices) == results['basis']['n_density_g'], 'density sphere count mismatch')
    validate_grid(grid, indices)
    return dft(density, grid, indices), meta, dict(integral_electrons=integral,
        expected_electrons=expected, charge_error_electrons=abs(integral-expected),
        sha256=meta['sha256'], grid_rule=GRID_RULE, density_g_count=len(indices))
