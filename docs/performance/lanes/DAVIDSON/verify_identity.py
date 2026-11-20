#!/usr/bin/env python3
"""Compare scientific outputs of completed CLI runs, excluding only clocks."""
import hashlib
import json
from pathlib import Path
import sys
import struct


def strip_clocks(x):
    if isinstance(x, dict):
        return {k: strip_clocks(v) for k, v in x.items() if k not in ('timings', 'wall_s')}
    if isinstance(x, list):
        return [strip_clocks(v) for v in x]
    if isinstance(x, float):
        return ("float64", struct.pack(">d", x).hex())
    return x


def read(directory):
    directory = Path(directory)
    r = json.loads((directory / 'results.json').read_text())
    if not r.get('run', {}).get('completed'):
        raise ValueError(f'{directory}: incomplete report')
    science = {k: r[k] for k in ('total_energy', 'eigenvalues', 'occupations', 'verified', 'verification')}
    science['history'] = r['convergence']['history']
    science['converged'] = r['convergence']['converged']
    density = (directory / 'density.bin').read_bytes()
    return strip_clocks(science), density, r


def main():
    if len(sys.argv) < 3:
        raise ValueError('usage: verify_identity.py baseline_dir candidate_dir [candidate_dir ...]')
    baseline, density, report = read(sys.argv[1])
    for path in sys.argv[2:]:
        candidate, new_density, _ = read(path)
        if baseline != candidate:
            raise ValueError(f'{path}: scientific record mismatch')
        if density != new_density:
            raise ValueError(f'{path}: density bytes mismatch')
    print(json.dumps(dict(status='PASS', directories=sys.argv[1:],
                          energy_ha=report['total_energy']['ha'],
                          absolute_energy_difference_ha=0,
                          density_sha256=hashlib.sha256(density).hexdigest(),
                          scientific_history_and_verification_bit_identical=True), indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        raise SystemExit(str(error))
