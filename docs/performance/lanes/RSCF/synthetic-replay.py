#!/usr/bin/env python3
"""Generate an analytic (NOT QE) three-step file-format oracle and run the CLI.
All output goes to --scratch. No external pseudopotential or network is used.
"""
import argparse
import json
import math
from pathlib import Path
import struct
import subprocess


def write(path, miller, rho):
    path.mkdir(parents=True, exist_ok=True)
    def record(data):
        marker = struct.pack('<i', len(data))
        return marker + data + marker
    content = record(struct.pack('<3i', 1, len(miller), 1))
    content += record(struct.pack('<9d', 1, 0, 0, 0, 1, 0, 0, 0, 1))
    content += record(b''.join(struct.pack('<3i', *m) for m in miller))
    content += record(b''.join(struct.pack('<2d', r, 0) for r in rho))
    (path / 'charge-density.dat').write_bytes(content)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--scratch', type=Path, required=True)
    args = p.parse_args()
    root = args.scratch.resolve()
    root.mkdir(parents=True, exist_ok=True)
    t = math.tau
    text = f'''[system]
name = "synthetic analytic replay, NOT QE"
lattice_unit = "bohr"
lattice_vectors = [[{t},0,0],[0,{t},0],[0,0,{t}]]
[[system.atoms]]
species = "Xx"
fractional = [0,0,0]
[pseudopotentials]
Xx = "unused-synthetic.UPF"
[basis]
ecutwfc_ry = 0.5
kpoints = [[0,0,0]]
[xc]
functional = "lda-pz"
[scf]
driver = "reference"
mixing_beta = 0.5
'''
    (root / 'input.toml').write_text(text)
    miller = [(i,j,k) for i in range(-1,2) for j in range(-1,2) for k in range(-1,2)
              if i*i+j*j+k*k <= 2 and (i > 0 or (i == 0 and (j > 0 or (j == 0 and k >= 0))))]
    points = [(0,0), (1,1.5), (1,93/97), (1,1)]
    volume = t**3
    def coefficients(x):
        return [2/volume if m == (0,0,0) else x[0]*1e-3 if m == (1,0,0)
                else x[1]*1e-3 if m == (0,1,0) else 0.0 for m in miller]
    for k in range(1,4):
        x = points[k-1]
        f = (2-2*x[0], 3-3*x[1])
        y = coefficients((x[0]+f[0], x[1]+f[1]))
        mix = coefficients(points[k])
        write(root / f'k{k}-out', miller, y)
        write(root / f'k{k}-mix', list(reversed(miller)), list(reversed(mix)))
        dr2 = 8*math.pi*volume*1e-6*sum(v*v for v in f)
        (root / f'k{k}-mix/stdout').write_text(f'estimated scf accuracy = {dr2:.16e} Ry\n')
    write(root / 'k1-in', miller, coefficients(points[0]))
    write(root / 'k1-beta1', miller, coefficients((2,3)))
    command = [str(args.binary.resolve()), 'refscf-check', str(root / 'input.toml'),
               '--qe', str(root), '--steps', '3', '--atol', '1e-14', '--rtol', '1e-12']
    good = subprocess.run(command, capture_output=True, text=True)
    (root / 'success.jsonl').write_text(good.stdout)
    (root / 'success.stderr').write_text(good.stderr)
    original = (root / 'k2-mix/stdout').read_text()
    (root / 'k2-mix/stdout').write_text('estimated scf accuracy = 99.000000000 Ry\n')
    bad = subprocess.run(command, capture_output=True, text=True)
    (root / 'failure.jsonl').write_text(bad.stdout)
    (root / 'failure.stderr').write_text(bad.stderr)
    (root / 'k2-mix/stdout').write_text(original)
    rows = [json.loads(line) for line in good.stdout.splitlines() if line.startswith('{')]
    evidence = {'source': 'synthetic analytic fixed point; NOT QE', 'command': command,
                'success_exit': good.returncode, 'injected_dr2_failure_exit': bad.returncode,
                'first_failing_term': bad.stderr.strip(), 'records': rows}
    (root / 'evidence.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(json.dumps({k:v for k,v in evidence.items() if k != 'records'}, indent=2))
    if good.returncode != 0 or bad.returncode == 0 or 'k=2 dr2' not in bad.stderr:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
