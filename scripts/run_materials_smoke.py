#!/usr/bin/env python3
"""Reproduce bulk NC smoke runs. No pseudopotential contents are copied.

DFT_NC_DIR=/path/to/ex python3 -I scripts/run_materials_smoke.py
Exit 0 requires every case to converge and preserve the valence electron count.
Input hashes identify the manifest and UPFs; results are not accuracy evidence.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time
import re

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / 'tests/reference/materials'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/dft-cli')
    parser.add_argument('--cases', default='si-20,si-30,si-40,mgo-40,mgo-50,mgo-60')
    parser.add_argument('--output', type=Path, default=MATERIALS / 'results.json')
    parser.add_argument('--work-dir', type=Path)
    parser.add_argument('--implicit-spin', action='store_true')
    parser.add_argument('--default-floor', action='store_true')
    args = parser.parse_args()
    nc = Path(os.environ['DFT_NC_DIR']).resolve()
    manifest_path = ROOT / 'docs/data/nc-datasets.json'
    manifest = json.loads(manifest_path.read_text())
    manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    datasets = {d['element']: d for d in manifest['datasets']}
    rows = []
    with tempfile.TemporaryDirectory(prefix='p13d-') as temp:
        work = args.work_dir or Path(temp)
        work.mkdir(parents=True, exist_ok=True)
        for name in args.cases.split(','):
            source = MATERIALS / (name + '.toml')
            text = source.read_text()
            species = re.findall(r'species = "([A-Za-z]+)"', text)
            cutoff_ry = float(re.search(r"ecutwfc_ry = ([0-9.]+)", text)[1])
            assert manifest_hash in text, 'manifest hash mismatch'
            for element in sorted(set(species)):
                filename = element + '.upf'
                path = nc / filename
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                assert digest == datasets[element]['sha256'] and digest in text
                text = text.replace(f'"{filename}"', json.dumps(str(path)))
            if args.implicit_spin:
                text = text.replace('[spin]\nmode = "unpolarized"\n', '')
            if args.default_floor:
                text = text.replace('density_floor = 1e-10\n', '')
            inp = work / (name + '.toml')
            inp.write_text(text)
            output = work / name
            start = time.monotonic()
            proc = subprocess.run([str(args.binary.resolve()), 'run', str(inp),
                                   '--output-dir', str(output)], capture_output=True, text=True)
            wall = time.monotonic() - start
            (work / (name + '.log')).write_text(proc.stdout + proc.stderr)
            expected = sum(datasets[e]['z_valence'] for e in species)
            row = dict(case=name, cutoff_ha=cutoff_ry/2,
                       exit_code=proc.returncode, wall_s=wall, expected_electrons=expected,
                       input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
            result_path = output / 'results.json'
            if result_path.exists() and proc.returncode in (0, 3):
                result = json.loads(result_path.read_text())
                density = json.loads((output / 'density.json').read_text())
                row.update(converged=result['convergence']['converged'],
                           iterations=result['convergence']['iterations'],
                           total_energy_ha=result['total_energy']['ha'],
                           band_gap=result.get('band_gap'),
                           electron_count=density['normalization']['integral_electrons'],
                           density_sha256=density['sha256'],
                           xc_floor=result['energies'].get('xc_floor'),
                           clipped_negative_charge_per_iteration=[[h['energies']['xc_floor'][k] for k in ['input_negative_charge_electrons','output_negative_charge_electrons']] for h in result['convergence']['history']])
                cert = result['convergence']['last_iteration'].get('fixed_point_certificate')
                row['certificate'] = None if cert is None else {k:cert[k] for k in
                    ['density_residual', 'occupied_eigenvalue_change_ha']}
                if cert is not None:
                    row['certificate']['max_strict_residual_ha'] = max(x for k in cert['strict_residuals_per_k_ha'] for x in k)
                row['passed'] = (proc.returncode == 0 and row['converged'] and cert is not None
                                 and abs(row['electron_count']-expected) < 1e-8)
            else:
                row.update(passed=False, error=proc.stderr[-4000:])
            rows.append(row)
            args.output.write_text(json.dumps({'manifest_sha256':manifest_hash,
                                              'cases': rows}, indent=2) + '\n')
            print(name, row.get('iterations'), row.get('total_energy_ha'), row['exit_code'], flush=True)
    return 0 if all(r['passed'] for r in rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
