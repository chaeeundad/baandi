#!/usr/bin/env python3
"""Mac correctness smoke only: no wall-time or performance analysis."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import subprocess

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('lane_parser', ROOT / 'scripts/lane_nonlocal_parse.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scientific(value):
    if isinstance(value, float):
        return {'f64_bits': struct.pack('>d', value).hex()}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items()
                if k not in {'timings', 'wall_s'}}
    return value


def signature(result):
    keys = ('energies', 'total_energy', 'eigenvalues', 'occupations',
            'band_accuracy', 'work', 'verified', 'verification')
    out = {key: result[key] for key in keys}
    out['history'] = result['convergence']['history']
    out['iterations'] = result['convergence']['iterations']
    out['converged'] = result['convergence']['converged']
    return scientific(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--upf', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    rows = []
    for case in ('si8', 'si32'):
        source = args.out / (case + '.toml')
        subprocess.run(['python3', '-I', 'scripts/fast_screen_dav.py',
                        '--upf', str(args.upf.resolve()), '--case', case,
                        '--seed', '601', '--amp', '0.001', '--threads', '1',
                        '--emit-input', str(source)], check=True, cwd=ROOT)
        reference = reference_density = counted = None
        for mode, profile, clocks in [('off', 'off', 'off'),
                                      ('counts', 'on', 'off'), ('timed', 'on', 'on')]:
            out = args.out / f'{case}-{mode}'
            env = dict(os.environ, BAANDI_PHASE_PROFILE=profile, BAANDI_STAGE_TIMERS=clocks)
            with (args.out / f'{case}-{mode}.log').open('w') as log:
                subprocess.run([str(args.binary.resolve()), 'run', str(source),
                                '--output-dir', str(out)], env=env, check=True,
                               stdout=log, stderr=subprocess.STDOUT, cwd=ROOT)
            result = json.loads((out / 'results.json').read_text())
            current = signature(result)
            density = digest(out / 'density.bin')
            assert result['convergence']['converged'] and result['verified']
            if mode == 'off':
                reference, reference_density = current, density
                assert 'phase_profile' not in result
            else:
                assert current == reference, (case, mode, 'scientific f64 bits/history differ')
                assert density == reference_density, (case, mode, 'density bytes differ')
                table = parser.parse_result(out / 'results.json', 'mac')
                assert len(table) == 48
                assert not any('seconds' in key for row in table for key in row)
                (out / 'nonlocal-table.json').write_text(json.dumps(table, indent=2)+'\n')
                counters = [{k: v for k, v in b['max'].items()
                             if k.startswith('nonlocal.') and k.split('.')[-1] not in parser.STAGES}
                            for b in result['phase_profile']['buckets']]
                if mode == 'counts':
                    counted = counters
                    for bucket in result['phase_profile']['buckets']:
                        for stage in parser.STAGES:
                            assert bucket['max']['nonlocal.'+stage] == 0
                else:
                    assert counters == counted, (case, 'count-only/timed work differs')
                    fine = next(b for b in result['phase_profile']['buckets'] if b['bucket']=='fine')
                    assert fine['max']['nonlocal.projection'] > 0
                    assert fine['max']['nonlocal.back'] > 0
            rows.append({'case': case, 'mode': mode, 'exit_code': 0,
                         'energy_ha': result['total_energy']['ha'], 'density_sha256': density,
                         'input_sha256': result['input_sha256'],
                         'iterations': result['convergence']['iterations'], 'verified': result['verified'],
                         'scientific_sha256': hashlib.sha256(json.dumps(current, sort_keys=True).encode()).hexdigest(),
                         'results_path': str(out / 'results.json')})
    report = {'purpose': 'correctness only; no Mac performance conclusion',
              'binary_sha256': digest(args.binary), 'upf_sha256': digest(args.upf),
              'max_energy_difference_ha': 0, 'runs': rows,
              'checks': ['scientific f64 bits', 'density bytes', 'history and work',
                         'verification attempts', 'count-only clocks zero',
                         'timed counters equal count-only', 'parser on actual CLI files']}
    (args.out / 'equivalence.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS: six CLI runs, four off/on scientific and density identities; four real-file parser checks')


if __name__ == '__main__':
    main()
