#!/usr/bin/env python3
"""Check evidence consistency, source identity and optional independent replay."""
import argparse
import copy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(run):
    return run['model'], run['policy'], run['arm']


def normalized(run):
    value = copy.deepcopy(run)
    for field in ('wall_seconds', 'elapsed_seconds'):
        value.pop(field, None)
    for h in value.get('history', []):
        if h['certificate_debug'] is not None:
            # FixedPointCertificate ends with PhaseTimings; retain every numerical field.
            assert ', timings: PhaseTimings {' in h['certificate_debug']
            h['certificate_debug'] = h['certificate_debug'].split(', timings: PhaseTimings {')[0]
    return value


def validate(directory):
    data = json.loads((directory / 'results.json').read_text())
    for path, digest in data['sources'].items():
        assert sha(ROOT / path) == digest, f'source changed: {path}'
    runs = {key(r): r for r in data['runs']}
    assert len(runs) == len(data['runs'])
    for k, r in runs.items():
        stem = '--'.join(k)
        assert json.loads((directory / (stem + '.json')).read_text()) == r
        inp = directory / (stem + '.input.json')
        assert sha(inp) == r['input_sha256']
        assert r['process_exit_code'] == 0 and not r['stderr'], k
        assert r['iterations'] == len(r['history'])
        assert r['converged'] == (r['state'] == 'Converged')
        assert (r['certified_energy_hartree'] is not None) == r['converged']
        assert r['charge_sloshing_hint'] == 'not_computed'
        if r['converged']:
            assert r['history'][-1]['certificate_debug'] is not None
            assert r['final_residual'] <= 1e-7
        for f in r['findings']:
            window = [h for h in r['history'] if f['first'] <= h['iteration'] <= f['last']]
            assert window and all(h['revision'] == f['revision'] for h in window)
            assert f['evidence_debug'] and f['version'] == '1.0'
            if f['rule'] == 'Stagnation':
                residuals = [h['residual'] for h in window]
                assert len(residuals) == 6 and min(residuals) > 1e-7
                assert max(residuals) / min(residuals) <= 1.2
                assert residuals[-1] / residuals[0] >= 0.95
            if f['rule'] == 'Oscillation':
                assert len(window) == 7
                deltas = [b['energy_hartree'] - a['energy_hartree'] for a, b in zip(window, window[1:])]
                assert all(abs(d) > 1e-8 for d in deltas)
                assert all(a*b < 0 for a, b in zip(deltas, deltas[1:]))
                assert window[-1]['residual'] / window[1]['residual'] >= 0.9
        if r['arm'] == 'intervention':
            baseline = runs[(r['model'], r['policy'], 'baseline')]
            reference = runs[(r['model'], 'reference', 'reference')]
            assert [a['outcome'] for a in r['audit']] == ['Accepted', 'Applied']
            a = r['audit'][-1]
            boundary = a['iteration']
            assert r['history'][:boundary] == baseline['history'][:boundary]
            assert a['history_reset'] and a['after_revision'] == a['before_revision'] + 1
            assert any(f['last'] == boundary and f['rule'] in ('Stagnation', 'Oscillation') for f in baseline['findings'])
            e, ref = r['certified_energy_hartree'], reference['certified_energy_hartree']
            delta = abs(e-ref) if e is not None and ref is not None else None
            assert r['energy_error_hartree'] == delta
            assert r['recovery_validated'] == (delta is not None and delta <= 1e-10)
    return data, runs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory', type=Path)
    p.add_argument('--repeat', type=Path)
    p.add_argument('--record', type=Path)
    a = p.parse_args()
    data, runs = validate(a.directory)
    report = {'runs_checked': len(runs), 'sources_checked': len(data['sources']),
              'validated_interventions': sum(r.get('recovery_validated', False) for r in runs.values()),
              'command': __import__('sys').argv, 'repeat_equal': None}
    if a.repeat:
        repeat, other = validate(a.repeat)
        assert data['binary_sha256'] == repeat['binary_sha256']
        assert runs.keys() == other.keys()
        for k in runs:
            assert normalized(runs[k]) == normalized(other[k]), f'replay mismatch: {k}'
        report['repeat_equal'] = True
        report['normalized_sha256'] = {
            '--'.join(k): hashlib.sha256(json.dumps(normalized(r), sort_keys=True).encode()).hexdigest()
            for k, r in runs.items()}
    if a.record:
        a.record.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'normalized_sha256'}))


if __name__ == '__main__':
    main()
