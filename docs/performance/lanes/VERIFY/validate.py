#!/usr/bin/env python3
"""Mac correctness/plumbing only. Never print or compare measured durations."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[4]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def science(r):
    def clean(x):
        if isinstance(x, dict):
            return {k: clean(v) for k, v in x.items() if k not in ('timings', 'wall_s')}
        if isinstance(x, list):
            return [clean(v) for v in x]
        return x
    return clean({k: r[k] for k in ('total_energy', 'energies', 'eigenvalues', 'occupations', 'certified', 'verified', 'verification') } | {
        'convergence': {k: r['convergence'][k] for k in ('converged', 'iterations', 'history', 'last_iteration')}})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--upf', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--cases', default='si8-gamma,si32-gamma,si-skew,si-slab,si-k222s')
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fast = load('screen', ROOT/'scripts/fast_screen.py')
    parser = load('lane', ROOT/'scripts/lane_verify_parse.py')
    rows = []
    for case in a.cases.split(','):
        directory = a.out/case
        directory.mkdir(exist_ok=True)
        src = directory/'input.toml'
        src.write_text(fast.case_input(fast.load_case(ROOT/f'benchmarks/cases/{case}.toml'),
                                      a.upf.resolve(), 601, 1e-3, 1, []))
        reference = None
        for arm, binary, lane, phase, clocks in (
            ('parent', a.baseline, 'off', 'off', 'off'),
            ('off', a.binary, 'off', 'off', 'off'),
            ('phase-only', a.binary, 'off', 'on', 'on'),
            ('count', a.binary, 'on', 'on', 'off'),
            ('clock', a.binary, 'on', 'on', 'on')):
            env = dict(os.environ, BAANDI_VERIFY_PROFILE=lane, BAANDI_PHASE_PROFILE=phase,
                       BAANDI_STAGE_TIMERS=clocks, BAANDI_PHASE_WAIT='off')
            dest = directory/arm
            proc = subprocess.run([str(binary.resolve()), 'run', str(src.resolve()), '--output-dir', str(dest.resolve())],
                                  env=env, capture_output=True, text=True)
            (directory/f'{arm}.stdout').write_text(proc.stdout)
            (directory/f'{arm}.stderr').write_text(proc.stderr)
            assert proc.returncode in (0, 3), (case, arm, proc.returncode, proc.stderr[-400:])
            result = json.loads((dest/'results.json').read_text())
            state = science(result), (dest/'density.bin').read_bytes(), proc.returncode
            if reference is None:
                reference = state
            assert state == reference, (case, arm, 'scientific identity failure')
            if lane == 'on':
                calls, attempts = parser.parse(dest/'results.json')
                if result['verification']['attempts']:
                    assert calls
                assert all(r['ranks_present'] == r['expected_ranks'] for r in attempts)
                assert all(r['timed'] == (clocks == 'on') for r in calls)
                assert len([r for r in calls if r['event'] == 'attempt' and r['cutoff_ha'] == result['scf']['cutoff_hartree']]) == len(result['verification']['attempts'])
            else:
                assert 'verification_lane' not in result.get('phase_profile', {})
            rows.append(dict(case=case, arm=arm, exit_code=proc.returncode, energy_ha=result['total_energy']['ha'],
                             density_sha256=hashlib.sha256(state[1]).hexdigest(),
                             iterations=result['convergence']['iterations'], attempts=len(result['verification']['attempts']),
                             identity=True))
            (a.out/'identity.json').write_text(json.dumps(rows, indent=2)+'\n')
        print(f'{case}: parent/off/phase-only/count/clock scientific state and density identical; exit {reference[2]}', flush=True)
    (a.out/'identity.json').write_text(json.dumps(rows, indent=2)+'\n')

if __name__ == '__main__':
    main()
