#!/usr/bin/env python3
"""Full-SCF correctness only; deliberately does not measure elapsed time."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('screen', ROOT / 'scripts/fast_screen.py')
screen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(screen)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def science(v):
    if isinstance(v, dict):
        return {k:science(x) for k,x in v.items() if k not in ('timings', 'wall_s')}
    return [science(x) for x in v] if isinstance(v, list) else v


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--upf', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--json-out', type=Path, required=True)
    p.add_argument('--cases', nargs='+', default=['si8-gamma', 'si32-gamma', 'si-skew'])
    a = p.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)
    records = []
    failures = []
    for case in a.cases:
        inp = a.work / (case + '.toml')
        inp.write_text(screen.case_input(screen.load_case(ROOT / f'benchmarks/cases/{case}.toml'),
                                         a.upf.resolve(), 601, 1e-3, 1, []))
        results = {}
        for mode in ['baseline', 'off', 'on', 'strict-off', 'strict-on']:
            out = (a.work / case / mode).resolve()
            binary = a.baseline if mode == 'baseline' else a.binary
            env = dict(os.environ, BAANDI_ATOMIC_LOCAL_BLOCK='on' if mode.endswith('on') else 'off',
                       BAANDI_REFERENCE_STRICT='on' if mode.startswith('strict') else 'off',
                       BAANDI_INIT_PROFILE='on', BAANDI_STAGE_TIMERS='off', BAANDI_PHASE_PROFILE='off')
            cmd = [str(binary.resolve()), 'run', str(inp.resolve()), '--output-dir', str(out)]
            run = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600)
            out.mkdir(parents=True, exist_ok=True)
            (out / 'stdout.log').write_text(run.stdout)
            (out / 'stderr.log').write_text(run.stderr)
            r = json.loads((out / 'results.json').read_text())
            results[mode] = r
            rec = dict(case=case, mode=mode, command=cmd, exit_code=run.returncode,
                       energy_ha=r['total_energy']['ha'], density_sha256=sha(out / 'density.bin'),
                       converged=r['convergence']['converged'], verified=r.get('verified'),
                       input_sha256=sha(inp), counters={k:v for k,v in r['run']['counters'].items()
                                                      if k.startswith('init_lane.atomic.')})
            records.append(rec)
            if run.returncode != 0 or not rec['converged'] or (not mode.startswith('strict') and not rec['verified']):
                failures.append(f'{case} {mode}: convergence/verification')
            print(f'{case} {mode}: exit {run.returncode}, E={rec["energy_ha"]:.15g} Ha', flush=True)
        subset = {r['mode']:r for r in records if r['case'] == case}
        # The disabled candidate must preserve the full scientific trajectory.
        for key in ['total_energy', 'energies', 'eigenvalues', 'kpoints', 'verification']:
            if science(results['baseline'].get(key)) != science(results['off'].get(key)):
                failures.append(f'{case}: default {key} changed')
        if subset['baseline']['density_sha256'] != subset['off']['density_sha256']:
            failures.append(f'{case}: default density changed')
        # Histories contain diagnostic wall fields: strip them recursively.
        if science(results['baseline']['convergence']['history']) != science(results['off']['convergence']['history']):
            failures.append(f'{case}: default history changed')
        ref = subset['strict-off']['energy_ha']
        for mode in ['off', 'on', 'strict-on']:
            error = abs(subset[mode]['energy_ha'] - ref)
            subset[mode]['error_vs_strict_ha'] = error
            if error > (1e-10 if mode.startswith('strict') else 5e-9):
                failures.append(f'{case} {mode}: accuracy {error} Ha')
        for mode, base in [('on', 'off'), ('strict-on', 'strict-off')]:
            subset[mode]['cross_arm_delta_ha'] = abs(subset[mode]['energy_ha'] - subset[base]['energy_ha'])
            subset[mode]['rounding_gate_1e11'] = subset[mode]['cross_arm_delta_ha'] <= 1e-11
        for mode in ['off', 'on']:
            c = subset[mode]['counters']
            accepted = c.get('init_lane.atomic.local_accepted', 0)
            if (accepted > 0) != (mode == 'on'):
                failures.append(f'{case} {mode}: local dispatch counter')
    a.json_out.write_text(json.dumps(dict(purpose='Mac correctness only; no timing evidence',
        binary_sha256=sha(a.binary), baseline_sha256=sha(a.baseline), upf_sha256=sha(a.upf),
        records=records, failures=failures), indent=2) + '\n')
    if failures:
        raise SystemExit('\n'.join(failures))


if __name__ == '__main__':
    main()
