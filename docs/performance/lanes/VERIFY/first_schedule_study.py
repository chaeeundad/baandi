#!/usr/bin/env python3
"""Correctness and work counters only; never interpret Mac durations."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


def prediction(history, threshold):
    below = lambda r: (r.get('standard_error_ha') or 0) > 0 and r['standard_error_ha'] < threshold
    if len(history) < 2 or not below(history[-1]):
        return 0
    crossings = [i for i, r in enumerate(history) if below(r)]
    if not crossings or crossings[0] == 0:
        return 0
    anchor = crossings[0]
    elapsed = history[-1]['iteration'] - history[anchor]['iteration']
    if elapsed >= 2 or any(r['mixer_reset'] or r['mixing_raw_recovery'] for r in history[anchor-1:]):
        return 0
    old, now = history[anchor-1:anchor+1]
    if (old['mixing_history_used'] or 0) < 2:
        return 0
    q = max(now['standard_error_ha']/old['standard_error_ha'], (now['density_residual']/old['density_residual'])**2)
    if not 0 < q < 1:
        return 0
    for n in range(3-elapsed):
        if old['standard_error_ha'] / (1-q**0.5)**2 * q**(elapsed+n) < threshold:
            return n
    return 2-elapsed


def features(r):
    a = r['verification']['attempts'][0]
    h = [x for x in r['convergence']['history'] if x['iteration'] <= a['iteration']]
    x, prev = h[-1], h[-2]
    e = x['standard_error_ha']
    return dict(iteration=a['iteration'], passed=a['passed'], threshold_ha=a['threshold_ha'],
                error_ha=e, verification_error_ha=a['coulomb_error_ha'], amplification=a['coulomb_error_ha']/e,
                density_residual=x['density_residual'], coulomb_ratio=e/prev['standard_error_ha'],
                density_squared_ratio=(x['density_residual']/prev['density_residual'])**2,
                eigen_tolerance_ha=x['requested_eigen_tolerance_ha'], electrons=r['electrons'],
                volume_bohr3=r.get('cell', {}).get('volume_bohr3'),
                predicted_maps=prediction(h, a['threshold_ha']))


def clean(x):
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items() if k not in ('timings', 'wall_s', 'schedule')}
    if isinstance(x, list):
        return [clean(v) for v in x]
    return x


def science(r):
    return clean({k: r.get(k) for k in ('total_energy', 'energies', 'eigenvalues', 'occupations', 'certified', 'verified', 'verification')} | {
        'history': r['convergence']['history'], 'converged': r['convergence']['converged']})


def counters(r):
    w = r['convergence']['work']
    attempts = r['verification']['attempts']
    # Recorded plus ladder's own total, as used by the existing screen.
    ladder = w.get('cutoff_ladder') or {}
    return dict(maps=r['convergence']['iterations'], attempts=len(attempts),
                failed=sum(not a['passed'] for a in attempts),
                h_recorded=w['h_vectors_recorded'],
                h_total=w['h_vectors_recorded'] + w['certificate_h_vectors'] + ladder.get('coarse_h_vectors', 0) + ladder.get('coarse_start_h_vectors', 0),
                h_failed=sum(a['h_vectors'] for a in attempts if not a['passed']),
                davidson_failed=sum(a['davidson_iterations'] for a in attempts if not a['passed']),
                energy_ha=r['total_energy']['ha'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--records', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--binary', type=Path)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--analyze-only', action='store_true')
    p.add_argument('--candidate-tag', default='', help='fresh trace/first output suffix after a revision')
    p.add_argument('--cases', default='si8,si16,si32,si-skew,mgo8,syn-gaussian-fcc,syn-p3,syn-spd')
    a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in a.cases.split(','):
        for seed in (601, 602, 603):
            strict = None
            for threshold in ('5e-09', '1e-10'):
                name = f'{case}-{seed}-{threshold}'
                source = a.records/name
                historical = json.loads((source/'on/results.json').read_text())
                row = dict(case=case, seed=seed, threshold_ha=float(threshold), historical=features(historical))
                if not a.analyze_only:
                    dest = a.out/name; dest.mkdir(exist_ok=True)
                    src = dest/'input.toml'; src.write_text((source/'input.toml').read_text())
                    def run(arm, binary, mode, reference=False):
                        target = dest/(arm + a.candidate_tag if arm in ('trace', 'first') else arm)
                        if not (target/'exit.json').exists():
                            target.mkdir(exist_ok=True)
                            env = {k:v for k,v in os.environ.items() if not k.startswith('BAANDI_')}
                            env.update(BAANDI_VERIFY_SCHEDULE=mode, BAANDI_STAGE_TIMERS='off', BAANDI_PHASE_PROFILE='off',
                                       BAANDI_REFERENCE_STRICT='on' if reference else 'off', TMPDIR=str(a.out.resolve()))
                            binary_hash = hashlib.sha256(binary.read_bytes()).hexdigest()
                            cmd = [str(binary.resolve()), 'run', str(src.resolve()), '--output-dir', str(target.resolve())]
                            with (target/'stdout').open('w') as out, (target/'stderr').open('w') as err:
                                proc = subprocess.run(cmd, env=env, stdout=out, stderr=err, timeout=300)
                            (target/'exit.json').write_text(json.dumps(dict(command=cmd, exit_code=proc.returncode,
                                binary_sha256=binary_hash)))
                        saved = json.loads((target/'exit.json').read_text())
                        assert saved['exit_code'] == 0, str(target)
                        if arm in ('trace', 'first'):
                            assert saved['binary_sha256'] == hashlib.sha256(binary.read_bytes()).hexdigest(), ('stale candidate cache; use a new --candidate-tag', str(target))
                        result = json.loads((target/'results.json').read_text())
                        assert result['convergence']['converged'], str(target)
                        assert result['single_path']['verification'] == 'strict_certificate' if reference else result['verified'], str(target)
                        if reference:
                            certificate = result['convergence']['last_iteration']['fixed_point_certificate']
                            assert certificate['occupied_eigenvalue_change_ha'] <= 1e-7
                            assert max(v for block in certificate['strict_residuals_per_k_ha'] for v in block) <= 1e-9
                        if not reference:
                            last = result['verification']['attempts'][-1]
                            assert last['passed'] and last['coulomb_error_ha'] < last['threshold_ha']
                            assert last['residual_bound_ha'] is None or last['residual_bound_ha'] < last['threshold_ha']
                        return result
                    if strict is None:
                        strict = run('strict', a.baseline, 'on', True)
                        strict_b = run('strict-first', a.binary, 'first', True)
                        assert science(strict) == science(strict_b), 'strict identity'
                    baseline = run('parent', a.baseline, 'on')
                    default = run('default', a.binary, 'on')
                    trace = run('trace', a.binary, 'first-trace')
                    candidate = run('first', a.binary, 'first')
                    for arm, result in [('default', default), ('trace', trace)]:
                        assert science(baseline) == science(result), (name, arm, 'science identity')
                        assert (dest/'parent/density.bin').read_bytes() == (dest/(arm+a.candidate_tag if arm == 'trace' else arm)/'density.bin').read_bytes()
                    row['current'] = features(baseline)
                    forecast = trace['verification']['schedule']['records'][0]
                    assert forecast['predicted_additional_maps'] == row['current']['predicted_maps'], name
                    row['skipped_failure'] = bool(forecast['would_wait'] and not row['current']['passed'])
                    row['wrong_skip'] = bool(forecast['would_wait'] and row['current']['passed'])
                    row['baseline'], row['candidate'] = counters(baseline), counters(candidate)
                    row['error_vs_strict_ha'] = abs(candidate['total_energy']['ha'] - strict['total_energy']['ha'])
                    row['baseline_error_vs_strict_ha'] = abs(baseline['total_energy']['ha'] - strict['total_energy']['ha'])
                    row['delta_energy_ha'] = candidate['total_energy']['ha'] - baseline['total_energy']['ha']
                    row['accuracy_gate'] = row['error_vs_strict_ha'] <= float(threshold)
                    row['trace_and_default_identity'] = True
                rows.append(row)
                (a.out/'summary.json').write_text(json.dumps(rows, indent=2)+'\n')
                assert row.get('accuracy_gate', True), (name, 'accuracy gate')
                print(name, 'historical fail', not row['historical']['passed'], 'predict', row['historical']['predicted_maps'],
                      'gate', row.get('accuracy_gate', 'not-run'), flush=True)

if __name__ == '__main__':
    main()
