#!/usr/bin/env python3
"""Count/correctness experiment only. Serial jobs, immutable seeds, no timing claims.

UPFs and large CLI artifacts stay in --work; retained evidence contains no UPF.
Run from any directory; --arms baseline,strict,coulomb,cell,cell-coulomb,strict-cell.
"""
import argparse
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('suite', ROOT / 'scripts/bench_suite.py')
suite = importlib.util.module_from_spec(spec)
spec.loader.exec_module(suite)
CASES = ['si-elongated', 'si-slab', 'si-long4', 'si8-gamma', 'si16-gamma',
         'si32-gamma', 'si-skew', 'syn-gaussian-fcc', 'syn-p3']
ARMS = {
    'baseline': {}, 'strict': {'BAANDI_REFERENCE_STRICT': 'on'},
    'coulomb': {'BAANDI_MIX_COULOMB_GUARD': 'on'},
    'cell': {'BAANDI_MIX_CELL_KERKER': 'on'},
    'cell-coulomb': {'BAANDI_MIX_CELL_KERKER': 'on', 'BAANDI_MIX_COULOMB_GUARD': 'on'},
    'strict-cell': {'BAANDI_REFERENCE_STRICT': 'on', 'BAANDI_MIX_CELL_KERKER': 'on'},
}
SWITCHES = sorted(set(k for a in ARMS.values() for k in a))


def case(name):
    if name != 'si-long4':
        return suite.load_case(ROOT / 'benchmarks/cases' / (name + '.toml'))
    c = copy.deepcopy(suite.load_case(ROOT / 'benchmarks/cases/si8-gamma.toml'))
    c['id'] = name
    c['lattice'][2][2] *= 4
    c['atoms'] = [[a[0], a[1], (a[2] + z) / 4] for z in range(4) for a in c['atoms']]
    c['purpose'] = 'Generated 1x1x4 conventional diamond supercell; 32 atoms'
    return c


def history_summary(rows):
    # A divergent stretch is a maximal consecutive run of increases in the
    # Coulomb standard error, attached to the preceding local minimum.
    stretches, start = [], None
    for i in range(1, len(rows)):
        a, b = rows[i-1].get('standard_error_ha'), rows[i].get('standard_error_ha')
        grows = a is not None and b is not None and b > a
        if grows and start is None:
            start = i - 1
        if start is not None and (not grows or i == len(rows)-1):
            end = i if grows else i - 1
            e0 = rows[start]['standard_error_ha']
            stretches.append({'start_map': start + 1, 'end_map': end + 1,
                              'growth': rows[end]['standard_error_ha'] / e0 if e0 else None})
            start = None
    return {'stretches': stretches,
            'resets': sum(bool(r.get('mixer_reset')) for r in rows),
            'feasibility_damping': sum(r.get('mixing_fraction', 1) < 1 for r in rows),
            'min_predicted_ratio': min((r['pulay']['predicted_ratio'] for r in rows
                                       if r.get('pulay')), default=None),
            'frozen_tolerance_intervals': [i+1 for i in range(1, len(rows))
                if rows[i].get('requested_eigen_tolerance_ha') ==
                rows[i-1].get('requested_eigen_tolerance_ha')],
            'max_density_growth_from_best': max((r['density_residual'] /
                min(x['density_residual'] for x in rows[:i+1])
                for i, r in enumerate(rows) if r['density_residual'] > 0), default=None)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--upf', required=True, type=Path)
    p.add_argument('--work', required=True, type=Path)
    p.add_argument('--evidence', required=True, type=Path)
    p.add_argument('--arms', default='baseline')
    p.add_argument('--cases', default=','.join(CASES))
    p.add_argument('--seeds', default='601,602,603')
    p.add_argument('--amp', type=float, default=1e-3)
    p.add_argument('--max-iterations', type=int, default=100)
    p.add_argument('--timeout', type=int, default=900)
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args()
    args.binary, args.upf = args.binary.resolve(), args.upf.resolve()
    args.work.mkdir(parents=True, exist_ok=True)
    args.evidence.mkdir(parents=True, exist_ok=True)
    binary_hash = suite.file_hash(args.binary)
    for name in args.cases.split(','):
        for seed in map(int, args.seeds.split(',')):
            c = case(name)
            rng = random.Random(seed)
            c['atoms'] = [[(x + rng.uniform(-args.amp, args.amp)) % 1 for x in a] for a in c['atoms']]
            c['max_iterations'] = args.max_iterations
            for arm in args.arms.split(','):
                dest = args.evidence / f'{name}-{seed}-{arm}.json'
                if dest.exists():
                    continue
                work = args.work / f'{name}-{seed}-{arm}'
                work.mkdir(parents=True, exist_ok=True)
                upf = args.upf
                if c['potential'] == 'synthetic':
                    upf = work / 'synthetic.UPF'
                    upf.write_text(suite.synthetic(c))
                inp = suite.case_input(c, upf.resolve(), seed, args.threads)
                if arm in ('strict', 'strict-cell'):
                    inp = inp.replace('energy_tolerance_ha = 1e-8', 'energy_tolerance_ha = 1e-12')
                    inp = inp.replace('density_residual_tolerance = 1e-7', 'density_residual_tolerance = 1e-10')
                    inp += 'eigen_tolerance = 1e-11\n'
                src = work / 'input.toml'
                src.write_text(inp)
                env = {k: v for k, v in os.environ.items() if not k.startswith('BAANDI_')}
                env.update({k: 'off' for k in SWITCHES})
                env.update(ARMS[arm])
                env.update(BAANDI_STOP_PREDICTOR='off', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
                cmd = [str(args.binary), 'run', str(src), '--output-dir', str(work / 'out')]
                try:
                    proc = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=args.timeout)
                    code, stderr = proc.returncode, proc.stderr
                except subprocess.TimeoutExpired:
                    code, stderr = 124, 'Experiment timeout (not a convergence result)'
                (work / 'stderr.txt').write_text(stderr)
                result = work / 'out/results.json'
                failure = work / 'out/failure.json'
                data = json.loads(result.read_text()) if result.exists() else (
                    json.loads(failure.read_text()) if failure.exists() else {})
                hist = data.get('convergence', {}).get('history', [])
                work_counts = data.get('convergence', {}).get('work', {})
                ladder = work_counts.get('cutoff_ladder') or {}
                verification = data.get('verification') or {}
                h_fine = sum(r.get('operator_applications', 0) for r in hist)
                h_cert = sum((r.get('fixed_point_certificate') or {}).get('operator_applications', 0) for r in hist)
                h_ver = sum(r.get('h_vectors', 0) for r in verification.get('attempts', []))
                record = dict(case=name, seed=seed, arm=arm, amplitude_fractional=args.amp,
                    command=cmd, environment={k:v for k,v in env.items() if k.startswith('BAANDI_')},
                    exit_code=code, binary_sha256=binary_hash, upf_sha256=suite.file_hash(upf),
                    input=inp, input_sha256=hashlib.sha256(inp.encode()).hexdigest(),
                    converged=data.get('convergence', {}).get('converged', False),
                    certified=data.get('certified', bool(data.get('convergence', {}).get('converged')
                        and data.get('single_path', {}).get('effective_profile') == 'strict'
                        and hist and hist[-1].get('fixed_point_certificate'))), verified=data.get('verified'),
                    energy_ha=data.get('total_energy', {}).get('ha'),
                    fine_maps=len(hist), fine_h=h_fine, certificate_h=h_cert,
                    verification_h=h_ver, coarse_h=ladder.get('coarse_h_vectors'),
                    total_h=(h_fine+h_cert+ladder.get('coarse_h_vectors', 0)+ladder.get('coarse_start_h_vectors', 0)+data.get('atomic_start', {}).get('initial_h_applications', 0)) if data else None,
                    verification=verification, history=hist, summary=history_summary(hist),
                    threads=args.threads, coarse_start_h=ladder.get('coarse_start_h_vectors', 0),
                    fine_start_h=data.get('atomic_start', {}).get('initial_h_applications', 0),
                    density_sha256=suite.file_hash(work/'out/density.bin') if (work/'out/density.bin').exists() else None,
                    failure=data.get('failure', stderr[-2000:] if code else None))
                dest.write_text(json.dumps(record, indent=2) + '\n')
                print(f'{name} seed={seed} {arm}: exit={code} maps={len(hist)} H={record["total_h"]} E={record["energy_ha"]}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
