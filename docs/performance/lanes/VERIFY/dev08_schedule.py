#!/usr/bin/env python3
"""Dev08-only crossed A/A or A/B, one case/model per job; no automatic adoption."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import statistics
import subprocess
import time


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('binary', 'input', 'upf', 'out'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--launcher', required=True)
    p.add_argument('--mode', choices=('aa', 'ab'), required=True)
    p.add_argument('--pairs', type=int, default=12)
    a = p.parse_args()
    if platform.system() != 'Linux':
        p.error('timing is dev08-only')
    if a.pairs < 2:
        p.error('at least two pairs required')
    a.out.mkdir(parents=True, exist_ok=False)
    src = a.out/'input.toml'; src.write_bytes(a.input.read_bytes())
    data = dict(binary_sha256=sha(a.binary), input_sha256=sha(src), upf_sha256=sha(a.upf),
                mode=a.mode, launcher=shlex.split(a.launcher), runs=[],
                inherited_environment={k:v for k,v in os.environ.items() if k.startswith(('BAANDI_', 'MKL_', 'OMP_'))})
    for pair in range(-1, a.pairs):
        for arm in ('A', 'B') if pair % 2 == 0 else ('B', 'A'):
            dest = a.out/f'{pair}-{arm}'
            mode = 'first' if arm == 'B' and a.mode == 'ab' else 'on'
            env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                       BAANDI_VERIFY_SCHEDULE=mode, BAANDI_STOP_PREDICTOR='off',
                       BAANDI_VERIFY_PROFILE='off', BAANDI_PHASE_PROFILE='off',
                       BAANDI_STAGE_TIMERS='off', BAANDI_DAVIDSON_STAGES='off', BAANDI_PHASE_WAIT='off')
            cmd = shlex.split(a.launcher)+[str(a.binary.resolve()), 'run', str(src.resolve()), '--output-dir', str(dest.resolve())]
            load = os.getloadavg(); start = time.perf_counter()
            proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
            wall = time.perf_counter()-start
            (a.out/f'{pair}-{arm}.stdout').write_text(proc.stdout)
            (a.out/f'{pair}-{arm}.stderr').write_text(proc.stderr)
            row = dict(pair=pair, arm=arm, retained=pair>=0, command=cmd, exit_code=proc.returncode,
                       external_wall_s=wall, load_before=load, load_after=os.getloadavg())
            if (dest/'results.json').exists():
                r = json.loads((dest/'results.json').read_text()); attempts = r['verification']['attempts']
                w = r['convergence']['work']; ladder = w.get('cutoff_ladder') or {}
                row.update(energy_ha=r['total_energy']['ha'], verified=r['verified'],
                           maps=r['convergence']['iterations'], attempts=len(attempts),
                           failed=sum(not x['passed'] for x in attempts),
                           failed_h_vectors=sum(x['h_vectors'] for x in attempts if not x['passed']),
                           failed_davidson_iterations=sum(x['davidson_iterations'] for x in attempts if not x['passed']),
                           total_h_vectors=w['h_vectors_recorded']+w['certificate_h_vectors']+ladder.get('coarse_h_vectors',0)+ladder.get('coarse_start_h_vectors',0),
                           density_sha256=sha(dest/'density.bin'), results_sha256=sha(dest/'results.json'))
            data['runs'].append(row)
            (a.out/'manifest.json').write_text(json.dumps(data, indent=2)+'\n')
            if proc.returncode or not row.get('verified'):
                raise SystemExit('Failed run; retained all evidence, no adoption verdict')
    deltas = []
    for pair in range(a.pairs):
        rows = {r['arm']:r for r in data['runs'] if r['pair'] == pair}
        deltas.append(rows['B']['external_wall_s']-rows['A']['external_wall_s'])
    mean = statistics.mean(deltas); se = statistics.stdev(deltas)/len(deltas)**0.5
    data['paired_wall'] = dict(deltas_s=deltas, mean_s=mean, se_s=se,
                               abs_mean_over_se=abs(mean)/se if se else None)
    (a.out/'manifest.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(data['paired_wall']))

if __name__ == '__main__':
    main()
