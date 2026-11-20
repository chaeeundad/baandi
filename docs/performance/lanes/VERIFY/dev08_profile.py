#!/usr/bin/env python3
"""Dev08-only Phase 0 runner: one case/model per invocation. No adoption verdict."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--upf', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--launcher', required=True, help='e.g. taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core')
    p.add_argument('--repeats', type=int, default=6)
    a = p.parse_args()
    if platform.system() != 'Linux':
        p.error('profiling is dev08-only; use validate.py for Mac correctness')
    if a.repeats < 1:
        p.error('repeats must be positive')
    a.out.mkdir(parents=True, exist_ok=False)
    src = a.out/'input.toml'
    src.write_bytes(a.input.read_bytes())
    manifest = dict(binary=str(a.binary.resolve()), binary_sha256=sha(a.binary),
                    input=str(a.input.resolve()), input_sha256=sha(src), upf_sha256=sha(a.upf),
                    launcher=shlex.split(a.launcher), machine=platform.node(),
                    inherited_environment={k:v for k,v in os.environ.items() if k.startswith(('BAANDI_', 'MKL_', 'OMP_'))},
                    runs=[])
    modes = {'off': ('off','off','off'), 'counts': ('on','on','off'), 'clocks': ('on','on','on')}
    for repeat in range(-1, a.repeats):
        order = ['off','counts','clocks'] if repeat % 2 == 0 else ['clocks','counts','off']
        for arm in order:
            lane,phase,clocks = modes[arm]
            dest = a.out/f'{"warmup" if repeat < 0 else repeat}-{arm}'
            env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                       BAANDI_VERIFY_PROFILE=lane, BAANDI_PHASE_PROFILE=phase,
                       BAANDI_STAGE_TIMERS=clocks, BAANDI_DAVIDSON_STAGES=clocks,
                       BAANDI_PHASE_WAIT='off')
            cmd = shlex.split(a.launcher) + [str(a.binary.resolve()), 'run', str(src.resolve()), '--output-dir', str(dest.resolve())]
            load = os.getloadavg()
            started = time.perf_counter()
            proc = subprocess.run(cmd, env=env, capture_output=True, text=True)
            wall = time.perf_counter()-started
            (a.out/f'{dest.name}.stdout').write_text(proc.stdout)
            (a.out/f'{dest.name}.stderr').write_text(proc.stderr)
            row = dict(repeat=repeat, arm=arm, retained=repeat>=0, command=cmd,
                       load_before=load, load_after=os.getloadavg(), external_wall_s=wall, exit_code=proc.returncode)
            if (dest/'results.json').exists():
                row['results_sha256']=sha(dest/'results.json')
            if (dest/'density.bin').exists():
                row['density_sha256']=sha(dest/'density.bin')
            manifest['runs'].append(row)
            (a.out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
            if proc.returncode:
                raise SystemExit(f'run failed: {dest}; exit {proc.returncode}; preserve all output')
    print(a.out/'manifest.json')

if __name__=='__main__':
    main()
