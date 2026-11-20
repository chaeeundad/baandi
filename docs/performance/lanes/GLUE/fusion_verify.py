#!/usr/bin/env python3
"""Correctness and logical-work checks only; never report Mac timing."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bits(v):
    if isinstance(v, float):
        return struct.pack('!d', v).hex()
    if isinstance(v, dict):
        return {k: bits(x) for k, x in v.items() if k not in ('timings', 'wall_s', 'coarse_setup_s', 'coarse_stage_s', 'coarse_wall_fine_step_equivalents')}
    if isinstance(v, list):
        return [bits(x) for x in v]
    return v


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--baseline', type=Path)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--json-out', type=Path, required=True)
    p.add_argument('--allow-nonconverged', action='store_true')
    p.add_argument('--expect-fallback', action='store_true')
    p.add_argument('inputs', type=Path, nargs='+')
    a = p.parse_args()
    records = []
    for src in a.inputs:
        reference = None
        off_counters = None
        for arm in (['baseline'] if a.baseline else []) + ['off', 'on']:
            binary = a.baseline if arm == 'baseline' else a.binary
            out = (a.work / src.stem / arm).resolve()
            env = dict(os.environ, BAANDI_SERIAL_DENSITY_FUSED='on' if arm == 'on' else 'off',
                       BAANDI_GLUE_PROFILE='on', BAANDI_PHASE_PROFILE='on', BAANDI_STAGE_TIMERS='off')
            cmd = [str(binary.resolve()), 'run', str(src.resolve()), '--output-dir', str(out)]
            run = subprocess.run(cmd, env=env, capture_output=True, text=True)
            out.mkdir(parents=True, exist_ok=True)
            (out / 'stdout.log').write_text(run.stdout)
            (out / 'stderr.log').write_text(run.stderr)
            assert run.returncode in ((0, 3) if a.allow_nonconverged else (0,)), (cmd, run.returncode, run.stderr[-1000:])
            r = json.loads((out / 'results.json').read_text())
            signature = bits({k: r.get(k) for k in ['total_energy', 'energies', 'eigenvalues', 'kpoints',
                'occupations', 'work', 'verified', 'verification', 'band_accuracy']})
            signature.update(density=sha(out / 'density.bin'), history=bits(r['convergence']['history']),
                             converged=r['convergence']['converged'], work=bits(r['convergence'].get('work')), exit_code=run.returncode)
            if reference is None:
                reference = signature
            assert signature == reference, (src, arm, 'scientific bits differ')
            if not a.allow_nonconverged:
                assert r.get('verified', True) and r['convergence']['converged']
            c = {k.removeprefix('glue_lane.density.'): v for k, v in r['run']['counters'].items()
                 if k.startswith('glue_lane.density.')}
            if arm == 'off':
                off_counters = c
                assert c['fused_calls'] == c['saved_bytes'] == 0
            if arm == 'on':
                assert c['paired_band_points'] == off_counters['paired_band_points']
                assert c['eligible_calls'] == off_counters['eligible_calls']
                assert c['saved_bytes'] == 16 * c['paired_band_points']
                assert off_counters['logical_bytes'] - c['logical_bytes'] == c['saved_bytes']
                assert (c['fused_calls'] == 0) == a.expect_fallback
                assert c['fused_calls'] == c['eligible_calls']
            records.append(dict(input=str(src), input_sha256=sha(src), arm=arm, command=cmd,
                exit_code=run.returncode, species=r.get('species'), energy_ha=r['total_energy']['ha'], delta_energy_ha=0.,
                density_sha256=signature['density'], work_bits_sha256=hashlib.sha256(json.dumps(signature['work'], sort_keys=True).encode()).hexdigest(), scientific_bits_equal=True,
                converged=r['convergence']['converged'], verified=r.get('verified'), counters=c))
            print(f'{src.stem} {arm}: scientific bits equal; {c}', flush=True)
    a.json_out.write_text(json.dumps(dict(purpose='Mac correctness only; no speed evidence',
        binary_sha256=sha(a.binary), baseline_sha256=sha(a.baseline) if a.baseline else None,
        inputs={str(src): src.read_text() for src in a.inputs}, records=records), indent=2) + '\n')


if __name__ == '__main__':
    main()
