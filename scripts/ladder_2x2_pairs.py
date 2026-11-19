#!/usr/bin/env python3
"""Full-price crossed pairs using seeds observed on this machine/build/geometry.

First run ladder_2x2_run.py on the target host. This driver then invokes the
original offline full-construction driver with its explicit test-only seed
injection. Each arm rebuilds both problems, runs coarse SCF, imports the state,
and converges with unchanged verification. No captured-state setup cost is free.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def arms_for(seeds, pairs):
    arms = ['0.5:1e-4:default']  # unsplit production identity oracle
    manifest = [dict(kind='natural', state=0, seed=0, pair=None)]
    for pair in range(pairs):
        cells = [('aa', 0, 0), ('aa', 0, 0), ('crossed', 0, 0), ('crossed', 0, 1), ('crossed', 1, 0), ('crossed', 1, 1)]
        if pair % 2:
            cells.reverse()
        for order, (kind, state, seed) in enumerate(cells):
            arms.append(f'{[0.5, 0.4][state]}:1e-4:first-seed:{seeds[seed]:.17g}')
            manifest.append(dict(kind=kind, pair=pair, order=order, state=state, seed=seed,
                                 discard_warmup=pair == 0))
    arms.append('0.45:1e-4:default')
    manifest.append(dict(kind='crosscheck', state=2, seed=2, pair=None))
    return ','.join(arms), manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture-root', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--pairs', type=int, default=9)
    p.add_argument('--pseudo', type=Path, default=ROOT.parent / 'pseudo')
    p.add_argument('--timeout', type=int, default=7200)
    args = p.parse_args()
    if args.pairs < 2:
        p.error('at least two pairs, including the discarded warmup, are required')
    job = args.capture_root / f'{args.case}-{args.seed}'
    seeds = [json.loads((job / f'state-{i}/state.json').read_text())['seed_ha'] for i in range(2)]
    provenance = json.loads((args.capture_root / 'provenance.json').read_text())
    import hashlib
    if hashlib.sha256(args.binary.read_bytes()).hexdigest() != provenance['binary_sha256']:
        p.error('capture and paired binary hashes must match; recapture on this host/build')
    arms, manifest = arms_for(seeds, args.pairs)
    cmd = [sys.executable, '-I', str(ROOT / 'scripts/ladder_handoff.py'),
        '--binary', str(args.binary), '--out', str(args.out), '--cases', args.case,
        '--seeds', str(args.seed), '--amp', str(provenance['amplitude_fractional']),
        '--pseudo', str(args.pseudo), '--arms', arms, '--timeout', str(args.timeout)]
    rc = subprocess.run(cmd, cwd=ROOT).returncode
    if rc:
        return rc
    (args.out / 'pair-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    rows = json.loads((args.out / 'rows.json').read_text())
    # The legacy runner sorts path names lexically; recover the integer arm index.
    measurements = [json.loads((args.out / f'{args.case}-{args.seed}' / f'arm-{i}/measurement.json').read_text())
                    for i in range(len(manifest))]
    assert len(measurements) == len(rows) == len(manifest)
    original_input = json.loads((job / 'input.json').read_text())
    measured_job = args.out / f'{args.case}-{args.seed}'
    repeated_input = json.loads((measured_job / 'input.json').read_text())
    assert original_input == repeated_input, 'geometry/UPF identity changed'
    for i, item in enumerate(manifest):
        state_dir = job / f'state-{item["state"]}'
        arm_dir = measured_job / f'arm-{i}'
        for name in ('density.bin', 'orbitals.bin'):
            assert (arm_dir / ('handoff-' + name)).read_bytes() == (state_dir / name).read_bytes(), (i, name)
    pairs = []
    for pair in range(args.pairs):
        data = [(m, r) for m, r in zip(manifest, measurements) if m['pair'] == pair]
        aa = [r for m, r in data if m['kind'] == 'aa']
        cells = {(m['state'],m['seed']):r for m,r in data if m['kind'] == 'crossed'}
        for state in (0,1):
            a, b = cells[state,0], cells[state,1]
            assert a['coarse_h'] == b['coarse_h']
            for row in (a,b):
                assert row['error_ha'] < 5e-9
            pairs.append(dict(pair=pair, state=state, discard_warmup=pair==0,
                seed0_s=a['wall_s'], seed1_s=b['wall_s'],
                seed1_vs_seed0_percent=100*(b['wall_s']/a['wall_s']-1),
                aa_percent=100*(aa[1]['wall_s']/aa[0]['wall_s']-1)))
    (args.out / 'pairs.json').write_text(json.dumps(pairs, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
