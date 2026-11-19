#!/usr/bin/env python3
"""Validate captured 2x2 evidence and emit complete cell/per-map tables."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def analyze(root):
    rows, maps, groups, controls = [], [], [], []
    max_error = 0.
    for job in sorted(root.glob('*-[0-9][0-9][0-9]')):
        assert read(job / 'exit.json')['exit_code'] == 0, job
        metadata = read(job / 'input.json')
        states = [read(job / f'state-{i}/state.json') for i in range(3)]
        assert len({tuple(s['fine_grid']) for s in states}) == 1
        assert len({s['fine_plane_waves'] for s in states}) == 1
        assert len({s['common_low_g_plane_waves'] for s in states}) == 1
        for i, state in enumerate(states):
            last = read(job / f'state-{i}/coarse-history.json')[-1]
            expected_seed = max(5e-14, min(5e-3, last['eigen_tolerance_ha'],
                0.1 * last['standard_error_ha'] / max(state['electrons'], 1.)))
            assert expected_seed == state['seed_ha']
            state['last_coarse_standard_error_ha'] = last['standard_error_ha']
            state['last_coarse_residual'] = last['density_residual']
            state['bands'] = (job / f'state-{i}/orbitals.bin').stat().st_size // (16 * state['fine_plane_waves'])
        cells = defaultdict(list)
        for path in sorted(job.glob('pair-*/measurement.json')):
            row = read(path)
            state, seed = row['state_id'], row['seed_id']
            assert row['seed_ha'] == states[seed]['seed_ha']
            assert row['coarse_h'] == states[state]['coarse_h']
            history = read(path.parent / 'history.json')
            assert len(history) == row['fine_maps'] == len(row['map_calls_h_verify_calls_h'])
            assert row['fine_h'] == sum(h['h'] for h in history)
            assert row['fine_h'] == sum(c[1] for c in row['map_calls_h_verify_calls_h'])
            assert row['verification_h'] == sum(c[3] for c in row['map_calls_h_verify_calls_h'])
            assert row['rejected_attempts'] == row['verification_passed'].count(False)
            assert row['verification_passed'][-1]
            assert abs(row['electrons'] - states[state]['electrons']) < 1e-8
            assert row['error_ha'] < 5e-9, path
            max_error = max(max_error, row['error_ha'])
            if state == seed:
                for name in ('density.bin', 'orbitals.bin', 'energy-eigenvalues.bin', 'history.json'):
                    assert (path.parent / name).read_bytes() == (job / f'natural-{state}' / name).read_bytes()
                assert row['diagonal_bits_equal']
            row.update(case=metadata['case'], seed=metadata['seed'], path=str(path.relative_to(root)),
                       total_h=row['coarse_h'] + row['fine_h'] + row['verification_h'])
            row['history_sha256'] = sha(path.parent / 'history.json')
            row['density_sha256'] = sha(path.parent / 'density.bin')
            row['handoff_density_sha256'] = sha(job / f'state-{state}/density.bin')
            row['handoff_orbitals_sha256'] = sha(job / f'state-{state}/orbitals.bin')
            rows.append(row)
            cells[state, seed].append(row)
            for h, calls, wall in zip(history, row['map_calls_h_verify_calls_h'], row['map_wall_s']):
                assert h['charge_relative_error'] < 1e-8
                assert h['input_charge_relative_error'] < 1e-8
                maps.append(dict(case=row['case'], seed=row['seed'], pair=row['pair'], order=row['order'],
                    state_id=state, seed_id=seed, **h, h_calls=calls[0], verification_calls=calls[2],
                    verification_h=calls[3], wall_s=wall,
                    charge_absolute_error_electrons=h['charge_relative_error'] * states[state]['electrons']))
        assert set(cells) == {(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)}
        summary = dict(case=metadata['case'], seed=metadata['seed'], states=states, cells={})
        for (state, seed), runs in sorted(cells.items()):
            for key in ('fine_maps', 'fine_h', 'verification_h', 'density_sha256', 'history_sha256'):
                assert len({r[key] for r in runs}) == 1, (job, state, seed, key)
            summary['cells'][f'{state}{seed}'] = {
                key: statistics.median(r[key] for r in runs) for key in (
                    'fine_maps', 'fine_h', 'verification_h', 'total_h', 'fine_s', 'full_price_s',
                    'fine_setup_s', 'state_construction_s', 'verification_s', 'error_ha', 'rejected_attempts')}
            summary['cells'][f'{state}{seed}']['fine_range_s'] = [min(r['fine_s'] for r in runs), max(r['fine_s'] for r in runs)]
        c = summary['cells']
        summary['interaction_maps'] = c['11']['fine_maps'] - c['10']['fine_maps'] - c['01']['fine_maps'] + c['00']['fine_maps']
        summary['interaction_fine_h'] = (c['11']['fine_h'] + c['11']['verification_h'] - c['10']['fine_h'] - c['10']['verification_h']
            - c['01']['fine_h'] - c['01']['verification_h'] + c['00']['fine_h'] + c['00']['verification_h'])
        aa = []
        for pair in sorted({r['pair'] for r in rows if r['case'] == metadata['case'] and r['seed'] == metadata['seed']}):
            rs = sorted((r for r in cells[0, 0] if r['pair'] == pair), key=lambda r:r['order'])
            assert len(rs) == 2
            aa.append((rs[1]['fine_s'] / rs[0]['fine_s'] - 1.) * 100)
        summary['aa_fine_percent'] = aa
        groups.append(summary)
        for i in range(3):
            controls.append(dict(case=metadata['case'], seed=metadata['seed'], state_id=i,
                **read(job / f'natural-{i}/measurement.json')))
    assert groups
    result = dict(groups=groups, controls=controls, max_error_ha=max_error,
                  rows=len(rows), maps=len(maps), diagonal_bit_comparisons=sum(r['state_id']==r['seed_id'] for r in rows))
    (root / 'analysis.json').write_text(json.dumps(result, indent=2) + '\n')
    (root / 'rows.json').write_text(json.dumps(rows, indent=2) + '\n')
    with (root / 'maps.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(maps[0]))
        writer.writeheader()
        writer.writerows(dict(m, pulay=json.dumps(m["pulay"])) for m in maps)
    lines = ['| System / seed | State / seed source | Fine maps | Hc / Hf / Hv | Total H | State cost (s) | Fine setup (s) | Fine incl. verify (s) | Full price (s) | Abs. error (Ha) |',
             '| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for g in groups:
        for key, c in g['cells'].items():
            state, seed = map(int, key)
            tag = ('tight' if g['states'][seed]['seed_ha'] == min(s['seed_ha'] for s in g['states'][:2]) else 'loose') if seed < 2 else 'natural'
            lines.append(f"| {g['case']}/{g['seed']} | {g['states'][state]['fraction']:.2f} / {g['states'][seed]['fraction']:.2f} ({tag}) | {c['fine_maps']:g} | {g['states'][state]['coarse_h']} / {c['fine_h']:g} / {c['verification_h']:g} | {c['total_h']:g} | {c['state_construction_s']:.4f} | {c['fine_setup_s']:.4f} | {c['fine_s']:.4f} | {c['full_price_s']:.4f} | {c['error_ha']:.3e} |")
    (root / 'cells.md').write_text('\n'.join(lines)+'\n')
    # One trajectory per cell; all A/A and reversal repetitions are in maps.csv.
    lines = ['| System/seed | State/seed source | Map | Effective tol (Ha) | Redo | H / blocks | Hv / blocks | E_H residual (Ha) | Charge error (e) | Reset / refit | Rejected verifies |',
             '| --- | --- | ---: | ---: | --- | --- | --- | ---: | ---: | --- | ---: |']
    seen = set()
    for row in rows:
        key = row['case'], row['seed'], row['state_id'], row['seed_id']
        if key in seen: continue
        seen.add(key)
        for h in (m for m in maps if (m['case'],m['seed'],m['state_id'],m['seed_id'],m['pair'],m['order']) == (*key,row['pair'],row['order'])):
            rejected = sum(n == h['map'] and not passed for n, passed in zip(row['verification_maps'],row['verification_passed']))
            refit = bool(h['pulay'] and h['pulay']['coulomb_refit'])
            lines.append(f"| {row['case']}/{row['seed']} | {row['state_id']}/{row['seed_id']} | {h['map']} | {h['eigen_tolerance_ha']:.8e} | {h['redo']} | {h['h']}/{h['h_calls']} | {h['verification_h']}/{h['verification_calls']} | {h['standard_error_ha']:.8e} | {h['charge_absolute_error_electrons']:.3e} | {h['reset']}/{refit} | {rejected} |")
    (root / 'maps.md').write_text('\n'.join(lines)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if not root.is_relative_to((ROOT / 'target/scratch').resolve()):
        parser.error('root must be in target/scratch')
    result = analyze(root)
    print(json.dumps({k:v for k,v in result.items() if k not in ('groups','controls')}, indent=2))
