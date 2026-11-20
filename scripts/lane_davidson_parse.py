#!/usr/bin/env python3
"""Validate Davidson phase-0 traces and print per-map, per-rank lane tables.

No performance verdict is inferred. Times are inclusive, not additive across
operations. JSON output preserves all calls and shape/floor inputs.
"""
import argparse
import collections
import json
import math
from pathlib import Path

FIELDS = ['call', 'map', 'solve', 'stage', 'kind', 'iteration', 'operation',
          'n', 'm', 'k', 'seconds', 'detail', 'tolerance_ha']


def floor_inputs(op, n, m, k, detail):
    """Lower-bound array accesses, NOT measured memory-controller traffic."""
    flops = None
    array_bytes = None
    if op in ('overlap_border', 'reduced_gram', 'upper_border_tile', 'upper_border_product'):
        flops = 2 * n * m * k
        array_bytes = 8 * (n * m + m * k) if k else 0
    elif op == 'reduced_transform':
        flops = 4 * n * n * k
        array_bytes = 8 * (2 * n * n + n * k) if k else 0
    elif op == 'ritz_coefficients':
        flops = 2 * n * n * k
        array_bytes = 8 * (n * n + 2 * n * k)
    elif op in ('rotate_residual', 'rotation_pair'):
        flops = 4 * n * m * k
        array_bytes = 8 * (2 * n * m + m * k + 2 * n * k)
    elif op == 'residual':
        array_bytes = (8 if detail == 1 else 24) * n * k
    return dict(gemm_flop_model=flops, minimum_array_bytes=array_bytes,
                eigen_order_cubed=n ** 3 if op == 'reduced_eigen' else None,
                reduction_payload_bytes=8 * n if op == 'tree_sum' else None)


def parse(path):
    source = json.loads(Path(path).read_text())
    if not source.get('run', {}).get('completed'):
        raise ValueError(f'{path}: incomplete run report')
    p = source.get('davidson_profile', {})
    if p.get('schema') != 'davidson-profile-1' or p.get('fields') != FIELDS:
        raise ValueError(f'{path}: missing or incompatible Davidson trace')
    ranks = p.get('per_rank', [])
    if not ranks or sorted(r['rank'] for r in ranks) != list(range(len(ranks))):
        raise ValueError(f'{path}: missing or duplicate rank')
    if source.get('run', {}).get('ranks', len(ranks)) != len(ranks):
        raise ValueError(f'{path}: rank count differs from run report')
    calls = []
    for rank in ranks:
        if rank['dropped'] != 0:
            raise ValueError(f'{path}: trace overflow on rank {rank["rank"]}')
        seen = set()
        if not rank['records']:
            raise ValueError(f'{path}: empty rank trace')
        for row in rank['records']:
            if len(row) != len(FIELDS) or any(
                    not isinstance(x, (int, float)) or isinstance(x, bool)
                    or not math.isfinite(x) or x < 0 for x in row):
                raise ValueError(f'{path}: malformed record')
            d = dict(zip(FIELDS, row))
            for key in FIELDS:
                if key not in ('seconds', 'tolerance_ha'):
                    if int(d[key]) != d[key] or d[key] > 2 ** 53:
                        raise ValueError(f'{path}: non-exact integer {key}')
                    d[key] = int(d[key])
            if d['call'] in seen or not d['call'] or not d['solve']:
                raise ValueError(f'{path}: duplicate call or missing solve')
            seen.add(d['call'])
            for key, labels in [('operation', 'operations'), ('stage', 'stages'), ('kind', 'kinds')]:
                try:
                    d[key] = p[labels][d[key]]
                except (KeyError, IndexError) as e:
                    raise ValueError(f'{path}: unknown {key}') from e
            d['rank'] = rank['rank']
            d['reduction_site'] = None
            if d['operation'] == 'tree_sum':
                try:
                    d['reduction_site'] = p['reduction_sites'][d['detail']]
                except IndexError as e:
                    raise ValueError(f'{path}: unknown reduction site') from e
            d.update(floor_inputs(d['operation'], d['n'], d['m'], d['k'], d['detail']))
            calls.append(d)
    # Every kernel is enclosed by a recorded solve on that same rank/map.
    solves = {(r['rank'], r['map'], r['solve']) for r in calls if r['operation'] == 'solve'}
    if any((r['rank'], r['map'], r['solve']) not in solves for r in calls):
        raise ValueError(f'{path}: incomplete solve trace')
    grouped = collections.defaultdict(list)
    for c in calls:
        grouped[(c['rank'], c['map'], c['stage'], c['kind'], c['iteration'],
                 c['operation'], c['reduction_site'])].append(c)
    tables = []
    for key, rows in sorted(grouped.items()):
        rank, map_id, stage, kind, iteration, op, site = key
        table = dict(rank=rank, map=map_id, stage=stage, kind=kind,
                     iteration=iteration, operation=op, reduction_site=site,
                     calls=len(rows), solve_tolerances_ha=sorted(set(r['tolerance_ha'] for r in rows if r['operation'] == 'solve')), inclusive_seconds=sum(r['seconds'] for r in rows),
                     seconds_per_call=sum(r['seconds'] for r in rows) / len(rows),
                     max_call_seconds=max(r['seconds'] for r in rows),
                     shapes=sorted(set((r['n'], r['m'], r['k'], r['detail']) for r in rows)))
        for metric in ('gemm_flop_model', 'minimum_array_bytes',
                       'eigen_order_cubed', 'reduction_payload_bytes'):
            values = [r[metric] for r in rows if r[metric] is not None]
            table[metric] = sum(values) if values else None
        tables.append(table)
    by_map = collections.defaultdict(list)
    for row in tables:
        by_map[(row['map'], row['stage'], row['kind'], row['iteration'],
                row['operation'], row['reduction_site'])].append(row)
    map_table = []
    for key, rows in sorted(by_map.items()):
        map_id, stage, kind, iteration, op, site = key
        times = [r['inclusive_seconds'] for r in rows]
        map_table.append(dict(map=map_id, stage=stage, kind=kind, iteration=iteration,
                              operation=op, reduction_site=site,
                              participating_ranks=[r['rank'] for r in rows],
                              min_rank_seconds=min(times), max_rank_seconds=max(times),
                              min_rank_calls=min(r['calls'] for r in rows),
                              max_rank_calls=max(r['calls'] for r in rows)))
    # Count full EVDs by actual map/solve identity, never by rank-summed time
    # or a guessed number of maps inferred from SCF/collective counts.
    eigen_groups = collections.defaultdict(list)
    for c in calls:
        if c['operation'] in ('solve', 'reduced_eigen_step', 'reduced_eigen', 'bordered_update'):
            eigen_groups[(c['rank'], c['map'], c['stage'], c['kind'], c['iteration'])].append(c)
    eigen_map_work = []
    for key, rows in sorted(eigen_groups.items()):
        rank, map_id, stage, kind, iteration = key
        full = [r for r in rows if r['operation'] == 'reduced_eigen']
        steps = sum(r['operation'] == 'reduced_eigen_step' for r in rows)
        solves = sorted({r['solve'] for r in rows if r['operation'] == 'solve'})
        eigen_map_work.append(dict(rank=rank, map=map_id, stage=stage, kind=kind,
                                   iteration=iteration, solve_ids=solves,
                                   reduced_steps=steps, full_evd_calls=len(full),
                                   bordered_attempts=sum(r['operation'] == 'bordered_update' for r in rows),
                                   full_evd_sum_m3=sum(r['n'] ** 3 for r in full),
                                   full_evd_order_counts=dict(sorted(collections.Counter(r['n'] for r in full).items())),
                                   full_evd_per_solve=len(full) / len(solves) if solves else None))
    return dict(source=str(path), binary_sha256=source.get('run', {}).get('binary', {}).get('sha256'),
                ranks=len(ranks), calls=calls, map_rank_table=tables, map_table=map_table,
                eigen_map_work=eigen_map_work,
                existing_lane_counters=source.get('parallel', {}),
                existing_serial_counters=source.get('convergence', {}).get('work', {}).get('davidson_lean'),
                phase_profile=source.get('phase_profile'),
                note='Inclusive nested timings; never sum operation times into a wall. '
                     'GEMM FLOP models exclude scalar glue. Array bytes are compulsory '
                     'array-access lower bounds, not measured DRAM bytes. Eigen m^3 is '
                     'a work proxy for full EVD calls, not exact FLOP. B5 and eigen-step '
                     'wrappers are separate inclusive operations.')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('results', nargs='+', type=Path)
    ap.add_argument('--json-out', type=Path)
    args = ap.parse_args()
    results = [parse(p) for p in args.results]
    for result in results:
        print(f'\n{result["source"]} ({result["ranks"]} ranks; inclusive timers)')
        print('| map | stage/kind | rank | operation | calls | ms | us/call | shapes |')
        print('|---|---|---|---|---|---|---|---|')
        for r in result['map_rank_table']:
            name = r['operation'] + (':' + r['reduction_site'] if r['reduction_site'] else '')
            print(f'| {r["map"]} | {r["stage"]}/{r["kind"]} | {r["rank"]} | {name} | '
                  f'{r["calls"]} | {1e3*r["inclusive_seconds"]:.6f} | '
                  f'{1e6*r["seconds_per_call"]:.3f} | {r["shapes"]} |')
    if args.json_out:
        args.json_out.write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, TypeError, OSError) as error:
        raise SystemExit(str(error))
