#!/usr/bin/env python3
"""Parse verify-profile-1 results into per-call/per-rank and critical-rank tables.
No optimization verdict or cross-parallel-model comparison is made.
"""
import argparse
import json
import math
from pathlib import Path
import sys

PARTS = ('owned_input', 'warm_start', 'map', 'coulomb', 'metrics', 'decision', 'record', 'restore', 'schedule')
SIZES = ('grid_x', 'grid_y', 'grid_z', 'full_density_elements', 'owned_density_elements',
         'warm_complex_bytes', 'k_blocks', 'lazy_blocks', 'bands_sum', 'plane_waves_sum', 'history_records', 'previous_attempts')
MAP_PARTS = ('potential', 'solve', 'density', 'energy')


def number(value, name, integer=False):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
        raise ValueError(f'{name}: expected a finite nonnegative number')
    if integer and (not isinstance(value, int)):
        raise ValueError(f'{name}: expected an integer')
    return value


def parse(path):
    data = json.loads(Path(path).read_text())
    lane = data.get('phase_profile', {}).get('verification_lane')
    if not isinstance(lane, dict) or lane.get('schema') != 'verify-profile-1':
        raise ValueError(f'{path}: missing verify-profile-1 (enable VERIFY_PROFILE and PHASE_PROFILE)')
    ranks = number(data['phase_profile']['ranks'], 'ranks', True)
    rows = []
    occurrences = {}
    for record in lane['records']:
        if record['event'] not in ('attempt', 'schedule'):
            raise ValueError('invalid event')
        rank = number(record['rank'], 'rank', True)
        if rank >= ranks:
            raise ValueError('rank exceeds reported communicator size')
        iteration = number(record['iteration'], 'iteration', True)
        cutoff = number(record['cutoff_ha'], 'cutoff_ha')
        if record['status'] not in (-1, 0, 1) or isinstance(record['status'], bool):
            raise ValueError('invalid attempt status')
        if not isinstance(record['timed'], bool):
            raise ValueError('timed must be boolean')
        sizes = {k: number(record['sizes'][k], k, True) for k in SIZES}
        seconds = {k: number(record['seconds'][k], k) for k in PARTS}
        nested = {k: number(record['map_seconds'][k], k) for k in MAP_PARTS}
        if not record['timed'] and any((*seconds.values(), *nested.values())):
            raise ValueError('count-only record contains clock values')
        key = (record['event'], rank, cutoff, iteration)
        occurrence = occurrences.get(key, 0)
        occurrences[key] = occurrence + 1
        rows.append(dict(source=str(path), event=record['event'], rank=rank, iteration=iteration, cutoff_ha=cutoff,
                         occurrence=occurrence, status=record['status'], timed=record['timed'],
                         sizes=sizes, seconds=seconds, map_seconds=nested,
                         total_s=sum(seconds.values()),
                         h_vectors=number(record['h_vectors'], 'h_vectors', True),
                         davidson_iterations=number(record['davidson_iterations'], 'davidson_iterations', True),
                         warm_copy_read_write_bytes=2 * sizes['warm_complex_bytes'],
                         owned_copy_read_write_bytes=16 * sizes['owned_density_elements'],
                         density_pair_payload_bytes=16 * (sizes['owned_density_elements'] or sizes['full_density_elements'])))
    groups = {}
    for row in rows:
        groups.setdefault((row['event'], row['cutoff_ha'], row['iteration'], row['occurrence']), []).append(row)
    attempts = []
    for (event, cutoff, iteration, occurrence), group in groups.items():
        critical = max(group, key=lambda r: r['total_s'])
        attempts.append(dict(source=str(path), event=event, cutoff_ha=cutoff, iteration=iteration, occurrence=occurrence,
                             ranks_present=len(group), expected_ranks=ranks, critical_rank=critical['rank'],
                             max_total_s=critical['total_s'], min_total_s=min(r['total_s'] for r in group),
                             critical_rank_seconds=critical['seconds'],
                             statuses=[r['status'] for r in group]))
    return rows, attempts


def aggregate(rows):
    """Close all lane events on each rank before selecting a run's critical rank."""
    sources = {}
    for row in rows:
        ranks = sources.setdefault(row['source'], {})
        total = ranks.setdefault(row['rank'], dict(rank=row['rank'], events=0,
                                seconds=dict.fromkeys(PARTS, 0.)))
        total['events'] += 1
        for key in PARTS:
            total['seconds'][key] += row['seconds'][key]
    runs = []
    for source, ranks in sources.items():
        totals = list(ranks.values())
        for total in totals:
            total['total_s'] = sum(total['seconds'].values())
        critical = max(totals, key=lambda r: r['total_s'])
        runs.append(dict(source=source, rank_totals=totals, critical_rank=critical['rank'],
                         max_lane_total_s=critical['total_s'], critical_rank_seconds=critical['seconds']))
    return runs


def render(rows, attempts):
    lines = ['# Verification lane profile', '',
             'Seconds are observations from the supplied files. Nested map phases are not added to map wall.',
             'Byte columns are payload/floor inputs, not measured DRAM traffic. Error attempts have incomplete attribution.', '',
             '| File | Event | Rank | Iteration | Cutoff Ha | Status | Timed | Total s | Warm bytes | Owned elements | Bands | H vectors |',
             '|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|']
    for r in rows:
        s = r['sizes']
        lines.append(f"| {r['source']} | {r['event']} | {r['rank']} | {r['iteration']} | {r['cutoff_ha']:g} | {r['status']} | {r['timed']} | {r['total_s']:.9g} | {s['warm_complex_bytes']} | {s['owned_density_elements']} | {s['bands_sum']} | {r['h_vectors']} |")
    lines += ['', '| File / rank / iteration | ' + ' | '.join(PARTS) + ' |', '|---|' + '---:|' * len(PARTS)]
    for r in rows:
        lines.append(f"| {r['source']} / {r['rank']} / {r['iteration']} | " + ' | '.join(f"{r['seconds'][k]:.9g}" for k in PARTS) + ' |')
    lines += ['', '| File / rank / iteration | Map potential s | Map solve s | Map density s | Map energy s |', '|---|---:|---:|---:|---:|']
    for r in rows:
        lines.append(f"| {r['source']} / {r['rank']} / {r['iteration']} | " + ' | '.join(f"{r['map_seconds'][k]:.9g}" for k in MAP_PARTS) + ' |')
    lines += ['', '| File | Iteration | Ranks present/expected | Critical rank | Max total s | Min total s |', '|---|---:|---|---:|---:|---:|']
    for r in attempts:
        lines.append(f"| {r['source']} | {r['iteration']} | {r['ranks_present']}/{r['expected_ranks']} | {r['critical_rank']} | {r['max_total_s']:.9g} | {r['min_total_s']:.9g} |")
    lines += ['', '| File | Run critical rank | Lane total s on that rank |', '|---|---:|---:|']
    for r in aggregate(rows):
        lines.append(f"| {r['source']} | {r['critical_rank']} | {r['max_lane_total_s']:.9g} |")
    lines += ['', 'Critical rank is chosen by the sum of exclusive parts on THAT rank; maxima of different parts are never summed.']
    return '\n'.join(lines) + '\n'


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('files', nargs='+', type=Path)
    p.add_argument('--json-out', type=Path)
    p.add_argument('--markdown-out', type=Path)
    args = p.parse_args(argv)
    rows, attempts = [], []
    try:
        for path in args.files:
            a, b = parse(path)
            rows.extend(a)
            attempts.extend(b)
    except (OSError, ValueError, KeyError, TypeError) as error:
        p.error(str(error))
    if args.json_out:
        args.json_out.write_text(json.dumps(dict(schema='verify-lane-tables-1', calls=rows, attempts=attempts, runs=aggregate(rows)), indent=2) + '\n')
    report = render(rows, attempts)
    if args.markdown_out:
        args.markdown_out.write_text(report)
    else:
        sys.stdout.write(report)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
