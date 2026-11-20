"""Parse production arrival-skew-2 from results.json or DAVBENCH_JSON; retain v1 design fixtures."""
import argparse
import json
import math
import statistics
from pathlib import Path


def nonnegative(x):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
        raise ValueError('invalid nonnegative metric')
    return x


def parse_legacy(doc):
    if doc.get('schema') != 'arrival-skew-1' or not doc.get('complete'):
        raise ValueError('missing schema or incomplete/overflowed trace')
    mode = doc.get('mode')
    if mode not in ('counts', 'timed'):
        raise ValueError('invalid mode')
    ranks = sorted(doc['ranks'], key=lambda r: r['rank'])
    if len(ranks) != doc['rank_count'] or [r['rank'] for r in ranks] != list(range(len(ranks))):
        raise ValueError('missing/duplicate rank')
    if len(ranks) < 2:
        raise ValueError('needs multiple ranks')
    error = nonnegative(doc['clock_error_ns']) if mode == 'timed' else 0
    if mode == 'timed' and doc.get('clock_valid') is not True:
        raise ValueError('uncalibrated clocks')
    keys = [(m['map_id'], m['bucket']) for m in ranks[0]['maps']]
    if len(set(keys)) != len(keys):
        raise ValueError('duplicate map')
    for r in ranks:
        if [(m['map_id'], m['bucket']) for m in r['maps']] != keys:
            raise ValueError('map sequence mismatch')
    out = []
    for mi, (map_id, bucket) in enumerate(keys):
        maps = [r['maps'][mi] for r in ranks]
        signatures = [[(e['seq'], e['op'], e['site']) for e in m['events']] for m in maps]
        if any(s != signatures[0] for s in signatures):
            raise ValueError('collective sequence mismatch')
        if [s[0] for s in signatures[0]] != list(range(len(signatures[0]))):
            raise ValueError('noncontiguous sequence')
        previous = [nonnegative(m['start_ns']) for m in maps] if mode == 'timed' else []
        for i, (_, op, site) in enumerate(signatures[0]):
            if op not in ('alltoall', 'tree_sum', 'allgatherv', 'bcast', 'barrier', 'reduce'):
                raise ValueError('unknown collective')
            events = [m['events'][i] for m in maps]
            payload = [nonnegative(e['payload_bytes']) for e in events]
            row = dict(map_id=map_id, bucket=bucket, seq=i, op=op, site=site,
                       payload_bytes_by_rank=payload, mode=mode)
            if mode == 'counts':
                row.update(last_rank=None, ranks=[])
                out.append(row)
                continue
            enter = [nonnegative(e['enter_ns']) for e in events]
            leave = [nonnegative(e['exit_ns']) for e in events]
            stage_totals = []
            for r, e in enumerate(events):
                if not previous[r] <= enter[r] <= leave[r]:
                    raise ValueError('invalid event interval')
                end = previous[r]
                totals = {}
                for span in e['preceding_stages']:
                    a, b = nonnegative(span['start_ns']), nonnegative(span['end_ns'])
                    if not end <= a <= b <= enter[r]:
                        raise ValueError('overlapping/out-of-window compute spans')
                    end = b
                    totals[span['stage']] = totals.get(span['stage'], 0) + b - a
                    for field in ('points', 'logical_bytes', 'flops'):
                        nonnegative(span[field])
                totals['unattributed'] = enter[r] - previous[r] - sum(totals.values())
                stage_totals.append(totals)
            latest = max(enter)
            candidates = [r for r, a in enumerate(enter) if latest - a <= 2 * error]
            last = candidates[0] if len(candidates) == 1 else None
            row.update(last_rank=last, last_rank_candidates=candidates,
                       arrival_spread_ns=latest-min(enter), clock_error_ns=error, ranks=[])
            for r in range(len(ranks)):
                gap = latest-enter[r]
                row['ranks'].append(dict(rank=r, arrival_gap_ns=gap,
                    arrival_gap_bounds_ns=[max(0, gap-2*error), gap+2*error],
                    collective_residence_ns=leave[r]-enter[r],
                    pre_last_entry_residence_ns=min(gap, leave[r]-enter[r]),
                    previous_release_ns=previous[r], preceding_compute_ns=stage_totals[r],
                    layout=maps[r]['layout'], preceding_work_spans=events[r]['preceding_stages']))
            if last is not None:
                names = set().union(*(s.keys() for s in stage_totals))
                row['last_rank_stage_excess_vs_other_median_ns'] = {
                    name: stage_totals[last].get(name, 0) - statistics.median(
                        s.get(name, 0) for r, s in enumerate(stage_totals) if r != last)
                    for name in sorted(names)}
                early = enter.index(min(enter))
                row['last_vs_earliest_decomposition_ns'] = dict(
                    earlier_rank=early,
                    previous_release=previous[last]-previous[early],
                    local_interval=(enter[last]-previous[last])-(enter[early]-previous[early]))
            previous = leave
            out.append(row)
    return out


def parse(doc):
    if 'phase_profile' in doc:
        doc = (doc['phase_profile'] or {}).get('arrival_trace')
    if not isinstance(doc, dict):
        raise ValueError('missing/disabled arrival trace')
    if doc.get('schema') == 'arrival-skew-1':
        return parse_legacy(doc)
    if doc.get('schema') != 'arrival-skew-2' or not doc.get('enabled') or not doc.get('complete'):
        raise ValueError('disabled, incomplete, overflowed or unknown trace')
    mode = doc['mode']
    if mode not in ('counts', 'timed'):
        raise ValueError('invalid mode')
    timed = mode == 'timed'
    if timed and not doc.get('clock_valid'):
        raise ValueError('invalid clock calibration')
    error = nonnegative(doc['clock_error_ns']) if timed else 0
    ranks = sorted(doc['ranks'], key=lambda r: r['rank'])
    if len(ranks) < 2 or len(ranks) != doc['rank_count'] or [r['rank'] for r in ranks] != list(range(len(ranks))):
        raise ValueError('missing/duplicate ranks')
    signature = lambda e: (e['map_id'], e['bucket'], e['seq'], e['op'], e.get('site'), e['call_index'])
    expected = [signature(e) for e in ranks[0]['events']]
    for r in ranks:
        if [signature(e) for e in r['events']] != expected:
            raise ValueError('collective sequence mismatch')
        if [e['call_index'] for e in r['events']] != list(range(len(expected))):
            raise ValueError('noncontiguous call index')
        previous_map, previous_seq, previous_end = None, None, None
        for e in r['events']:
            if e['seq'] != (previous_seq + 1 if previous_map == e['map_id'] else 0):
                raise ValueError('map call index mismatch')
            previous_map, previous_seq = e['map_id'], e['seq']
            nonnegative(e['payload_bytes'])
            if timed:
                a, b, prev = (nonnegative(e[k]) for k in ('enter_ns', 'exit_ns', 'previous_exit_ns'))
                if not prev <= a <= b or (previous_end is not None and a < previous_end):
                    raise ValueError('invalid collective interval')
                if e['seq'] and prev != previous_end:
                    raise ValueError('previous exit mismatch')
                previous_end = b
                if nonnegative(e['previous_stage_end_ns']) > a:
                    raise ValueError('stage finish after collective entry')
                for key in ('previous_stage_ns', 'previous_stage_end_ns'):
                    nonnegative(e[key])
            elif any(e[k] is not None for k in ('enter_ns', 'exit_ns', 'previous_exit_ns', 'previous_stage_ns', 'previous_stage_end_ns')):
                raise ValueError('count-only trace contains timestamps')
    out = []
    for events in zip(*(r['events'] for r in ranks)):
        e0 = events[0]
        if e0['bucket'] != 'fine':
            continue
        if e0['op'] not in ('alltoall', 'allgatherv', 'tree_sum', 'bcast', 'barrier'):
            raise ValueError('unknown collective')
        row = {k: e0[k] for k in ('map_id', 'bucket', 'seq', 'op', 'site', 'call_index')}
        row.update(mode=mode, ranks=[], last_rank=None)
        for rank, e in enumerate(events):
            row['ranks'].append({k: e[k] for k in
                ('payload_bytes', 'previous_stage', 'previous_stage_points',
                 'previous_stage_logical_bytes', 'previous_stage_shape')})
            row['ranks'][-1]['rank'] = rank
        if timed:
            enter = [e['enter_ns'] for e in events]
            latest = max(enter)
            candidates = [r for r, a in enumerate(enter) if latest-a <= 2*error]
            row.update(last_rank_candidates=candidates,
                       last_rank=candidates[0] if len(candidates) == 1 else None,
                       arrival_spread_ns=latest-min(enter), clock_error_ns=error)
            for r, e in enumerate(events):
                gap, residence = latest-e['enter_ns'], e['exit_ns']-e['enter_ns']
                row['ranks'][r].update(arrival_gap_ns=gap,
                    arrival_gap_bounds_ns=[max(0, gap-2*error), gap+2*error],
                    collective_residence_ns=residence,
                    pre_last_entry_residence_ns=min(gap, residence),
                    preceding_interval_ns=e['enter_ns']-e['previous_exit_ns'],
                    previous_stage_ns=e['previous_stage_ns'],
                    previous_stage_end_ns=e['previous_stage_end_ns'],
                    previous_stage_is_exclusive=False)
            if row['last_rank'] is not None:
                last, early = row['last_rank'], enter.index(min(enter))
                release = events[last]['previous_exit_ns']-events[early]['previous_exit_ns']
                row['last_vs_earliest_decomposition_ns'] = dict(
                    earlier_rank=early, previous_release=release,
                    local_interval=latest-min(enter)-release)
        out.append(row)
    if not out:
        raise ValueError('no ordinary fine-map collectives')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('trace', type=Path)
    args = ap.parse_args()
    print(json.dumps(parse(json.loads(args.trace.read_text())), indent=2, allow_nan=False))
