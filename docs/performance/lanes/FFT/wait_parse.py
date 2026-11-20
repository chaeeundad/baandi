"""Parse the proposed arrival-skew-1 trace. No production trace producer exists yet."""
import argparse
import json
import math
import statistics
from pathlib import Path


def nonnegative(x):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0:
        raise ValueError('invalid nonnegative metric')
    return x


def parse(doc):
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


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('trace', type=Path)
    args = ap.parse_args()
    print(json.dumps(parse(json.loads(args.trace.read_text())), indent=2, allow_nan=False))
