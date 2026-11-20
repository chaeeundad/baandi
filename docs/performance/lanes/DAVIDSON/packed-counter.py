#!/usr/bin/env python3
"""Validate border payload predictions; never infer performance from these counters."""
import json
import sys
from pathlib import Path


def counters(path):
    report = json.loads(Path(path).read_text())
    if not report['run']['completed']:
        raise ValueError('incomplete report')
    trace = report['davidson_profile']
    ranks = trace['per_rank']
    if not ranks or sorted(r['rank'] for r in ranks) != list(range(len(ranks))):
        raise ValueError('missing or duplicate rank')
    result = []
    for rank in ranks:
        if rank['dropped']:
            raise ValueError('dropped records')
        rows = sorted((dict(zip(trace['fields'], r)) for r in rank['records']),
                      key=lambda r: r['call'])
        original = predicted = actual = calls = 0
        pending = None
        for row in rows:
            op = trace['operations'][int(row['operation'])]
            if op == 'border_payload':
                if pending is not None:
                    raise ValueError('missing border collective')
                m, ks, kh, norms = (int(row[k]) for k in ('n', 'm', 'k', 'detail'))
                full = m * (ks + kh + 1)
                packed = sum(j + 1 for j in range(m - ks, m))
                packed += sum(j + 1 for j in range(m - kh, m)) + norms
                pending = (full, packed)
            elif op == 'tree_sum' and pending is not None:
                full, packed = pending
                sent = int(row['n'])
                if sent not in (full, packed):
                    raise ValueError('unexpected collective payload')
                original += 8 * full
                predicted += 8 * packed
                actual += 8 * sent
                calls += 1
                pending = None
        if pending is not None or not calls:
            raise ValueError('no complete border reduction trace')
        result.append(dict(rank=rank['rank'], calls=calls,
                           rectangle_bytes=original, predicted_packed_bytes=predicted,
                           sent_bytes=actual, saved_bytes=original-actual))
    return dict(path=str(path), ranks=result)


if __name__ == '__main__':
    print(json.dumps([counters(p) for p in sys.argv[1:]], indent=2))
