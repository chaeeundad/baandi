#!/usr/bin/env python3
"""Validate rotation routes and model MKL calls/copies. No timing verdict."""
import collections
import json
from pathlib import Path
import sys


def counters(path):
    r = json.loads(Path(path).read_text())
    if not r['run']['completed']:
        raise ValueError('incomplete run')
    p = r['davidson_profile']
    result = []
    for rank in p['per_rank']:
        if rank['dropped']:
            raise ValueError('dropped records')
        maps = collections.defaultdict(lambda: dict(pairs=0, column_pairs=0, flop=0,
                                                    row_calls_model=0, group_calls_model=0,
                                                    chosen_calls_model=0, extra_copy_bytes_model=0))
        for row in rank['records']:
            d = dict(zip(p['fields'], row))
            if p['operations'][int(d['operation'])] != 'rotation_pair':
                continue
            n, m, k, route = [int(d[x]) for x in ('n', 'm', 'k', 'detail')]
            if min(n, m, k) <= 0 or route not in (0, 1):
                raise ValueError('invalid rotation shape or route')
            v = maps[int(d['map'])]
            v['pairs'] += 1
            v['column_pairs'] += route
            v['flop'] += 4 * n * m * k
            v['row_calls_model'] += 2 * ((n + 255) // 256)
            v['group_calls_model'] += 2 * ((k + 31) // 32)
            v['chosen_calls_model'] += 2 * ((k + 31) // 32) if route else 2 * ((n + 255) // 256)
            v['extra_copy_bytes_model'] += 0 if route else 32 * n * k
        total = {key: sum(v[key] for v in maps.values()) for key in
                 ('pairs', 'column_pairs', 'flop', 'row_calls_model', 'group_calls_model',
                  'chosen_calls_model', 'extra_copy_bytes_model')}
        result.append(dict(rank=rank['rank'], maps=dict(maps), total=total,
                           active=bool(total['pairs'])))
    return dict(path=str(path), ranks=result,
                note='Call/copy models apply to MKL DGEMM routes, not pure kernels; '
                     'confirm gemm_dispatch route=mkl and compare rotate/other FLOP.')


if __name__ == '__main__':
    print(json.dumps([counters(path) for path in sys.argv[1:]], indent=2))
