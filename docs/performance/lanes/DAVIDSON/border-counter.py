#!/usr/bin/env python3
"""Shape-derived border work only; no timing or performance verdict."""
import json
import sys
from pathlib import Path


def counters(path):
    report = json.loads(Path(path).read_text())
    trace = report['davidson_profile']
    fields, operations = trace['fields'], trace['operations']
    result = []
    for rank in trace['per_rank']:
        if rank['dropped']:
            raise ValueError('trace dropped records')
        logical = actual = predicted = products = tiles = 0
        for row in rank['records']:
            d = dict(zip(fields, row))
            op = operations[int(d['operation'])]
            n, m, k, start = (int(d[x]) for x in ('n', 'm', 'k', 'detail'))
            if op == 'upper_border_product':
                if m - start != k:
                    raise ValueError('invalid border shape')
                products += 1
                logical += 2 * n * m * k
                predicted += 2 * n * sum(min(a + 64, m) * (min(a + 64, m) - a)
                                        for a in range(start, m, 64))
            elif op == 'upper_border_tile':
                tiles += 1
                actual += 2 * n * m * k
        if not products or not logical:
            raise ValueError('unusable: no instrumented border work')
        if tiles and actual != predicted:
            raise ValueError('tile work disagrees with predicted shapes')
        result.append(dict(rank=rank['rank'], products=products, tiles=tiles,
                           full_rectangle_flop=logical,
                           predicted_tiled_flop=predicted,
                           executed_flop=actual if tiles else logical))
    return dict(path=str(path), ranks=result)


if __name__ == '__main__':
    print(json.dumps([counters(p) for p in sys.argv[1:]], indent=2))
