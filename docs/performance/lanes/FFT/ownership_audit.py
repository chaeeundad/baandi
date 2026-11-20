#!/usr/bin/env python3
"""Arithmetic reconstruction, not profiling: cubic si64 full sphere, mirror sticks.
Assumes the relayed 23871-vector cubic Gamma basis and default mirror policy.
No dev08 runtime ownership is asserted; compare with returned trace shapes.
"""
import json
from collections import Counter


def audit():
    dims, ranks = 72, 8
    points = [(x % dims, y % dims, z % dims)
              for x in range(-18, 19) for y in range(-18, 19)
              for z in range(-18, 19) if x*x + y*y + z*z <= 320]
    assert len(points) == 23871
    lengths = Counter((z, y) for x, y, z in points)
    sticks = sorted(lengths)
    owners, loads = {}, [0] * ranks
    for s in sorted(sticks, key=lambda s: (-lengths[s], s)):
        if s in owners:
            continue
        rank = min(range(ranks), key=lambda r: (loads[r], r))
        for t in (s, ((-s[0]) % dims, (-s[1]) % dims)):
            if t in lengths and t not in owners:
                owners[t] = rank
                loads[rank] += lengths[t]
    counts = Counter(owners.values())
    assert sum(loads) == 23871 and sum(counts.values()) == len(sticks)
    return dict(assumption='cubic Gamma sphere |integer G|^2 <= 320; mirror sticks',
                plane_waves=len(points), sticks=len(sticks), active_z_planes=len({z for z,y in sticks}),
                ranks=[dict(rank=r, sticks=counts[r], coefficients=loads[r],
                            orc_rows=loads[r], transpose_planes=dims//ranks)
                       for r in range(ranks)])


if __name__ == '__main__':
    print(json.dumps(audit(), indent=2))
