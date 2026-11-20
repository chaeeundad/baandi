#!/usr/bin/env python3
"""Paired external-wall decision from fast_screen JSON; never measures time."""
import argparse
import json
import math
import statistics
from pathlib import Path


def differences(path):
    r = json.loads(path.read_text())['records']
    assert len(r['A']) == len(r['B']) >= 3
    for arm in r.values():
        assert all(x['exit'] == 0 and x['converged'] and x['verified'] for x in arm)
        assert len({x['density_sha256'] for x in arm}) == 1, 'non-reproducible arm'
    return [b['wall_s'] - a['wall_s'] for a, b in zip(r['A'], r['B'])]


def summarize(aa, ab):
    floor = max(abs(x) for x in aa)
    mean = statistics.mean(ab)
    se = statistics.stdev(ab) / math.sqrt(len(ab))
    z = abs(mean) / se if se else None
    significant = abs(mean) > 0 if se == 0 else z >= 3
    return dict(aa_differences_s=aa, ab_differences_s=ab, floor_s=floor,
                mean_difference_s=mean, standard_error_s=se, abs_mean_over_se=z,
                zero_se_nonzero_mean=se == 0 and mean != 0,
                passes_timing=mean < 0 and significant and -mean > 3 * floor,
                note='Timing screen only: counter, accuracy and generality gates are separate.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('aa', type=Path)
    p.add_argument('ab', type=Path)
    a = p.parse_args()
    print(json.dumps(summarize(differences(a.aa), differences(a.ab)), indent=2, allow_nan=False))
