#!/usr/bin/env python3
"""Analyze external-wall paired records from fast_screen.py; never run DFT."""
import argparse
import json
import math
import statistics


def summarize(data):
    a, b = data['records']['A'], data['records']['B']
    assert len(a) == len(b) and len(a) >= 3
    assert all(r['exit'] == 0 and r['converged'] and r['verified'] for r in a + b)
    assert len({r['density_sha256'] for r in a + b}) == 1
    assert len({r['energy_ha'] for r in a + b}) == 1
    delta = [x['wall_s'] - y['wall_s'] for x, y in zip(a, b)]
    mean = statistics.mean(delta)
    se = statistics.stdev(delta) / math.sqrt(len(delta))
    z = abs(mean) / se if se else None
    return dict(pairs=len(delta), paired_savings_s=delta, mean_saving_s=mean,
                se_s=se, abs_mean_over_se=z,
                resolved_positive=mean > 0 and se > 0 and z >= 3)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--aa', required=True)
    p.add_argument('--ab')
    a = p.parse_args()
    with open(a.aa) as f:
        aa = summarize(json.load(f))
    ab = None
    if a.ab:
        with open(a.ab) as f:
            ab = summarize(json.load(f))
    # A resolved A/A difference invalidates the window; do not subtract drift.
    aa_resolved = aa['se_s'] == 0 or aa['abs_mean_over_se'] >= 3
    print(json.dumps(dict(aa=aa, ab=ab, aa_window_suspect=aa_resolved,
                         candidate=bool(ab and ab['resolved_positive'] and not aa_resolved)), indent=2))


if __name__ == '__main__':
    main()
