#!/usr/bin/env python3
"""Summarize crossed dev08 pairs; statistical screen, not correctness approval."""
import json
import math
from pathlib import Path
import statistics
import sys


def summarize(job, mode):
    job = Path(job)
    def sample(stem):
        if mode == 'solve':
            docs = [json.loads(line.split(' ', 1)[1]) for line in
                    Path(str(stem) + '.stdout').read_text().splitlines()
                    if line.startswith('DAVBENCH_JSON ')]
            if len(docs) != 1 or len(docs[0]['kinds']) != 1:
                raise ValueError('expected one warm solve record')
            k = docs[0]['kinds'][0]
            if not k['repeats_identical'] or not k['counters_deterministic']:
                raise ValueError('non-deterministic solve')
            value = k['solve_ms_median'] / 1000
        else:
            r = json.loads((stem / 'results.json').read_text())
            if not r['run']['completed']:
                raise ValueError('incomplete run')
            value = float(Path(str(stem) + '.external-s').read_text().strip().splitlines()[-1])
        if not math.isfinite(value) or value <= 0:
            raise ValueError('invalid timing')
        return value
    phases = {}
    for phase in ('aa', 'ab'):
        pairs = []
        for i in range(1, 4):
            a, b = (sample(job / f'{phase}-{i}-{arm}') for arm in ('a', 'b'))
            pairs.append(dict(a_s=a, b_s=b, relative_difference=(b - a) / a))
        d = [p['relative_difference'] for p in pairs]
        mean, se = statistics.mean(d), statistics.stdev(d) / math.sqrt(len(d))
        phases[phase] = dict(pairs=pairs, mean=mean, se=se,
                             abs_mean_over_se=abs(mean) / se if se else None,
                             min=min(d), max=max(d))
    floor = max(abs(p['relative_difference']) for p in phases['aa']['pairs'])
    effect = phases['ab']
    resolved = abs(effect['mean']) >= 3 * effect['se'] and effect['mean'] != 0
    return dict(job=str(job), metric='warm solve s' if mode == 'solve' else 'external s',
                phases=phases, aa_floor=floor,
                statistical_gain=effect['mean'] < 0 and resolved,
                exceeds_three_times_aa_floor=abs(effect['mean']) > 3 * floor,
                warning='Statistics only: verify counters and science before a verdict.')


if __name__ == '__main__':
    print(json.dumps(summarize(sys.argv[1], sys.argv[2]), indent=2))
