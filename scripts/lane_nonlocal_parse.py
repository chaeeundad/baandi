#!/usr/bin/env python3
"""Parse NONLOCAL phase-profile counters. Mac mode deliberately excludes timings.

Rank extrema are separate observations: a maximum time and maximum work need
not belong to the same rank. Derived averages are bounds, never their ratio.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

STAGES = ('projection', 'reduce', 'coupling', 'back', 'add', 'audit')
SHAPES = ('rows_sum', 'projectors_sum', 'bands_sum', 'width1', 'width2',
          'width3_8', 'width9_32', 'width33_plus')


def number(values, name):
    value = values.get(name)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f'missing/invalid counter {name}: {value!r}')
    return value


def parse_result(path, machine):
    result = json.loads(Path(path).read_text())
    profile = result.get('phase_profile')
    if not isinstance(profile, dict) or not profile.get('buckets'):
        raise ValueError(f'{path}: PHASE_PROFILE data missing')
    rows = []
    for bucket in profile['buckets']:
        lo, hi = bucket['min'], bucket['max']
        for stage in STAGES:
            key = f'nonlocal.{stage}'
            calls_lo = number(lo, key + '.calls')
            calls_hi = number(hi, key + '.calls')
            row = dict(file=str(path), machine=machine, ranks=profile['ranks'],
                       bucket=bucket['bucket'], stage=stage,
                       status='observed' if calls_hi else 'unusable_zero_calls',
                       calls_min=calls_lo, calls_max=calls_hi,
                       maps_min=number(lo, 'maps'), maps_max=number(hi, 'maps'))
            for metric in ('flops', 'bytes'):
                row[metric + '_min'] = number(lo, key + '.' + metric)
                row[metric + '_max'] = number(hi, key + '.' + metric)
            if calls_lo > calls_hi:
                raise ValueError(f'{path}: reversed rank extrema for {key}')
            if stage == 'projection':
                for shape in SHAPES:
                    row[shape + '_min'] = number(lo, 'nonlocal.' + shape)
                    row[shape + '_max'] = number(hi, 'nonlocal.' + shape)
            if machine == 'dev08' and profile.get('stage_timers') == 'on':
                seconds_lo, seconds_hi = number(lo, key), number(hi, key)
                row.update(seconds_min=seconds_lo, seconds_max=seconds_hi)
                row['seconds_per_call_lower'] = seconds_lo / calls_hi if calls_hi else None
                row['seconds_per_call_upper'] = seconds_hi / calls_lo if calls_lo else None
                row['seconds_per_map_lower'] = seconds_lo / row['maps_max'] if row['maps_max'] else None
                row['seconds_per_map_upper'] = seconds_hi / row['maps_min'] if row['maps_min'] else None
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', nargs='+', type=Path)
    parser.add_argument('--machine', required=True, choices=('mac', 'dev08'))
    parser.add_argument('--format', choices=('json', 'csv'), default='csv')
    args = parser.parse_args()
    try:
        rows = [row for path in args.results for row in parse_result(path, args.machine)]
    except (ValueError, KeyError, OSError) as error:
        parser.exit(2, str(error) + '\n')
    if args.format == 'json':
        print(json.dumps(rows, indent=2, allow_nan=False))
    else:
        keys = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(sys.stdout, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
