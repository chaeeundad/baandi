#!/usr/bin/env python3
"""GLUE Phase 0 tables from results.json; no cross-rank or cross-model sums.

New spans are inclusive rank-local totals and means per entered call. Existing
history supplies fine-map envelopes; phase_profile supplies coarse/verification
buckets and FFT density detail. No substage trace is inferred from totals.
"""
import argparse
import json
import math
from pathlib import Path

PREFIX = 'glue_lane.'
PHASES = ('potential_s', 'density_s', 'total_energy_s', 'mixing_s')


def number(value, name, integer=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError(f'{name}: expected finite nonnegative number')
    if integer and not isinstance(value, int):
        raise ValueError(f'{name}: expected integer')
    return value


def parse(path, correctness_only=False):
    report = json.loads(path.read_text())
    run = report.get('run', {})
    if run.get('completed') is not True:
        raise ValueError(f'{path}: incomplete/missing run report')
    rank, ranks = run.get('rank'), run.get('ranks')
    number(rank, 'rank', True)
    number(ranks, 'ranks', True)
    if not 0 <= rank < ranks:
        raise ValueError(f'{path}: invalid rank/ranks')
    lane = {k[len(PREFIX):]: number(v, k, True)
            for k, v in run.get('counters', {}).items() if k.startswith(PREFIX)}
    if not lane or 'serial.mix.calls' not in lane or 'grid.potential.calls' not in lane:
        raise ValueError(f'{path}: missing GLUE counters (wrong binary/schema)')
    switches = report.get('scf', {}).get('switches', {})
    mode = lambda key: switches.get('BAANDI_' + key, {}).get('value')
    if mode('GLUE_PROFILE') == 'off' or mode('PHASE_PROFILE') == 'off':
        if any(lane.values()):
            raise ValueError(f'{path}: counters populated with profiler disabled')
    if mode('STAGE_TIMERS') == 'off' and any(v for k, v in lane.items() if k.endswith('.ns')):
        raise ValueError(f'{path}: new timers populated in count-only mode')
    spans = []
    timed = any(v for k, v in lane.items() if k.endswith('.ns'))
    for key, calls in lane.items():
        if not key.endswith('.calls'):
            continue
        name = key[:-6]
        if name + '.ns' not in lane:
            raise ValueError(f'{path}: missing nanoseconds for {name}')
        ns = lane[name + '.ns']
        if ns and not calls:
            raise ValueError(f'{path}: nanoseconds without calls for {name}')
        seconds = ns / 1e9 if timed and not correctness_only else None
        spans.append(dict(span=name, calls=calls, seconds=seconds,
                          seconds_per_call=seconds/calls if seconds is not None and calls else None))
    for key in lane:
        if key.endswith('.ns') and key[:-3] + '.calls' not in lane:
            raise ValueError(f'{path}: timer without call counter {key}')
    work = {k: v for k, v in lane.items() if not k.endswith(('.calls', '.ns'))}
    maps = []
    for row in report.get('convergence', {}).get('history', []):
        for kind, obj in [('ordinary', row), ('strict_certificate', row.get('fixed_point_certificate'))]:
            if not obj:
                continue
            values = {k: number(obj['timings'][k], k) for k in PHASES}
            maps.append(dict(iteration=row['iteration'], kind=kind,
                             seconds=None if correctness_only else values))
    profile = report.get('phase_profile') or report.get('parallel', {}).get('phases') or {}
    envelope = []
    for b in profile.get('buckets', []):
        envelope.append(dict(bucket=b['bucket'], **{
            bound: {k: v for k, v in b.get(bound, {}).items()
                    if k in ('maps', 'potential', 'density', 'energy', 'mixing', 'wall')}
            for bound in ('min', 'max')}))
    payload = {k + '_bytes': 8*v for k, v in work.items() if k.endswith(('send_values', 'recv_values'))}
    payload['dist_dot_flops'] = 2*work.get('dist.dot_elements', 0)
    payload['dist_update_flops'] = 4*work.get('dist.update_elements', 0)
    payload['dist_update_logical_bytes'] = 32*work.get('dist.update_elements', 0)
    return dict(file=str(path), run_id=run.get('run_id'), rank=rank, ranks=ranks,
                source_commit=run.get('source_commit'), binary=run.get('binary'),
                state='timed' if timed else 'counts_only' if any(lane.values()) else 'disabled_or_not_reached',
                spans=spans, work=work, logical_work=payload, fine_maps=maps,
                phase_envelopes=None if correctness_only else envelope,
                fft_lane=None if correctness_only else profile.get('fft_lane'),
                phase_per_rank=None if correctness_only else profile.get('per_rank'),
                existing_owned_counters=report.get('parallel', {}).get('m2', {}).get('counters'),
                shapes={k: report.get(k) for k in ('basis', 'bands', 'threading')},
                energy_ha=report.get('total_energy', {}).get('ha'),
                converged=report.get('convergence', {}).get('converged'), verified=report.get('verified'))


def markdown(data):
    out = ['# GLUE lane profile', '', data['label'], '',
           'Local inclusive spans; never add parent/child times or rank maxima.',
           'Zero calls mean not reached. A dash means no usable timing.', '']
    fmt = lambda v: '-' if v is None else f'{v:.9g}'
    for r in data['runs']:
        out += [f"## {r['file']}", '', f"Rank {r['rank']}/{r['ranks']}; {r['state']}.", '',
                '| span | calls | seconds | seconds/call |', '| --- | ---: | ---: | ---: |']
        out += [f"| {s['span']} | {s['calls']} | {fmt(s['seconds'])} | {fmt(s['seconds_per_call'])} |" for s in r['spans']]
        out += ['', '| work / size | count |', '| --- | ---: |']
        out += [f'| {k} | {v} |' for k, v in r['work'].items()]
        out += ['', '| fine iteration | map kind | potential s | density s | energy s | mixing s |',
                '| ---: | --- | ---: | ---: | ---: | ---: |']
        for m in r['fine_maps']:
            out.append(f"| {m['iteration']} | {m['kind']} | " + ' | '.join(fmt((m['seconds'] or {}).get(k)) for k in PHASES) + ' |')
        out += ['', 'Coarse/fine/verification envelopes (already reduced over ranks):', '',
                '| bucket | bound | potential s | density s | energy s | mixing s |',
                '| --- | --- | ---: | ---: | ---: | ---: |']
        for b in r['phase_envelopes'] or []:
            for bound in ('min', 'max'):
                out.append(f"| {b['bucket']} | {bound} | " + ' | '.join(fmt(b[bound].get(k)) for k in ('potential','density','energy','mixing')) + ' |')
        out.append('')
    return '\n'.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('paths', type=Path, nargs='+')
    ap.add_argument('--json-out', type=Path, required=True)
    ap.add_argument('--markdown-out', type=Path, required=True)
    ap.add_argument('--correctness-only', action='store_true', help='suppress timing values for Mac smoke checks')
    args = ap.parse_args()
    files = sorted({p.resolve() for arg in args.paths for p in (arg.rglob('results.json') if arg.is_dir() else [arg])})
    if not files:
        ap.error('no results.json files')
    try:
        data = dict(schema='lane-glue-profile-1',
                    label='Mac correctness only; timings suppressed.' if args.correctness_only else 'dev08 diagnostics; no optimization/adoption claim.',
                    runs=[parse(p, args.correctness_only) for p in files])
    except (ValueError, KeyError, TypeError, OSError) as e:
        ap.error(str(e))
    args.json_out.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    args.markdown_out.write_text(markdown(data)+'\n')
    print(f'Parsed {len(files)} result files, preserving ranks and runs.')


if __name__ == '__main__':
    main()
