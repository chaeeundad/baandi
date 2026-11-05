#!/usr/bin/env python3
"""Recompute frozen-cap comparisons and a compact table from retained traces."""
import gzip
import importlib.util
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('scope',HERE/'scope.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
for p in HERE.rglob('*.json'):
    if p.stem == 'spin':continue
    trace=p.with_name(p.stem+'-traces.json.gz')
    if not trace.exists():continue
    d=json.loads(p.read_text());solves=json.load(gzip.open(trace))
    for row in d['comparisons']:
        solve=solves[row['solve']]
        # All scope cases use fixed spin-unpolarized occupations. The native
        # strict prefix is occupied bands plus the first empty state.
        count=solve['0']['strict_bands']-1
        occupations=[2.]*count+[0.]*(len(solve['0']['steps'][-1]['eigenvalues_ha'])-count)
        dense=solve['dense_eigenvalues_ha']
        row.update(s.compare(solve['0']['steps'][-1]['eigenvalues_ha'],dense,occupations))
        for variant,comparison in row.get('variants',{}).items():
            comparison.update(s.compare(solve[variant]['steps'][-1]['eigenvalues_ha'],dense,occupations))
    p.write_text(json.dumps(d,indent=2)+'\n')
rows=[]
for p in sorted(HERE.glob('*.json')):
    d=json.loads(p.read_text())
    if not isinstance(d,dict) or 'comparisons' not in d or p.stem=='spin':continue
    final=[x for x in d['comparisons'] if x['iteration']==d['iterations'] and x['phase']=='ordinary']
    if final:
        rows.append(dict(case=d['case'],converged=d['converged'],iterations=d['iterations'],sampled_k=[x['k'] for x in final],
            occupied_max_error_ha=max(x['occupied_max_error_ha'] for x in final),
            first_empty_max_error_ha=max(x['first_empty_error_ha'] for x in final),
            buffer_max_error_ha=max(x['buffer_max_error_ha'] for x in final),
            maximum_absolute_per_k_band_sum_excess_ha=max(abs(x['band_sum_excess_ha']) for x in final)))
(HERE/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
lines=['| Case | Converged | SCF steps | Sampled k indices | Occupied max (Ha) | First empty max (Ha) | Buffer max (Ha) | Max absolute per-k band-sum excess (Ha) |',
       '| --- | --- | ---: | --- | ---: | ---: | ---: | ---: |']
for r in rows:
    lines.append('| {case} | {converged} | {iterations} | {sampled_k} | {occupied_max_error_ha:.6e} | {first_empty_max_error_ha:.6e} | {buffer_max_error_ha:.6e} | {maximum_absolute_per_k_band_sum_excess_ha:.6e} |'.format(**r))
(HERE/'table.md').write_text('\n'.join(lines)+'\n')
