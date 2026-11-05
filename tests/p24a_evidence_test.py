#!/usr/bin/env python3
"""Validate recorded gate outcomes, including expected scientific failures."""
import gzip
import importlib.util
import json
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tests/reference/perf/p24a'
def read(p):return json.loads(p.read_text())
spec=importlib.util.spec_from_file_location('validate',ROOT/'tests/p24a_validate.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
s=read(OUT/'summary.json')
assert s['default']=='off'
assert len(s['suite'])==19
assert sum(x['after_status']=='pass' for x in s['suite'])==16
assert [x['case'] for x in s['suite'] if x['energy_gate_pass'] is False]==['syn-s2']
assert len(s['materials'])==6
for x in s['materials']:
    assert abs(x['energy_delta_ha']) <= 1e-12
    assert x['all_dense_roots_pass']
assert [x['case'] for x in s['materials'] if x['before_counters']['iterations']!=x['after_counters']['iterations']]==['mgo-60']
assert len(s['determinism'])==21
failures=[]
for x in s['determinism']:
    assert x['threads']==[1,2,4,8]
    if x['state_identity'][0] is not None:assert all(x['state_identity'])
    else:
        assert x['exits']==[4]*4
        traces=[]
        for t in x['threads']:
            raw=gzip.open(OUT/x['group']/f"{x['case']}-t{t}-stderr.txt.gz",'rt').read()
            traces.append(re.sub(r'timings: PhaseTimings \{[^}]*\}', 'timings: omitted',raw))
        assert all(raw==traces[0] for raw in traces),x['case']
        failures.append(dict(case=x['case'],exit_code=4,scientific_error_history_identical=True))
    if x['group']=='final-off':assert all(x['off_baseline_density_identity'])
for p in (OUT/'reviewed-state').glob('*.json'):
    x=read(p);y=read(OUT/'final-threads'/p.name)
    a=x['runs'][0];b=y['runs'][0]
    assert a['exit_code']==b['exit_code']
    assert v.scientific(a['report'],a['density_sha256'])==v.scientific(b['report'],b['density_sha256']),x['case']
si=read(OUT/'final-diagnostic/si8.json')['dense']['comparisons']
assert si[0]['all_roots_within_cap'] is False
assert si[0]['unmatched_dense_indices_zero_based']==[17,19]
assert si[0]['occupied_max_error_ha']<1e-9
assert si[1]['all_roots_within_cap']
guard=read(OUT/'final-guard/si8.json')['runs'][0]['report']
assert guard['spectral_completeness']['passed']
assert guard['spectral_completeness']['retries']==1
assert guard['convergence']['iterations']==16
spin=read(OUT/'spin/results.json')
assert len(spin['runs'])==4 and all(x['thread_identity'] for x in spin['runs'])
assert abs(spin['energy_delta_ha'])<=1e-12
assert all(x['max_all_root_error_ha']<=1e-9 for x in spin['comparisons'])
(OUT/'evidence-validation.json').write_text(json.dumps(dict(
    status='passed', meaning='Evidence is internally consistent; adoption gates intentionally fail.',
    deterministic_failures=failures,thread_combinations=88),indent=2)+'\n')
print('P24A evidence validated: 88 thread combinations; scientific adoption failures retained.')
