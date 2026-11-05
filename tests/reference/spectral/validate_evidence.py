#!/usr/bin/env python3
"""Check factual scope claims against stored measured records, including negatives."""
import json
from pathlib import Path
P=Path(__file__).resolve().parent
rows=json.loads((P/'summary.json').read_text())
for row in rows:
    case=row['case']; d=json.loads((P/(case+'.json')).read_text())
    assert d['energy_identical'] and d['density_identical'] and d['history_identical'],case
    assert d['eigenvalue_cap_ha']==d['band_sum_cap_ha']==1e-9
    if not case.startswith('syn'):
        assert row['converged'],case
        for k in ['occupied_max_error_ha','first_empty_max_error_ha','buffer_max_error_ha','maximum_absolute_per_k_band_sum_excess_ha']:
            assert row[k] <= 1e-9,(case,k,row[k])
    if case=='syn20-skew':
        assert row['converged'] and row['occupied_max_error_ha']>1
        assert row['maximum_absolute_per_k_band_sum_excess_ha']>20
    if case in ('syn4-slab','syn16-k211'):
        assert not row['converged']
base=json.loads((P/'baseline-t1.json').read_text())['cases']
for t in [1,2,4,8]:
    current=json.loads((P/f'default-t{t}.json').read_text())['cases']
    for case,b in base.items():
        for key in ['total_energy_ha','density_sha256','iterations','operator_applications','certificate_applications']:
            assert current[case][key]==b[key],(t,case,key)
spin=json.loads((P/'spin.json').read_text())
assert spin['exact_p13b_fixture'] and spin['energy_identical'] and spin['density_identical']
for row in spin['final_four']:
    assert row['all_roots_within_cap'] and abs(row['band_sum_excess_ha']) <= 1e-9
print(f'{len(rows)} scope rows and 20 default thread/case combinations verified; negative outcomes retained.')
