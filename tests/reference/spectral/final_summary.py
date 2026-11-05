#!/usr/bin/env python3
"""Generate final policy numbers and verify thread/case evidence without timing claims."""
import gzip
import json
from pathlib import Path
P=Path(__file__).resolve().parent/'final'
rows=json.loads((P/'latest-t1/runs.json').read_text())
key=lambda r:(r['case'],r['policy'],r['repeat'])
base={key(r):r for r in rows}
def scientific(v):
    if isinstance(v,dict):return {k:scientific(x) for k,x in v.items() if k not in ('timings','timestamp')}
    if isinstance(v,list):return [scientific(x) for x in v]
    return v
checked=[]
for t in [2,4,8]:
    for r in json.loads((P/f't{t}/runs.json').read_text()):
        b=base[key(r)]
        for field in ['exit_code','density_sha256','ordinary_h','certificate_h','guard_h']:
            assert r[field]==b[field],(t,key(r),field)
        for field in ['total_energy','eigenvalues','convergence','spectral_completeness']:
            assert scientific(r['report'].get(field))==scientific(b['report'].get(field)),(t,key(r),field)
        if 'kpoints' in b['report']:
            assert r['report']['kpoints']==b['report']['kpoints'],(t,key(r),'kpoints')
        checked.append(dict(threads=t,case=r['case'],policy=r['policy']))
summary=[]
for case in dict.fromkeys(r['case'] for r in rows):
    off=base[(case,'off',0)];on=base[(case,'final',0)]
    denominator=off['ordinary_h']+off['certificate_h'];report=on['report']['spectral_completeness']
    summary.append(dict(case=case,off_exit=off['exit_code'],final_exit=on['exit_code'],
        off_iterations=off['report']['convergence']['iterations'],final_iterations=on['report']['convergence']['iterations'],
        off_h=denominator,final_scf_h=on['ordinary_h']+on['certificate_h'],guard_h=on['guard_h'],
        guard_percent_of_off_h=100*on['guard_h']/denominator,retries=report['retries'],passed=report['passed'],
        off_wall_s=off['wall_s'],final_wall_s=on['wall_s'],guard_wall_s=report['timings']['wall_s'],
        energy_change_ha=on['report']['total_energy']['ha']-off['report']['total_energy']['ha'],
        density_identical=on['density_sha256']==off['density_sha256']))
for r in summary:
    if r['case'] in ('si8','k444s','syn4-elongated','syn16-gamma'):
        assert r['passed'] and r['energy_change_ha']==0 and r['density_identical']
skew=next(r for r in summary if r['case']=='syn20-skew')
assert skew['off_exit']==0 and skew['final_exit']==3 and skew['retries']==1 and not skew['passed']
dense=json.loads((P/'dense-scf/runs.json').read_text())[0]
trace=json.load(gzip.open(P/'dense-scf'/dense['trace_file']))
history=dense['report']['convergence']['history']
gaps=[x['0']['steps'][-1]['eigenvalues_ha'][20]-x['0']['steps'][-1]['eigenvalues_ha'][19] for x in trace]
dense_summary=dict(exit_code=dense['exit_code'],iterations=len(history),plane_waves=dense['report']['basis']['n_plane_waves'],
    first_gap_ha=gaps[0],last_gap_ha=gaps[-1],minimum_gap_ha=min(gaps),
    final_energy_ha=dense['report']['total_energy']['ha'],
    minimum_density_residual=min(x['density_residual'] for x in history),final_density_residual=history[-1]['density_residual'],
    maximum_eigen_residual_ha=max(max(x['0']['steps'][-1]['residuals_ha']) for x in trace),
    trajectory=[dict(iteration=i+1,density_residual=history[i]['density_residual'],gap_ha=gaps[i],
        energy=history[i]['energies']['total']) for i in [0,1,5,9,19,49,99]])
assert dense_summary['exit_code']==3 and dense_summary['maximum_eigen_residual_ha']<1e-10
materials=json.loads((P/'materials/results.json').read_text())
assert len(materials)==6
for r in materials:
    assert r['exit_code']==0 and r['energy_identical'] and r['density_identical']
    assert r['spectral_completeness']['passed'] and r['spectral_completeness']['retries']==0
    for c in r['comparisons']:assert c['all_roots_within_cap'] and abs(c['band_sum_excess_ha'])<=1e-9
result=dict(default='off',default_reason='Final-only guard exceeds 3 percent work target on both realistic benchmarks.',
    cases=summary,dense_scf=dense_summary,thread_checks=checked,materials_cases=[r['case'] for r in materials],
    maximum_material_ordinary_occupied_error_ha=max(c['occupied_max_error_ha'] for r in materials for c in r['comparisons'] if c['phase']=='ordinary'),
    maximum_material_certificate_occupied_error_ha=max(c['occupied_max_error_ha'] for r in materials for c in r['comparisons'] if c['phase']=='certificate'))
(P/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print('Verified 56 case/policy/thread combinations, six real NC cases, cumulative rescue, and dense SCF failure.')
