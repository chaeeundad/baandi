#!/usr/bin/env python3
"""Summarize stored P24A measurements; MACs are complex multiply-add counts."""
import gzip
import json
import math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tests/reference/perf/p24a'
def read(path):return json.loads(path.read_text())
def work(shapes):
    rows=json.loads(gzip.decompress(shapes.read_bytes()))
    steps={};maxstep={};sites={};widths={}
    for r in rows:
        key=(r['solve'],r['step'])
        steps[key]=steps.get(key,0)+r['multiply_adds']
        maxstep[r['solve']]=max(maxstep.get(r['solve'],0),r['step'])
        sites[r['site']]=sites.get(r['site'],0)+r['multiply_adds']
        if r['site']=='ritz_and_residual':widths[key]=r['n']
    tail={k:v for k,v in steps.items() if k[1]>0 and k[1]>=max(1,math.ceil(0.75*maxstep[k[0]]))}
    return dict(total_multiply_adds=sum(steps.values()),sites=sites,
        tail_definition='last quarter of numbered steps in each solve, including terminal rotation/certificate; setup step zero excluded',
        tail_steps=len(tail),tail_multiply_adds=sum(tail.values()),
        mean_tail_multiply_adds=sum(tail.values())/len(tail) if tail else None,
        mean_tail_ritz_width=sum(widths[k] for k in tail if k in widths)/sum(k in widths for k in tail) if any(k in widths for k in tail) else None)
def counters(report):
    if not report:return None
    h=report['convergence']['history']
    cert=[x['fixed_point_certificate'] for x in h if x.get('fixed_point_certificate')]
    return dict(iterations=report['convergence']['iterations'],
        steps=sum(x['davidson_iterations'] for x in h+cert),
        h_vectors=sum(x['operator_applications'] for x in h+cert))
def main():
    result=dict(default='off',reason='Generality energy gate and full-root dense gate failed; work targets also missed.',
        mac_units='complex multiply-adds (one half of historical P14B 2mnk)',cases=[],ablations=[],suite=[],materials=[])
    for path in sorted((OUT/'final-diagnostic').glob('*.json')):
        x=read(path);r=x['runs'][0];case=x['case']
        before=work(path.with_name(case+'-baseline-gemm.json.gz'))
        after=work(path.with_name(case+'-t1-gemm.json.gz'))
        result['cases'].append(dict(case=case,exit_code=r['exit_code'],energy_delta_ha=r['energy_delta_ha'],
            before_work=before,after_work=after,before_counters=counters(x['baseline_report']),after_counters=counters(r['report']),
            mac_reduction_fraction=1-after['total_multiply_adds']/before['total_multiply_adds'] if before['total_multiply_adds'] and r['report'] else None,
            dense=x.get('dense',{}).get('comparisons',[])))
    for directory in sorted((OUT/'ablations').iterdir()):
        if not directory.is_dir():continue
        for path in sorted(directory.glob('*.json')):
            x=read(path);r=x['runs'][0]
            result['ablations'].append(dict(variant=directory.name,case=x['case'],exit_code=r['exit_code'],energy_delta_ha=r['energy_delta_ha'],
                work=r['work'],counters=counters(r['report'])))
    b=read(OUT/'final-suite-off.json')['cases'];a=read(OUT/'final-suite-auto.json')['cases']
    for name in sorted(a):
        x=a[name];y=b[name]
        delta=x.get('total_energy_ha',0)-y.get('total_energy_ha',0)
        result['suite'].append(dict(case=name,before_status=y['status'],after_status=x['status'],
            before_iterations=y.get('iterations'),after_iterations=x.get('iterations'),energy_delta_ha=delta,
            energy_gate_pass=(abs(delta)<=1e-12) if x['status']==y['status']=='pass' else None))
    for path in sorted((OUT/'materials').glob('*.json')):
        x=read(path);r=x['runs'][0];d=x['dense']['comparisons']
        result['materials'].append(dict(case=x['case'],energy_delta_ha=r['energy_delta_ha'],before_counters=counters(x['baseline_report']),
            after_counters=counters(r['report']),occupied_dense_max_ha=max(v['occupied_max_error_ha'] for v in d),
            all_dense_roots_pass=all(v['all_roots_within_cap'] for v in d)))
    result['determinism']=[]
    for directory in ['final-threads','material-threads','final-off']:
        for path in sorted((OUT/directory).glob('*.json')):
            x=read(path)
            result['determinism'].append(dict(group=directory,case=x['case'],threads=[r['threads'] for r in x['runs']],
                exits=[r['exit_code'] for r in x['runs']],state_identity=[r['thread_identity'] for r in x['runs']],
                off_baseline_density_identity=[r['density_sha256']==x['baseline_density'] for r in x['runs']] if directory=='final-off' else None))
    (OUT/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    lines=['# P24A measured work and validation','', 'Default remains **off**. MACs count complex multiply-adds, not real FLOPs. Timings are concurrent Mac observations, not acceptance.', '',
        '| Case | MAC before | MAC after | Reduction | Tail mean before | Tail mean after | H columns before/after | SCF before/after |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for x in result['cases']:
        b=x['before_work'];a=x['after_work'];bc=x['before_counters'];ac=x['after_counters']
        if not ac:continue
        lines.append(f"| {x['case']} | {b['total_multiply_adds']} | {a['total_multiply_adds']} | {100*x['mac_reduction_fraction']:.2f}% | {b['mean_tail_multiply_adds']:.1f} | {a['mean_tail_multiply_adds']:.1f} | {bc['h_vectors']}/{ac['h_vectors']} | {bc['iterations']}/{ac['iterations']} |")
    lines+=['','Tail = final quarter of numbered steps per solve, including final rotation/certificate. Raw per-step/site shapes are the adjacent `*-gemm.json.gz` artifacts.','',
        '## Generality suite','', '| Case | Status before/after | SCF before/after | Energy change (Ha) |', '|---|---|---:|---:|']
    for x in result['suite']:
        lines.append(f"| {x['case']} | {x['before_status']}/{x['after_status']} | {x['before_iterations']}/{x['after_iterations']} | {x['energy_delta_ha']:.9e} |")
    (OUT/'tables.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':main()
