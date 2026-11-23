#!/usr/bin/env python3
"""Recompute stage-frequency estimates and a deliberately optimistic MPI8 model.
Counts are local Mac replay observations, never Mac-to-dev08 timing ratios.
"""
import csv
import json
from pathlib import Path

lane = Path(__file__).resolve().parent
index = json.loads((lane / 'capture-index.json').read_text())
rows = []
combined = []
for case in index['cases']:
    records = []
    for sample in case['samples']:
        v = json.loads((lane / f"{case['case']}-{sample['phase']}.json").read_text())['result']
        records.append((sample, v))
        for c in v['candidates']:
            w = c['filter_work']
            charged = sum(i * n for i, n in enumerate(c['charged_h_calls_by_width']))
            rows.append(dict(case=case['case'], phase=sample['phase'], map=v['map'],
                degree=c['degree'] or 'auto', frequency=sample['weight_actual'], weight_11=sample['weight_11'],
                davidson_h_columns=v['baseline']['h_columns'], filter_h_columns=w['h_columns'],
                charged_h_columns=charged, sweeps=len(c['sweeps']), fallback=c['fallback'] or '',
                occupied_eigenvalue_error_ha=c['accuracy']['occupied_eigenvalue_max_error_ha'],
                band_energy_error_ha=c['accuracy']['frozen_band_energy_error_ha'],
                projector_spectral_error=c['accuracy']['projector_spectral_error'],
                density_max_error_e_per_bohr3=c['accuracy']['density_max_error_e_per_bohr3']))
    for i in range(5):
        a = dict(case=case['case'], degree=records[0][1]['candidates'][i]['degree'],
                 sampled_davidson_h_11=0., charged_h_11=0., filter_h_11=0., bounds_h_11=0.,
                 terminal_h_11=0., qr_flops_global_11=0., qr_bytes_global_11=0.,
                 gram_rotate_flops_global_11=0., recurrence_flops_global_11=0., copy_bytes_global_11=0.,
                 allocation_bytes_global_model_11=0., scalar_flops_global_model_11=0., eigen_sum_m3_11=0., eigen_calls_11=0.,
                 solver_reductions_11=0., fallback_maps_11=0., fallback_flops_global_11=0., fallback_eigen_m3_11=0.)
        for sample, v in records:
            f = sample['weight_11']; c = v['candidates'][i]; w = c['filter_work']; n = v['global_rows']; b = v['bands']
            a['sampled_davidson_h_11'] += f*v['baseline']['h_columns']
            a['charged_h_11'] += f*sum(i*n for i,n in enumerate(c['charged_h_calls_by_width']))
            for key, field in [('filter_h_11','h_columns'),('bounds_h_11','bound_h_columns'),
                ('terminal_h_11','terminal_h_columns'),('qr_flops_global_11','qr_flops_local'),
                ('recurrence_flops_global_11','recurrence_flops_local'),('copy_bytes_global_11','copy_bytes_local'),
                ('solver_reductions_11','reduction_calls')]: a[key] += f*w[field]
            a['gram_rotate_flops_global_11'] += f*(w['gram_flops_local']+w['rotate_flops_local'])
            a['qr_bytes_global_11'] += f*w['qr_calls']*(40*n*b*(b-1)+48*n*b)
            # Explicit zero-filled matrices/vectors, excluding opaque backend packing.
            a['allocation_bytes_global_model_11'] += f*8*(5*n*b+4*n+3*b*b+4*b+2+(n+255)//256*b*b+((b+3)//4)*4*b)
            a['scalar_flops_global_model_11'] += f*((10*(w['bound_h_columns']/2)+4)*n + 4*n*b*len(w['eigen_orders']))
            a['eigen_sum_m3_11'] += f*sum(m**3 for m in w['eigen_orders'])
            a['eigen_calls_11'] += f*len(w['eigen_orders'])
            if c['fallback']:
                a['fallback_maps_11'] += f
                counters = c['fallback_work']['counters']
                a['fallback_flops_global_11'] += f*counters['davidson_flop']
                a['fallback_eigen_m3_11'] += f*counters['reduced_eigen_full_m3']
        combined.append(a)
models=[]
h_price=2.782/6025
for a in (a for a in combined if a['case']=='si64'):
    b=132
    # Full 2.953s is an upper credit although only fine maps can change.
    # A fallback reruns Davidson: its dense algebra is NOT avoided.
    credit=2.953*(1-a['fallback_maps_11']/11)
    prices=[]
    for name,gflops,bw,h_factor,latency,max_order in [
        ('optimistic',108e9,12.97e9,0.74+0.26*64/108,1e-6,3*b),
        ('pessimistic',64e9,6.485e9,1.2,1e-5,b)]:
        pieces={
            'credit_s':credit,
            'incremental_h_s':a['charged_h_11']*h_price*h_factor-3504*h_price,
            'qr_s':(a['qr_flops_global_11']/8/gflops if name=='optimistic' else max(a['qr_flops_global_11']/8/gflops,a['qr_bytes_global_11']/8/bw)),
            'gram_rotate_s':a['gram_rotate_flops_global_11']/8/gflops,
            'recurrence_and_scalar_s':(a['recurrence_flops_global_11']+a['scalar_flops_global_model_11'])/8/gflops,
            'eigen_s':0.988*a['eigen_sum_m3_11']/(183*max_order**3),
            'copies_and_allocation_s':(a['copy_bytes_global_11']+a['allocation_bytes_global_model_11'])/8/bw,
            'solver_reductions_s':a['solver_reductions_11']*latency,
        }
        # Fallback H is charged in charged_h, fallback dense cost is retained
        # by removing its fraction from the credit; never grant it a saving.
        pieces['net_saving_s']=credit-sum(v for k,v in pieces.items() if k!='credit_s')
        prices.append(dict(scenario=name,**pieces))
    models.append(dict(degree=a['degree'],scenarios=prices))
output=dict(schema='filter-summary-1',combined_11_map=combined,si64_np8_model=models,
    h_price_s_per_column=h_price,decision='STOP integration; all sampled candidates lose even under optimistic assumptions',
    limitations=['Three-snapshot frequency extrapolation, not a full candidate SCF.',
                 'MPI8 projection uses Mac counts and relayed dev08 prices; no measured MPI8 filter wall.',
                 'Reduced EVD priced by cubic work between allowed orders 132 and 396, not by a DGEMM roof.',
                 'QR memory traffic is a streaming model; cache and latency require dev08 pricing.',
                 'Opaque backend allocations and production terminal refinement are not measured; zero credit for these nonnegative costs favors the candidate. Scalar residual and bound work is modeled from loop dimensions.',
                 '1/10 microsecond solver-reduction prices are sensitivity assumptions, not measured latencies.'])
(lane/'summary.json').write_text(json.dumps(output,indent=2)+'\n')
with (lane/'replay-summary.csv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
for m in models: print(m['degree'],[(s['scenario'],round(s['net_saving_s'],4))for s in m['scenarios']])
