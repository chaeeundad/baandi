#!/usr/bin/env python3
"""Regenerate the Mac evidence tables from collected results; never infer timing gains."""
import hashlib
import json
import math
from pathlib import Path
import statistics

LANE = Path(__file__).resolve().parent
ROOT = LANE.parents[3]
CASES = ['si-elongated','si-slab','si-long4','si8-gamma','si16-gamma','si32-gamma',
         'si-skew','syn-gaussian-fcc','syn-p3']
SEEDS = [601,602,603]


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k,v in value.items() if k not in ('timings','wall_s')}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def load(group, name, seed, arm):
    p = LANE / group / f'{name}-{seed}-{arm}.json'
    return json.loads(p.read_text()) if p.exists() else None


def counts(row):
    if row is None:
        return 'missing'
    return str(row['fine_maps']) + ('' if row['converged'] else 'F')


def reference_error(row, ref):
    if not row or not ref or not row['converged'] or not ref['converged'] or not ref['certified']:
        return None
    return abs(row['energy_ha'] - ref['energy_ha'])


def main():
    # Refresh complete work accounting from the CLI evidence retained under target.
    # Never rerun a calculation or invent a missing result. Standalone analysis
    # still works after cleanup: these fields are retained in the evidence JSON.
    for group in ['baseline','candidate','reference','combined','threads','strict-cell','final-off']:
        for p in (LANE / group).glob('*.json'):
            r = json.loads(p.read_text())
            raw = ROOT / r['command'][-1] / 'results.json'
            if raw.exists():
                d = json.loads(raw.read_text())
                r['certified'] = d.get('certified', bool(d.get('convergence',{}).get('converged')
                    and d.get('single_path',{}).get('effective_profile') == 'strict'
                    and r['history'] and r['history'][-1].get('fixed_point_certificate')))
                w = d.get('convergence',{}).get('work',{})
                ladder = w.get('cutoff_ladder') or {}
                r['coarse_start_h'] = ladder.get('coarse_start_h_vectors', 0)
                r['fine_start_h'] = d.get('atomic_start',{}).get('initial_h_applications', 0)
                r['verification_h'] = sum(x.get('h_vectors',0) for x in (d.get('verification') or {}).get('attempts',[]))
                r['coarse_h'] = ladder.get('coarse_h_vectors')
                if r['total_h'] is not None:
                    r['total_h'] = r['fine_h'] + r['certificate_h'] + (r['coarse_h'] or 0) + r['coarse_start_h'] + r['fine_start_h']
                r['effective_cell_filter'] = d.get('scf',{}).get('switches',{}).get('BAANDI_MIX_CELL_KERKER')
                r['raw_results_sha256'] = hashlib.sha256(raw.read_bytes()).hexdigest()
                # Human-readable tables accompany compact machine-readable traces.
                p.write_text(json.dumps(r, separators=(',',':')) + '\n')
    rows, identity = [], []
    out = ['# Mac work and correctness evidence', '',
        'Seeds are 601 / 602 / 603 in order; F means not converged. All H counts are vectors,',
        'including coarse/fine/certificate work and recorded atomic-start work.',
        'No wall-time inference. Blank references never pass an accuracy gate.', '',
        '| System | Baseline fine maps | Cell filter fine maps | Coulomb guard fine maps | Baseline H | Cell H | Max cell error vs strict (Ha) |',
        '|---|---|---|---|---|---|---|']
    for name in CASES:
        aa, cc, gg, errs = [], [], [], []
        for seed in SEEDS:
            a=load('baseline',name,seed,'baseline'); c=load('candidate',name,seed,'cell')
            g=load('candidate',name,seed,'coulomb'); b=load('candidate',name,seed,'baseline')
            ref=load('reference',name,seed,'strict')
            error=reference_error(c,ref)
            if error is not None: errs.append(error)
            aa.append(a); cc.append(c); gg.append(g)
            row=dict(case=name,seed=seed,baseline_maps=a['fine_maps'] if a else None,
                candidate_maps=c['fine_maps'] if c else None,candidate_converged=c['converged'] if c else False,
                baseline_h=a['total_h'] if a else None,candidate_h=c['total_h'] if c else None,
                cell_error_ha=error,baseline_error_ha=reference_error(a,ref),
                coulomb_error_ha=reference_error(g,ref),reference_converged=ref['converged'] if ref else False)
            rows.append(row)
            if a and b:
                identity.append(dict(case=name,seed=seed,energy=a['energy_ha']==b['energy_ha'],
                    density=a['density_sha256']==b['density_sha256'],
                    history=scientific(a['history'])==scientific(b['history']),
                    work=a['total_h']==b['total_h']))
        values=lambda rs,key: ' / '.join(str(r[key]) if r else 'missing' for r in rs)
        out.append(f'| {name} | {" / ".join(map(counts,aa))} | {" / ".join(map(counts,cc))} | {" / ".join(map(counts,gg))} | {values(aa,"total_h")} | {values(cc,"total_h")} | {max(errs):.6e} |' if errs else
                   f'| {name} | {" / ".join(map(counts,aa))} | {" / ".join(map(counts,cc))} | {" / ".join(map(counts,gg))} | {values(aa,"total_h")} | {values(cc,"total_h")} | unavailable |')
    out += ['', '## Baseline control signals', '',
        'The stretch column gives the first consecutive increase in standard error (inclusive map endpoints),',
        'not a claim of every later relapse. The predicted ratio minimum spans the full fine history.', '',
        '| System/seed | First growth stretch (factor) | Resets | Feasibility dampings | Min projected ratio | Largest L2 / best L2 | Frozen tolerance intervals |',
        '|---|---|---:|---:|---:|---:|---|']
    for name in CASES:
        for seed in SEEDS:
            a=load('baseline',name,seed,'baseline')
            if not a: continue
            s=a['summary']; stretches=s['stretches']; first=stretches[0] if stretches else None
            label=f'{first["start_map"]}-{first["end_map"]} ({first["growth"]:.3g}x)' if first else 'none'
            out.append(f'| {name}/{seed} | {label} | {s["resets"]} | {s["feasibility_damping"]} | {s["min_predicted_ratio"]} | {s["max_density_growth_from_best"]} | {s["frozen_tolerance_intervals"]} |')
    failed=[x for x in identity if not all(x[k] for k in ['energy','density','history','work'])]
    out += ['', '## Default-off identity', '',
        f'{len(identity)} paired pre-change / candidate-off runs inspected; {len(failed)} mismatches.',
        'Compared energy values, available density SHA-256, all scientific fine-history fields,',
        'and total recorded H. Timing/provenance fields are excluded. Failed runs have no',
        'certified density identity claim; their entire numerical histories still match.', '',
        '## Accuracy scope', '',
        f'{sum(r["reference_converged"] for r in rows)}/27 tight baseline references converged.',
        'A standard result is qualified here only when its unchanged verification passed AND',
        'an available tight reference differs by at most 5e-9 Ha. Failing stress cases remain',
        'unqualified; their last energies are not reference energies.']
    out += ['', '## Independent combination and strict/thread checks', '',
        '| System | C + Coulomb guard fine maps (601/602/603) | Max strict-C vs strict-A error (Ha) | C thread-1/thread-4 identities |',
        '|---|---|---:|---:|']
    additional=[]
    for name in CASES:
        combo=[]; strict_errors=[]; thread_checks=[]
        for seed in SEEDS:
            c=load('candidate',name,seed,'cell')
            t=load('threads',name,seed,'cell')
            both=load('combined',name,seed,'cell-coulomb')
            sc=load('combined',name,seed,'strict-cell')
            ref=load('reference',name,seed,'strict')
            error=reference_error(sc,ref)
            if error is not None: strict_errors.append(error)
            combo.append(both)
            same=None
            if c and t:
                same=(c['energy_ha']==t['energy_ha'] and c['density_sha256']==t['density_sha256']
                      and scientific(c['history'])==scientific(t['history']) and c['total_h']==t['total_h'])
                thread_checks.append(same)
            additional.append(dict(case=name,seed=seed,strict_error_ha=error,thread_identity=same,
                combined_maps=both['fine_maps'] if both else None,
                combined_error_ha=reference_error(both,ref)))
        e=f'{max(strict_errors):.6e}' if strict_errors else 'unavailable'
        out.append(f'| {name} | {" / ".join(map(counts,combo))} | {e} | {sum(thread_checks)}/{len(thread_checks)} |')
    final_identity=[]
    for name in CASES:
        a=load('baseline',name,601,'baseline'); b=load('final-off',name,601,'baseline')
        if a and b:
            final_identity.append(dict(case=name,identity=a['energy_ha']==b['energy_ha']
                and a['density_sha256']==b['density_sha256']
                and scientific(a['history'])==scientific(b['history']) and a['total_h']==b['total_h']))
    out += ['', f'Final-binary OFF confirmation: {sum(r["identity"] for r in final_identity)}/{len(final_identity)} systems (seed 601).',
        'Thread comparisons exclude timings and retain failed-run histories; unavailable delivered densities',
        'on failed calculations do not count as a certified density result.']
    (LANE/'evidence.md').write_text('\n'.join(out)+'\n')
    (LANE/'summary.json').write_text(json.dumps(dict(rows=rows,identity=identity,additional=additional,final_identity=final_identity),indent=2)+'\n')
    print(json.dumps(dict(pairs=len(identity),identity_failures=failed,
        references=sum(r['reference_converged'] for r in rows),
        max_cell_error_ha=max((r['cell_error_ha'] for r in rows if r['cell_error_ha'] is not None),default=None))))
    return bool(failed or any(not r['identity'] for r in final_identity) or any(r['thread_identity'] is False for r in additional))


if __name__ == '__main__':
    raise SystemExit(main())
