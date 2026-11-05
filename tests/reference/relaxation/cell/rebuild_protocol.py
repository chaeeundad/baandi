"""Run P13D material rebuild/polish oracles without changing frozen-only evidence.

A polish cycle is a NEW ionic optimizer run at the previous accepted geometry,
with the nominal cutoff selecting a new frozen basis. It is not --resume or an
exact electronic continuation. No changes to convergence tolerances are made.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent
PRESSURE_TOLERANCE = 1e-7  # Existing real-file stress tolerance, Ha/Bohr^3.


def input_text(template, nc_dir, every, frame=None):
    text = template
    for element in ("Si", "Mg", "O"):
        text = re.sub(r'^' + element + r'\s*=.*$',
                      element + ' = ' + json.dumps(str(nc_dir / (element + '.upf'))),
                      text, flags=re.MULTILINE)
    text = text.replace('[relax]\n', '[relax]\nrebuild_every = %d\n' % every)
    if frame is not None:
        columns = frame['cell_columns_bohr']
        rows = [[columns[i][j] for i in range(3)] for j in range(3)]
        text = re.sub(r'^lattice_vectors\s*=.*$', 'lattice_vectors = ' + json.dumps(rows),
                      text, flags=re.MULTILINE)
        fractional = iter(frame['fractional'])
        text = re.sub(r'^fractional\s*=.*$', lambda _: 'fractional = ' + json.dumps(next(fractional)),
                      text, flags=re.MULTILINE)
    return text


def run_stage(binary, work, name, text):
    input_path = work / (name + '.toml')
    input_path.write_text(text)
    print('START', name, flush=True)
    start = time.monotonic()
    completed = subprocess.run([str(binary), 'relax', str(input_path), '--output-dir', str(work / name)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.monotonic() - start
    (work / (name + '.stdout')).write_bytes(completed.stdout)
    (work / (name + '.stderr')).write_bytes(completed.stderr)
    result_path = work / name / 'results.json'
    result = json.loads(result_path.read_text()) if result_path.exists() else {}
    r = result.get('relaxation', {})
    row = dict(stage=name, elapsed_s=elapsed, exit_code=completed.returncode,
               status=r.get('status'), input=str(input_path), provenance=result.get('provenance'),
               failures=r.get('failures'))
    frozen, rebuilt = r.get('frozen_final'), r.get('rebuilt_final')
    if frozen is not None and rebuilt is not None:
        row.update(lattice_bohr=2*frozen['frame']['cell_columns_bohr'][0][1],
                   frozen_pressure_ha_bohr3=frozen['pressure_ha_bohr3'],
                   rebuilt_pressure_ha_bohr3=rebuilt['pressure_ha_bohr3'],
                   frozen_final=frozen, rebuilt_final=rebuilt,
                   rebuilt_pressure_check=r.get('rebuilt_pressure_check'),
                   steps=r['steps'],
                   basis_rebuild_count=sum(f['provenance']=='basis_rebuild' for f in r['frames']))
        print('RESULT %s a=%.12f Bohr Pf=%.9e Pr=%.9e Ha/Bohr^3 t=%.3f s exit=%d' %
              (name, row['lattice_bohr'], frozen['pressure_ha_bohr3'], rebuilt['pressure_ha_bohr3'],
               elapsed, completed.returncode), flush=True)
    else:
        print('FAILED', name, completed.returncode, completed.stderr.decode(errors='replace'), flush=True)
    return row


def compare(rows, bounds):
    pairs = []
    for material, cutoff in [('MgO', 50), ('Si', 30)]:
        for every in (3, 5):
            group = sorted([r for r in rows if r['material']==material and r['rebuild_every']==every],
                           key=lambda r: r['start_factor'])
            if len(group)!=2 or any('lattice_bohr' not in r['stages'][-1] for r in group):
                continue
            bound = next(r['estimated_combined_optimizer_tolerance_bohr'] for r in bounds
                         if r['material']==material and r['cutoff_ha']==cutoff)
            final = [r['stages'][-1] for r in group]
            initial = [r['stages'][0] for r in group]
            delta = final[1]['lattice_bohr'] - final[0]['lattice_bohr']
            pressure_ok = all(abs(r['rebuilt_pressure_ha_bohr3']) < PRESSURE_TOLERANCE for r in final)
            stages_ok = all(s['exit_code']==0 and s['status']=='converged' for r in group for s in r['stages'])
            pairs.append(dict(material=material, cutoff_ha=cutoff, rebuild_every=every,
                              before_polish_difference_bohr=initial[1]['lattice_bohr']-initial[0]['lattice_bohr'],
                              before_polish_rebuilt_pressures_ha_bohr3=[r['rebuilt_pressure_ha_bohr3'] for r in initial],
                              final_difference_bohr=delta, unchanged_estimated_tolerance_bohr=bound,
                              final_rebuilt_pressures_ha_bohr3=[r['rebuilt_pressure_ha_bohr3'] for r in final],
                              pressure_threshold_ha_bohr3=PRESSURE_TOLERANCE,
                              start_agreement=abs(delta)<=bound, pressure_agreement=pressure_ok,
                              passed=abs(delta)<=bound and pressure_ok and stages_ok))
    return pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--nc-dir', type=Path, default=os.environ.get('DFT_NC_DIR'))
    parser.add_argument('--work-dir', required=True, type=Path)
    parser.add_argument('--materials', default='MgO,Si')
    args=parser.parse_args()
    if args.nc_dir is None:
        parser.error('set DFT_NC_DIR or --nc-dir')
    binary=args.binary.resolve()
    nc_dir=Path(args.nc_dir).resolve()
    work=args.work_dir.resolve()
    work.mkdir(parents=True,exist_ok=True)
    materials=args.materials.split(',')
    if not materials or any(m not in ('MgO','Si') for m in materials):
        parser.error('materials must be MgO,Si or either one')
    baseline=json.loads((HERE/'real-followup/real-validation.json').read_text())
    rows=[]
    def save():
        pairs=compare(rows,baseline['start_comparisons'])
        data=dict(protocol='Periodic rebuild then 2-3 fresh-reference polish cycles; frozen-only records unchanged',
                  original_frozen_evidence='real-followup/real-validation.json',
                  pressure_tolerance_ha_bohr3=PRESSURE_TOLERANCE,
                  tolerance_note='Same previously reported combined lattice tolerance estimates, without widening or refitting.',
                  cases=rows, comparisons=pairs)
        (work/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
        return pairs
    for material in materials:
        cutoff=50 if material=='MgO' else 30
        for every in (3,5):
            for factor in (0.98,1.02):
                label='%s-%d-N%d-%s' % (material,cutoff,every,factor)
                template=(HERE/'real-followup'/('%s-%d-%s.toml'%(material,cutoff,factor))).read_text()
                row=dict(case=label,material=material,cutoff_ha=cutoff,rebuild_every=every,start_factor=factor,stages=[])
                rows.append(row)
                stage=run_stage(binary,work,label+'-main',input_text(template,nc_dir,every))
                row['stages'].append(stage)
                save()
                if stage['exit_code']!=0:
                    continue
                for cycle in range(1,4):
                    stage=run_stage(binary,work,label+'-polish%d'%cycle,
                                    input_text(template,nc_dir,0,stage['frozen_final']['frame']))
                    row['stages'].append(stage)
                    save()
                    if stage['exit_code']!=0:
                        break
                    if cycle>=2 and abs(stage['rebuilt_pressure_ha_bohr3'])<PRESSURE_TOLERANCE:
                        break
    pairs=save()
    for pair in pairs:
        print('PAIR',json.dumps(pair),flush=True)
    passed=len(pairs)==2*len(materials) and all(p['passed'] for p in pairs)
    print('PASS rebuild start agreement and rebuilt pressure' if passed else 'INCONCLUSIVE rebuild acceptance',flush=True)
    return 0 if passed else 1


if __name__=='__main__':
    raise SystemExit(main())
