#!/usr/bin/env python3
"""Measure actual reported-state final/off policy, retaining every failure.
Wall samples on a shared Mac are observations; H vector counts are exact.
"""
import argparse
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import time
ROOT=Path(__file__).resolve().parents[3]
spec=importlib.util.spec_from_file_location('replay',ROOT/'tests/p15a_replay.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--binary',type=Path,required=True)
ap.add_argument('--upf',type=Path,default=ROOT.parent/'pseudo/Si.pz-vbc.UPF')
ap.add_argument('--cases',default='si8,k444s,syn20-skew,syn4-elongated,syn4-slab,syn16-gamma,syn16-k211')
ap.add_argument('--threads',type=int,default=1)
ap.add_argument('--repeats',type=int,default=1)
ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
rows=[]
for case in a.cases.split(','):
    with tempfile.TemporaryDirectory(prefix='p15c-final-') as temp:
        work=Path(temp)
        text,pseudo=r.inputs(case,work,a.upf.resolve())
        text=text.replace('threads = 1',f'threads = {a.threads}')
        for repeat in range(a.repeats):
            for policy in ['off','final']:
                inp=work/'input.toml';inp.write_text(text.replace('[scf]',f'[scf]\nspectral_guard = "{policy}"'))
                start=time.monotonic()
                proc,result,density,command=r.run(a.binary.resolve(),inp,work/f'{policy}-{repeat}')
                wall=time.monotonic()-start
                guard=result.get('spectral_completeness') if result else None
                if result:
                    history=result['convergence']['history']
                    ordinary=sum(h['operator_applications'] for h in history)
                    certificate=sum((h.get('fixed_point_certificate') or {}).get('operator_applications',0) for h in history)
                else:ordinary=certificate=None
                row=dict(case=case,policy=policy,repeat=repeat,threads=a.threads,command=command,
                    input=inp.read_text(),upf_sha256=r.sha(pseudo),binary_sha256=r.sha(a.binary),exit_code=proc.returncode,
                    wall_s=wall,report=result,density_sha256=density,ordinary_h=ordinary,certificate_h=certificate,
                    guard_h=guard['extra_operator_applications'] if guard else 0,
                    stderr='\n'.join(l for l in proc.stderr.splitlines() if not l.startswith('P15A_')))
                solves,_=r.parse(proc.stderr)
                if solves:
                    trace=f'{case}-{policy}-{repeat}-traces.json.gz'
                    (a.output/trace).write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
                    row['trace_file']=trace
                rows.append(row)
                (a.output/'runs.json').write_text(json.dumps(rows,indent=2)+'\n')
                print(case,policy,proc.returncode,result['convergence']['iterations'] if result else None,
                    'guard H',row['guard_h'],'retries',guard['retries'] if guard else None,flush=True)
                # A runner exit code cannot hide a failed certificate or missing output.
                assert proc.returncode in [0,3] and result is not None,proc.stderr[-1500:]
                if policy=='final':
                    assert guard is not None
                    assert guard['retries']==sum(e['retries'] for e in guard['events'])
                    assert guard['extra_operator_applications']==sum(e['extra_operator_applications'] for e in guard['events'])
                    assert guard['passed']==result['convergence']['converged']
                    if proc.returncode==0:
                        events=[e for e in guard['events'] if e['iteration']==result['convergence']['iterations']]
                        assert len(events)==result.get('kpoints',{}).get('count',1)
                        assert all(e['passed'] for e in events)
                else: assert guard is None
