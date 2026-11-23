#!/usr/bin/env python3
"""One archived input, one MPI8 job. Master-only dev08 runner; never runs QE.

AA is the noise floor. AB/BA are crossed pairs. Each invocation has a total
budget, defaults to 600 s; --job-name final-* permits at most 7200 s.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import signal
import time


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True,type=Path)
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--upf',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--job-name',required=True)
    p.add_argument('--sequence',choices=['AA','AB','BA'],required=True)
    p.add_argument('--budget',type=int,default=600)
    p.add_argument('--strict-result',required=True,type=Path)
    args=p.parse_args()
    cap=7200 if args.job_name.startswith('final-') else 600
    if not 30 < args.budget <= cap: p.error('budget exceeds job-class cap or is too small')
    args.out.mkdir(parents=True,exist_ok=False)
    ref=json.loads(args.strict_result.read_text())
    strict_history=ref.get('convergence',{}).get('history',[])
    certified=(ref.get('convergence',{}).get('converged') is True
        and ref.get('single_path',{}).get('effective_profile') == 'strict'
        and strict_history and strict_history[-1].get('fixed_point_certificate') is not None)
    if not certified:
        p.error('strict reference must have converged and certified; check input identity separately')
    ref_e=ref['total_energy']['ha']
    manifest=dict(job_name=args.job_name,sequence=args.sequence,budget_s=args.budget,
        binary_sha256=sha(args.binary),input_sha256=sha(args.input),upf_sha256=sha(args.upf),
        strict_result_sha256=sha(args.strict_result),reference_energy_ha=ref_e,runs=[])
    start=time.monotonic()
    for index,arm in enumerate(args.sequence):
        remain=args.budget-(time.monotonic()-start)-10
        if remain<=0: break
        env=dict(os.environ,OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
                 BAANDI_MIX_CELL_KERKER='off' if arm=='A' else 'on',
                 BAANDI_MIX_COULOMB_GUARD='off',BAANDI_STOP_PREDICTOR='off',BAANDI_REFERENCE_STRICT='off')
        out=args.out/f'{index}-{arm}'
        cmd=['taskset','-c','0-7','mpirun','-np','8','--bind-to','core','--map-by','core',
             str(args.binary.resolve()),'run',str(args.input.resolve()),'--output-dir',str(out.resolve())]
        before=time.monotonic()
        with (args.out/f'{index}-{arm}.stdout').open('w') as stdout, (args.out/f'{index}-{arm}.stderr').open('w') as stderr:
            child=subprocess.Popen(cmd,env=env,stdout=stdout,stderr=stderr,start_new_session=True)
            try: code=child.wait(timeout=remain)
            except subprocess.TimeoutExpired:
                # Own process group only; never kill by a process name.
                os.killpg(child.pid,signal.SIGTERM)
                try: child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid,signal.SIGKILL)
                    child.wait()
                code=124
        entry=dict(arm=arm,command=cmd,exit_code=code,wall_s=time.monotonic()-before,
                   environment={k:v for k,v in env.items() if k.startswith('BAANDI_') or k.endswith('NUM_THREADS')})
        result=out/'results.json'
        if result.exists():
            d=json.loads(result.read_text()); h=d.get('convergence',{}).get('history',[])
            energy=d.get('total_energy',{}).get('ha')
            verified=d.get('verified') is True
            entry.update(results_sha256=sha(result),fine_maps=len(h),energy_ha=energy,verified=verified,
                standard_error_trace_ha=[r.get('standard_error_ha') for r in h],
                density_trace=[r.get('density_residual') for r in h],
                predicted_ratio_trace=[(r.get('pulay') or {}).get('predicted_ratio') for r in h],
                resets=sum(bool(r.get('mixer_reset')) for r in h),
                accepted=code==0 and verified and energy is not None and abs(energy-ref_e)<=5e-9,
                density_sha256=sha(out/'density.bin') if (out/'density.bin').exists() else None)
        manifest['runs'].append(entry)
        (args.out/'job.json').write_text(json.dumps(manifest,indent=2)+'\n')
        if code!=0: break
    complete=len(manifest['runs'])==2 and all(r.get('accepted',False) for r in manifest['runs'])
    manifest['complete_and_correct']=complete
    (args.out/'job.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return 0 if complete else 3


if __name__=='__main__': raise SystemExit(main())
