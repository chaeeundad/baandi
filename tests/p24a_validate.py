#!/usr/bin/env python3
"""P24A actual SCF, deterministic-state, GEMM and final dense-oracle evidence.
No timing acceptance, no pseudopotential contents are persisted.
"""
import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import tempfile

ROOT = Path(__file__).resolve().parents[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m
r = load('p15a', ROOT/'tests/p15a_replay.py')
s = load('scope', ROOT/'tests/reference/spectral/scope.py')
scf = r.load_module('bench_scf')

def inputs(case, work):
    if case.startswith('suite-'):
        suite=r.load_module('bench_suite')
        fixture=suite.load_case(ROOT/'benchmarks/cases'/f'{case[6:]}.toml')
        upf=work/'synthetic.UPF';upf.write_text(suite.synthetic(fixture))
        return suite.case_input(fixture,upf,1,1), {'Si':r.sha(upf)}
    if case in ('gamma','k222','k222s'):
        upf=ROOT.parent/'pseudo/Si.pz-vbc.UPF'
        return scf.case_input(str(upf),case,threads=1), {'Si':r.sha(upf)}
    if case.startswith(('si-','mgo-')):
        text=(ROOT/'tests/reference/materials'/f'{case}.toml').read_text()
        hashes={}
        for element in sorted(set(re.findall(r'species = "([A-Za-z]+)"',text))):
            path=Path(os.environ['DFT_NC_DIR'])/f'{element}.upf'
            hashes[element]=r.sha(path)
            assert hashes[element] in text
            text=text.replace(f'"{element}.upf"',json.dumps(str(path)))
        return re.sub(r'threads = \d+','threads = 1',text), hashes
    text,upf=r.inputs(case,work,ROOT.parent/'pseudo/Si.pz-vbc.UPF')
    return text, {'Si':r.sha(upf)}

def scientific(report,density):
    if not report: return None
    return dict(energy=report['total_energy'],density=density,
                eigenvalues=report['eigenvalues'],
                history=r.scientific(report['convergence']),
                locking=report.get('locking'),
                kpoints=report.get('kpoints'))

def gemms(stderr):
    rows=[]
    for line in stderr.splitlines():
        if line.startswith('P24A_GEMM,'):
            _,solve,step,site,m,n,k,calls,macs=line.split(',')
            rows.append(dict(solve=int(solve),step=int(step),site=site,m=int(m),n=int(n),k=int(k),calls=int(calls),multiply_adds=int(macs)))
    totals={}
    steps={}
    for row in rows:
        site=row['site']; totals[site]=totals.get(site,0)+row['multiply_adds']
        key=(row['solve'],row['step']);steps[key]=steps.get(key,0)+row['multiply_adds']
    return dict(total_multiply_adds=sum(totals.values()),sites=totals,
                steps=[dict(solve=k[0],step=k[1],multiply_adds=v) for k,v in sorted(steps.items())]),rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--cases',default='gamma,k222s,k222,k444s,si8,syn20-skew,syn4-elongated,syn4-slab,syn16-gamma,syn16-k211')
    ap.add_argument('--threads',default='1,2,4,8')
    ap.add_argument('--dense',action='store_true')
    ap.add_argument('--locking',choices=['off','auto'],default='auto')
    ap.add_argument('--guard',choices=['off','final'],default='off')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    a.binary=a.binary.resolve();a.baseline=a.baseline.resolve()
    for case in a.cases.split(','):
        with tempfile.TemporaryDirectory(prefix='p24a-') as td:
            work=Path(td);text,hashes=inputs(case,work)
            inp=work/'input.toml';inp.write_text(text)
            os.environ['P15C_DENSE_IDS']=''
            bp,br,bd,bc=r.run(a.baseline,inp,work/'baseline')
            base_work,base_shapes=gemms(bp.stderr)
            row=dict(case=case,locking=a.locking,input=text,upf_sha256=hashes,
                binary_sha256=r.sha(a.binary),baseline_binary_sha256=r.sha(a.baseline),
                baseline_exit=bp.returncode,baseline_command=bc,baseline_report=br,
                baseline_density=bd,baseline_work=base_work,runs=[])
            reference=None
            for threads in map(int,a.threads.split(',')):
                candidate=re.sub(r'threads\s*=\s*\d+',f'threads = {threads}',text).replace('[scf]',f'[scf]\nlocking="{a.locking}"\nspectral_guard="{a.guard}"')
                inp.write_text(candidate)
                proc,report,density,command=r.run(a.binary,inp,work/f't{threads}')
                (a.output/f'{case}-t{threads}-stderr.txt.gz').write_bytes(gzip.compress(proc.stderr.encode(),mtime=0))
                state=scientific(report,density)
                # Diagnostic solve ordering can change under k parallelism.
                if state: state['history'].pop('davidson_work',None)
                stats,shapes=gemms(proc.stderr)
                if reference is None:reference=state
                run=dict(threads=threads,exit_code=proc.returncode,command=command,
                    report=report,density_sha256=density,work=stats,
                    thread_identity=(state==reference) if report else None,
                    energy_delta_ha=(report['total_energy']['ha']-br['total_energy']['ha']) if report and br else None,
                    stderr_tail='\n'.join(x for x in proc.stderr.splitlines() if not x.startswith(('P14','P15','P24'))) [-3000:])
                if not report and a.dense and 'P15A_SOLVE' in proc.stderr:
                    solves,_=r.parse(proc.stderr)
                    run['failed_solve_spectra']=[]
                    for sid,solve in enumerate(solves):
                        dense=solve.get('dense_eigenvalues_ha')
                        steps=solve[0]['steps']
                        if dense and steps and br:
                            run['failed_solve_spectra'].append(dict(solve=sid,step=steps[-1]['step'],**s.compare(steps[-1]['eigenvalues_ha'],dense,br['occupations'])))
                    (a.output/f'{case}-failed-dense.json.gz').write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
                row['runs'].append(run)
                (a.output/f'{case}-t{threads}-gemm.json.gz').write_bytes(gzip.compress(json.dumps(shapes,separators=(',',':')).encode(),mtime=0))
                if threads==1 and a.dense and report:
                    nk=report.get('kpoints',{}).get('count',1);ks=sorted(set([0,nk//2,nk-1]));sid=0;ids=[];locations={}
                    for h in report['convergence']['history']:
                        final=h['iteration']==report['convergence']['iterations']
                        if final:
                            for k in ks:ids.append(sid+k);locations[sid+k]=dict(phase='ordinary',k=k)
                        sid+=nk
                        if h.get('fixed_point_certificate'):
                            if final:
                                for k in ks:ids.append(sid+k);locations[sid+k]=dict(phase='certificate',k=k)
                            sid+=nk
                    os.environ['P15C_DENSE_IDS']=','.join(map(str,ids))
                    dp,dr,dd,dc=r.run(a.binary,inp,work/'dense')
                    solves,_=r.parse(dp.stderr);comparisons=[]
                    for i in ids:
                        dense=solves[i].get('dense_eigenvalues_ha')
                        if dense:
                            values=solves[i][0]['steps'][-1]['eigenvalues_ha']
                            comparisons.append(dict(solve=i,**locations[i],**s.compare(values,dense,report['occupations'])))
                    row['dense']=dict(exit_code=dp.returncode,command=dc,same_state=scientific(dr,dd)==state,comparisons=comparisons)
                    (a.output/f'{case}-dense.json.gz').write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
                    os.environ['P15C_DENSE_IDS']=''
                print(case,threads,proc.returncode,run['energy_delta_ha'],run['thread_identity'],flush=True)
                (a.output/f'{case}.json').write_text(json.dumps(row,indent=2)+'\n')
            (a.output/f'{case}-baseline-gemm.json.gz').write_bytes(gzip.compress(json.dumps(base_shapes,separators=(',',':')).encode(),mtime=0))
if __name__=='__main__':main()
