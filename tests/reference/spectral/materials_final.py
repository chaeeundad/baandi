#!/usr/bin/env python3
"""External P13D fixtures, final/off identity, same-H full-space final samples.
Requires a P15C_ORACLE=1 release binary and DFT_NC_DIR. Does not modify fixtures.
"""
import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
ROOT=Path(__file__).resolve().parents[3]
spec=importlib.util.spec_from_file_location('scope',Path(__file__).with_name('scope.py'))
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
r=s.replay
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--binary',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--cases',default='si-20,si-30,si-40,mgo-40,mgo-50,mgo-60')
a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
# Freeze the executable while other Cargo gates may rebuild target/release.
binary_temp=tempfile.TemporaryDirectory(prefix='p15c-oracle-bin-')
frozen=Path(binary_temp.name)/'dft-cli'
shutil.copy2(a.binary,frozen)
a.binary=frozen
nc=Path(os.environ['DFT_NC_DIR']).resolve()
rows=[]
for case in a.cases.split(','):
    with tempfile.TemporaryDirectory(prefix='p15c-material-final-') as td:
        work=Path(td)
        text=(ROOT/'tests/reference/materials'/f'{case}.toml').read_text().replace('threads = 4','threads = 1')
        hashes={}
        for e in sorted(set(re.findall(r'species = "([A-Za-z]+)"',text))):
            path=nc/f'{e}.upf';hashes[e]=r.sha(path);assert hashes[e] in text
            text=text.replace(f'"{e}.upf"',json.dumps(str(path)))
        inp=work/'input.toml';inp.write_text(text)
        os.environ['P15C_DENSE_IDS']=''
        bp,br,bd,bc=r.run(a.binary.resolve(),inp,work/'off')
        assert bp.returncode==0 and br is not None,bp.stderr[-1000:]
        nk=br.get('kpoints',{}).get('count',1)
        ks=sorted(set([0,nk//2,nk-1]));ids=[];locations={};sid=0
        for h in br['convergence']['history']:
            if h['iteration']==br['convergence']['iterations']:
                for k in ks:ids.append(sid+k);locations[sid+k]=dict(phase='ordinary',k=k)
            sid+=nk
            if h.get('fixed_point_certificate'):
                if h['iteration']==br['convergence']['iterations']:
                    for k in ks:ids.append(sid+k);locations[sid+k]=dict(phase='certificate',k=k)
                sid+=nk
        os.environ['P15C_DENSE_IDS']=','.join(map(str,ids))
        inp.write_text(text.replace('[scf]','[scf]\nspectral_guard="final"'))
        proc,result,density,command=r.run(a.binary.resolve(),inp,work/'final')
        assert proc.returncode==0 and result is not None,proc.stderr[-1000:]
        assert result['spectral_completeness']['passed']
        assert result['total_energy']==br['total_energy'] and density==bd
        solves,_=r.parse(proc.stderr)
        comparisons=[]
        for i in ids:
            dense=solves[i]['dense_eigenvalues_ha']
            values=solves[i][0]['steps'][-1]['eigenvalues_ha']
            comparison=s.compare(values,dense,result['occupations'])
            assert comparison['all_roots_within_cap'] and abs(comparison['band_sum_excess_ha'])<=1e-9,(case,i,comparison)
            comparisons.append(dict(solve=i,**locations[i],**comparison))
        row=dict(case=case,input=inp.read_text(),upf_sha256=hashes,command=command,baseline_command=bc,
            binary_sha256=r.sha(a.binary),dense_ids=ids,sampled_k=ks,exit_code=proc.returncode,
            iterations=result['convergence']['iterations'],energy_identical=True,density_identical=True,
            spectral_completeness=result['spectral_completeness'],comparisons=comparisons)
        rows.append(row)
        (a.output/f'{case}-traces.json.gz').write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
        (a.output/'results.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(case,'occupied',max(x['occupied_max_error_ha'] for x in comparisons),
            'all-root',max(max(x['occupied_max_error_ha'],x['first_empty_error_ha'],x['buffer_max_error_ha']) for x in comparisons),flush=True)
