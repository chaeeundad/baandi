#!/usr/bin/env python3
"""Reuse the exact P15C/P13B spin fixture; no pseudopotential is persisted."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tests/reference/perf/p24a/spin'
OUT.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('r',ROOT/'tests/p15a_replay.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
reference=json.loads((ROOT/'tests/reference/spectral/spin.json').read_text())
lib=max((ROOT/'target/release/deps').glob('libdft_physics-*.rlib'),key=lambda p:p.stat().st_mtime)
with tempfile.TemporaryDirectory(prefix='p24a-spin-') as td:
    work=Path(td);generator=work/'generator'
    command=['rustc','--edition=2024',str(ROOT/'tests/reference/spectral/spin_fixture.rs'),'-L',f'dependency={ROOT}/target/release/deps','--extern',f'dft_physics={lib}','-o',str(generator)]
    subprocess.run(command,check=True)
    pseudo=work/'syn.upf';pseudo.write_bytes(subprocess.check_output([str(generator)]))
    assert r.sha(pseudo)==reference['upf_sha256']
    text=re.sub("H='[^']+'",f"H='{pseudo}'",reference['input'])
    rows=[];expected=None;inp=work/'input.toml'
    for t in [1,2,4,8]:
        inp.write_text(text.replace('[scf]',f'[scf]\nlocking="auto"\nspectral_guard="final"\nthreads={t}'))
        p,report,density,cmd=r.run(Path('/tmp/p24a-reviewed-cli'),inp,work/f't{t}')
        assert p.returncode==0,p.stderr
        assert report['spectral_completeness']['passed']
        state=[report['total_energy'],report['band_data'],report['density'],r.scientific(report['spectral_completeness'])]
        if expected is None:expected=state
        assert state==expected
        rows.append(dict(threads=t,exit_code=p.returncode,command=cmd,report=report,thread_identity=True))
    inp.write_text(text)
    p,baseline,bd,bc=r.run(Path('/tmp/p24a-baseline-cli'),inp,work/'baseline')
    assert p.returncode==0,p.stderr
    inp.write_text(text.replace('[scf]','[scf]\nlocking="auto"\nthreads=1'))
    os.environ['P15C_DENSE_IDS']=''
    p,report,dd,dc=r.run(Path('/tmp/p24a-final-diagnostic-cli'),inp,work/'oracle')
    assert p.returncode==0,p.stderr
    solves,_=r.parse(p.stderr);comparisons=[]
    for sid in range(len(solves)-4,len(solves)):
        solve=solves[sid];dense=solve['dense_eigenvalues_ha'];values=solve[0]['steps'][-1]['eigenvalues_ha']
        error=max(abs(x-y) for x,y in zip(values,dense));assert error<=1e-9
        comparisons.append(dict(solve=sid,max_all_root_error_ha=error))
    delta=report['total_energy']['ha']-baseline['total_energy']['ha'];assert abs(delta)<=1e-12
    result=dict(input=text,upf_sha256=r.sha(pseudo),generator_command=command,
        baseline_command=bc,oracle_command=dc,binary_sha256=r.sha(Path('/tmp/p24a-reviewed-cli')),
        oracle_sha256=r.sha(Path('/tmp/p24a-final-diagnostic-cli')),runs=rows,
        baseline_report=baseline,oracle_report=report,energy_delta_ha=delta,comparisons=comparisons)
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/'traces.json.gz').write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
    print('spin: T=1,2,4,8 identical; final four dense spectra pass; energy delta Ha',delta)
