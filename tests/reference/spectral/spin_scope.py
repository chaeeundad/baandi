#!/usr/bin/env python3
"""Exact P13B one-electron spin fixture through the existing Rust test generator."""
import importlib.util
import json
import gzip
import subprocess
from pathlib import Path
import sys
import tempfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('scope',HERE/'scope.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
r=s.replay
with tempfile.TemporaryDirectory(prefix='p15c-spin-') as td:
    work=Path(td)
    pseudo=work/'syn.upf'
    generator=Path('/tmp/p15c-spin-fixture')
    pseudo.write_bytes(subprocess.check_output([str(generator)]))
    text=f'''[system]
name='synthetic spin'
lattice_unit='bohr'
lattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
[[system.atoms]]
species='H'
fractional=[0.0,0.0,0.0]
[pseudopotentials]
H='{pseudo}'
[basis]
ecutwfc_ry=5.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional='lda-pw'
[scf]
nbands=4
eigen_tolerance=1e-10
density_residual_tolerance=1e-10
energy_tolerance_ha=1e-11
[spin]
mode='collinear'
total_moment=1.0
axis='z'
initial_moment_seed=0.5
'''
    inp=work/'input.toml';inp.write_text(text)
    bp,br,bd,bc=r.run(Path('/tmp/p15c-baseline-cli'),inp,work/'baseline')
    binary=Path(sys.argv[1])
    proc,report,density,command=r.run(binary,inp,work/'oracle')
    solves,_=r.parse(proc.stderr)
    assert report and solves,proc.stderr
    rows=[]
    for i,solve in enumerate(solves):
        dense=solve.get('dense_eigenvalues_ha')
        if dense:
            values=solve[0]['steps'][-1]['eigenvalues_ha']
            rows.append(dict(solve=i,spin_channel=i%2,**s.compare(values,dense,[1,0,0,0] if i%2==0 else [0,0,0,0])))
    result=dict(exact_p13b_fixture=True,generator_sha256=r.sha(generator),input=text,upf_sha256=r.sha(pseudo),binary_sha256=r.sha(binary),command=command,baseline_command=bc,
        exit_code=proc.returncode,baseline_exit_code=bp.returncode,converged=report['convergence']['converged'],
        iterations=report['convergence']['iterations'],energy_identical=report['total_energy']==br['total_energy'],density_identical=density==bd,
        comparisons=rows,final_four=rows[-4:])
    (HERE/'spin-traces.json.gz').write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
    (HERE/'spin.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['final_four'],indent=2))
