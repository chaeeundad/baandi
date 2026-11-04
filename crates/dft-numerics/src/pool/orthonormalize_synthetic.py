#!/usr/bin/env python3
"""Bounded six-step synthetic Gamma regression sweep, not convergence evidence.
Run at repo root: python3 -I crates/dft-numerics/src/pool/orthonormalize_synthetic.py
BINARY NEW_OUTPUT_DIR [REFERENCE_SUMMARY_JSON]. Outputs generated UPFs only in
NEW_OUTPUT_DIR, never the repository. Explicit Xx/He mismatch is intentional.
The synthetic fixtures have no transferability or atomic-reference claim.
"""
import importlib.util, pathlib, subprocess, json, hashlib, sys, time
s=importlib.util.spec_from_file_location('g','scripts/gen_synthetic_upf.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
base=pathlib.Path(sys.argv[2]).resolve();base.mkdir(exist_ok=False);binary=str(pathlib.Path(sys.argv[1]).resolve());summary=[]
for n in [4,8,20]:
 for skew in [False,True]:
  label=f'syn{n}-'+('skew' if skew else 'cubic');upf=base/(label+'.upf')
  upf.write_text(g.generate(2, (0,0,1) if skew else (0,), points=101))
  length={4:8.,8:10.,20:14.}[n]
  cell=[[length,0.,0.],[1.5 if skew else 0.,length,0.],[.5 if skew else 0.,.7 if skew else 0.,length]]
  for t in [1,2,4,8]:
   d=base/f'{label}-T{t}';d.mkdir(exist_ok=True)
   text=f'[system]\nname="synthetic"\nlattice_unit = "bohr"\nlattice_vectors = {cell}\n'
   for i in range(n):
    site=[(i%4)/4+.07,(i//4)/((n+3)//4)+.03,.35 if i%2==0 else .65]
    text+=f'[[system.atoms]]\nspecies="Xx"\nfractional={site}\n'
   text+=f'[pseudopotentials]\nXx={{path="{upf}",element="He",allow_element_mismatch=true}}\n[basis]\nkpoints=[[0.0,0.0,0.0]]\necutwfc_ry=20.0\n[xc]\nfunctional="lda-pz"\n[scf]\nthreads={t}\nmax_iterations=6\nnbands={n+4}\n'
   inp=d/'input.toml';inp.write_text(text)
   start=time.perf_counter();r=subprocess.run([binary,'run',str(inp),'--output-dir',str(d/'out')],capture_output=True,text=True)
   (d/'run.log').write_text(r.stdout+r.stderr)
   assert r.returncode in (0,3), (label,t,r.returncode,r.stderr[-400:])
   result=d/'out/results.json'
   if not result.exists(): raise RuntimeError((label,t,r.returncode,r.stderr[-400:]))
   q=json.loads(result.read_text());h=q['convergence']['history'];cs=[x['fixed_point_certificate'] for x in h if x.get('fixed_point_certificate')]
   row=dict(case=label,threads=t,exit_code=r.returncode,wall_s=time.perf_counter()-start,energy_ha=q['total_energy']['ha'],density_sha256=hashlib.sha256((d/'out/density.bin').read_bytes()).hexdigest(),davidson_iterations=sum(x['davidson_iterations'] for x in h),certificate_iterations=sum(x['davidson_iterations'] for x in cs),iterations=len(h),convergence=q['convergence'].get('converged'))
   summary.append(row);(base/'summary.json').write_text(json.dumps(summary,indent=2));print(row,flush=True)

assert len(summary) == 24
for case in {r['case'] for r in summary}:
    records = [r for r in summary if r['case'] == case]
    assert len({(r['energy_ha'],r['density_sha256'],r['davidson_iterations']) for r in records}) == 1
if len(sys.argv) > 3:
    reference = {(r['case'],r['threads']):r for r in json.loads(pathlib.Path(sys.argv[3]).read_text())}
    for row in summary:
        for key in ['energy_ha','density_sha256','davidson_iterations','certificate_iterations']:
            assert row[key] == reference[row['case'],row['threads']][key], (row['case'],row['threads'],key)
print('PASS bounded synthetic cross-thread identity; six SCF steps are not a convergence certificate.')
