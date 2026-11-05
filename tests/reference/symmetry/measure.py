#!/usr/bin/env python3
"""Mac observation, paired full/TR-only/crystal timings at identical discretization."""
import argparse, hashlib, json, os, platform, runpy, statistics, struct, subprocess, tempfile, time
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--binary',required=True);p.add_argument('--upf',required=True);p.add_argument('--output',required=True)
a=p.parse_args()
repo=Path(__file__).resolve().parents[3]
bench=runpy.run_path(str(repo/'scripts/bench_scf.py'))
base=bench['case_input'](str(Path(a.upf).resolve()),'k444s',1,[])
base=base.replace('ecutwfc_ry = 30.0','ecutwfc_ry = 30.0\nfft_grid=[32,32,32]')
base=base.replace('energy_tolerance_ha = 1e-8','energy_tolerance_ha = 1e-12\neigen_tolerance = 1e-11').replace('density_residual_tolerance = 1e-7','density_residual_tolerance = 1e-11')
inputs={name:base.replace('[xc]',f'ktime_reversal = {tr}\n[xc]')+f'\n[symmetry]\nmode="{mode}"\n' for name,tr,mode in [('full','false','off'),('time_reversal','true','off'),('crystal','true','auto')]}
record={'kind':'Mac observation, no dev08 acceptance','platform':platform.platform(),'load_before':os.getloadavg(),'binary_sha256':hashlib.sha256(Path(a.binary).read_bytes()).hexdigest(),'upf_sha256':hashlib.sha256(Path(a.upf).read_bytes()).hexdigest(),'inputs':inputs,'cases':{name:[] for name in inputs}}
densities={}
with tempfile.TemporaryDirectory() as temp:
 root=Path(temp)
 for repeat in range(3):
  for name in (list(inputs) if repeat%2==0 else list(reversed(inputs))):
   input_path=root/f'{name}.toml';input_path.write_text(inputs[name]);out=root/f'{name}-{repeat}'
   command=[str(Path(a.binary).resolve()),'run',str(input_path),'--output-dir',str(out)]
   start=time.perf_counter();run=subprocess.run(command,capture_output=True,text=True);wall=time.perf_counter()-start
   if run.returncode: raise RuntimeError(run.stdout+run.stderr)
   result=json.loads((out/'results.json').read_text());h=result['convergence']['history']
   record['cases'][name].append({'wall_s':wall,'engine_s':result['timings']['engine_total_s'],'energy_ha':result['total_energy']['ha'],'iterations':len(h),'operator_applications':sum(x['operator_applications'] for x in h),'certificate_applications':sum((x.get('fixed_point_certificate') or {}).get('operator_applications',0) for x in h),'irreducible_k_count':len(result['kpoints']['points']),'density_sha256':hashlib.sha256((out/'density.bin').read_bytes()).hexdigest()})
   raw=(out/'density.bin').read_bytes();densities[name]=struct.unpack('<'+'d'*(len(raw)//8),raw)
   volume=result['cell']['volume_bohr3'] if 'cell' in result else result['volume_bohr3']
   if name=='crystal':record['symmetry']=result['symmetry']
 for name,rows in record['cases'].items():
  assert len({r['density_sha256'] for r in rows})==1
  delta=abs(rows[0]['energy_ha']-record['cases']['full'][0]['energy_ha'])
  # Match engine density_distance: sqrt(integral(delta n)^2) / N.
  distance=(sum((x-y)**2 for x,y in zip(densities[name],densities['full']))*volume/len(densities[name]))**0.5/8.0
  record.setdefault('comparison',{})[name]={'delta_energy_ha':delta,'density_distance_bohr_minus_3_over_2':distance,'wall_median_s':statistics.median(r['wall_s'] for r in rows)}
  assert delta<=1e-10 and distance<=1e-9
record['load_after']=os.getloadavg()
Path(a.output).write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record['comparison'],indent=2))
