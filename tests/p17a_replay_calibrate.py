"""Calibrate the proposed root policy on the existing frozen-H P15A traces."""
import gzip
import importlib.util
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('replay',ROOT/'tests/p15a_analyze.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
rows=[]
for name in ['si8','k444s','syn20-skew','syn4-elongated','syn4-slab','syn16-gamma','syn16-k211']:
 data=json.loads((r.OUT/(name+'.json')).read_text())
 data['solves']=json.loads(gzip.decompress((r.OUT/(name+'-traces.json.gz')).read_bytes()))
 schedule={v['iteration']:v for v in r.schedule(data['metrics'])}
 steps=vectors=native=missing=0
 for h,offset,group in r.groups(data):
  epsilon=schedule[h['iteration']]['ethr_ha']
  for solve in group:
   t=solve['0']; n=t['strict_bands']; count=len(t['steps'][0]['eigenvalues_ha'])
   # Occupied+first-empty required; first empty is loose but must converge.
   limits=[epsilon if j<n-1 else max(5*epsilon,5e-6) for j in range(count)]
   hit=r.stop(t,'ritz',epsilon,limits_override=limits)
   native+=len(t['steps'])
   if hit: steps+=hit['step'];vectors+=hit['h_vectors']
   else: missing+=1
 rows.append(dict(case=name,native_steps=native,ratcheted_ritz_prefix_steps=steps,h_vectors=vectors,censored=missing))
out=ROOT/'tests/reference/perf/p17a/replay-calibration.json'
out.write_text(json.dumps(dict(command='python3 -I tests/p17a_replay_calibrate.py',rows=rows,limits='Frozen trajectories only; not SCF speedups. No first-step change is invented. All roots must meet the rule in this calibration.'),indent=2)+'\n')
print(json.dumps(rows,indent=2))
