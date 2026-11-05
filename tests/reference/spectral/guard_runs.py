#!/usr/bin/env python3
"""Run guarded synthetic trajectories; retain failures, never relabel exhaustion."""
import importlib.util
import json
from pathlib import Path
import tempfile
import time
import argparse
ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('replay', ROOT/'tests/p15a_replay.py')
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
ap=argparse.ArgumentParser()
ap.add_argument("--binary",type=Path,default=Path("/tmp/p15c-guard-cli"))
ap.add_argument("--output",type=Path,default=Path(__file__).parent/"guard-synthetic.json")
args=ap.parse_args()
rows = []
for name in ['syn20-skew','syn4-elongated','syn4-slab','syn16-gamma','syn16-k211']:
    with tempfile.TemporaryDirectory(prefix='p15c-guard-') as temp:
        work = Path(temp)
        text,pseudo = r.inputs(name,work,ROOT.parent/'pseudo/Si.pz-vbc.UPF')
        inp = work/'input.toml'
        inp.write_text(text)
        start = time.monotonic()
        proc, report, density, command = r.run(args.binary,inp,work/'out')
        row = dict(case=name,command=command,input=text,upf_sha256=r.sha(pseudo),
                   binary_sha256=r.sha(args.binary),exit_code=proc.returncode,
                   wall_s=time.monotonic()-start,report=report,density_sha256=density,stderr=proc.stderr)
        rows.append(row)
        args.output.write_text(json.dumps(rows,indent=2)+'\n')
        print(name,proc.returncode,report['convergence']['iterations'] if report else None,proc.stderr[-300:],flush=True)
