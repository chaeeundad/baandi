"""Bitwise standard-profile thread checks, including failed trajectories."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import equal_accuracy as ea


def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items() if k!='timings'}
    if isinstance(v,list):return [clean(x) for x in v]
    return v


rows=[]
with tempfile.TemporaryDirectory(prefix='p17a-threads-') as td:
    td=Path(td)
    for name in ea.CASES:
        directory=td/name;directory.mkdir()
        base,hashes=ea.case_input(name,directory,ROOT.parent/'pseudo/Si.pz-vbc.UPF',1)
        for metric in ['l2','hartree']:
            records=[]
            for thread in [1,8]:
                inp=directory/f'{metric}-{thread}.toml'
                inp.write_text(ea.scf_settings(base,dict(convergence_profile='standard',mixing_metric=metric,threads=thread)))
                out=directory/f'{metric}-{thread}'
                command=[str(ROOT/'target/release/dft-cli'),'run',str(inp),'--output-dir',str(out)]
                p=subprocess.run(command,capture_output=True,text=True,timeout=1800)
                result=json.loads((out/'results.json').read_text()) if (out/'results.json').exists() else {}
                scientific=clean({k:result.get(k) for k in ['convergence','total_energy','eigenvalues','band_accuracy','kpoints','work','standard_threshold_ha','convergence_profile','certified']})
                density=ea.sha(out/'density.bin') if (out/'density.bin').exists() else None
                records.append(dict(threads=thread,command=command,exit_code=p.returncode,
                    scientific_sha256=hashlib.sha256(json.dumps(scientific,sort_keys=True).encode()).hexdigest(),density_sha256=density,
                    iterations=result.get('convergence',{}).get('iterations')))
            same=all(records[0][k]==records[1][k] for k in ['exit_code','scientific_sha256','density_sha256'])
            rows.append(dict(case=name,metric=metric,identical=same,runs=records,upf_sha256=hashes))
            print(name,metric,same,records[0]['exit_code'],flush=True)
            assert same,rows[-1]
(ROOT/'tests/reference/perf/p17a/thread-identity.json').write_text(json.dumps(dict(binary_sha256=ea.sha(ROOT/'target/release/dft-cli'),command='DFT_NC_DIR=/path/to/ex python3 -I tests/p17a_threads.py',rows=rows),indent=2)+'\n')
