"""Final serial verification: no overlapping benchmark processes; no UPFs stored."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'tests/reference/perf/p17a'
p=argparse.ArgumentParser();p.add_argument('--work-dir',type=Path,required=True);a=p.parse_args()
a.work_dir.mkdir(parents=True,exist_ok=False)
commands=[]
env=dict(os.environ);env['PATH']='/opt/homebrew/opt/rustup/bin:'+env.get('PATH','')

def run(command, expected=0):
    log=a.work_dir/f'command-{len(commands)}.log'
    with log.open('w') as stream:
        result=subprocess.run(command,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT)
    commands.append(dict(command=command,exit_code=result.returncode,expected_exit_code=expected,log=str(log),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest()))
    (OUT/'final-verification.json').write_text(json.dumps(dict(commands=commands),indent=2)+'\n')
    print('exit',result.returncode,':',' '.join(command),flush=True)
    assert result.returncode==expected,log.read_text()[-5000:]

run([sys.executable,'-I','scripts/equal_accuracy.py','--output-dir',str(a.work_dir/'accuracy')],3)
for name in ['report.json','table.md']:
    shutil.copyfile(a.work_dir/'accuracy'/name,OUT/('final-'+name))
for threads in [1,8]:
    run([sys.executable,'-I','scripts/bench_scf.py','--binary','target/release/dft-cli',
         '--upf',str(ROOT.parent/'pseudo/Si.pz-vbc.UPF'),'--cases','gamma,k222s,k222,k444s,si8',
         '--repeats','1','--threads',str(threads),'--baseline',str(OUT/f'strict-t{threads}-before.json'),
         '--require-identical','--record',str(OUT/f'feature-strict-t{threads}.json')])
run([sys.executable,'-I','tests/p17a_threads.py'])
run([sys.executable,'-I','tests/p17a_cli_test.py'])
run([sys.executable,'-I','tests/p17a_equal_accuracy_test.py'])
run([sys.executable,'-I','scripts/tests/test_qe.py'])
run([sys.executable,'-I','scripts/tests/test_qe_real.py'])
run([sys.executable,'-I','scripts/tests/test_materials_qe.py'])
run(['scripts/check.sh'])
run([sys.executable,'-I','scripts/make_readme_figures.py'])
run([sys.executable,'-I','scripts/tests/test_readme_figures.py'])
source_files=[]
for directory in ['crates/dft-engine/src','crates/dft-solvers/src','crates/dft-cli/src']:
    source_files.extend(sorted((ROOT/directory).rglob('*.rs')))
source_files += [ROOT/'scripts/equal_accuracy.py',ROOT/'target/release/dft-cli']
record=json.loads((OUT/'final-verification.json').read_text())
record.update(status='passed_with_documented_reference_failures',baseline_commit='84f06a478549ce39151ba6095f5980bb652d8d26',
    sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in source_files},
    limits=['No remote QE runs; local timing samples are observations.','Accuracy harness exit 3 is expected: three stress fixtures have no tight reference; inspect every failed row.','Standard is not certified; no claim of spectral completeness when guard is off.'])
(OUT/'final-verification.json').write_text(json.dumps(record,indent=2)+'\n')
