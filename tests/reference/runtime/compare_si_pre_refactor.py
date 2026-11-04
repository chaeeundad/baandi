"""Compare archived pre-P11A Si outputs against an existing candidate bench_scf report.

Run from the repository root; requires the read-only local Si UPF and offline vendor.
No pseudopotential is copied. See README.md for prerequisites.
"""
import os,json,tempfile,subprocess,shutil
from pathlib import Path
root=Path(__file__).resolve().parents[3]
os.chdir(root)
env=dict(os.environ,PATH='/opt/homebrew/opt/rustup/bin:'+os.environ['PATH'])
with tempfile.TemporaryDirectory(prefix='dft-pre-si-') as temp:
 t=Path(temp); archive=t/'source.tar'; source=t/'source';source.mkdir()
 with archive.open('wb') as f:
  subprocess.run(['git','archive','9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec'],stdout=f,check=True)
 subprocess.run(['tar','-xf',str(archive),'-C',str(source)],check=True)
 subprocess.run(['cargo','build','--offline','--release','-p','dft-cli','--target-dir',str(root/'target'),'--config','source.vendored-sources.directory='+json.dumps(str((root/'../rust-dft-vendor').resolve()))],cwd=source,env=env,check=True)
 binary=t/'pre-refactor-cli';shutil.copy2(root/'target/release/dft-cli',binary)
 args=['python3','-I','scripts/bench_scf.py','--binary',str(binary),'--upf','/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF','--cases','gamma,k222','--repeats','3','--threads','1','--baseline','tests/reference/runtime/si-baseline-comparison.json','--require-identical','--record','tests/reference/runtime/si-pre-refactor.json','--note','P11J independent pre-refactor 9495d9f on same shared macOS host; compare to current 43868d9']
 result=subprocess.run(args,env=env,capture_output=True,text=True)
 (root/'tests/reference/runtime/si-pre-refactor.log').write_text(result.stdout+result.stderr)
 import hashlib
 metadata={'reference_revision':'9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec','candidate_revision':'43868d9','reference_cli_sha256':hashlib.sha256(binary.read_bytes()).hexdigest(),'upf_sha256':hashlib.sha256(Path('/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF').read_bytes()).hexdigest(),'command':args,'exit_code':result.returncode}
 (root/'tests/reference/runtime/si-provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
 print(result.stdout)
 raise SystemExit(result.returncode)
