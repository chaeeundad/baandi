#!/usr/bin/env python3
"""P27D real SCF score. Explicit SKIP without external UPFs. Evidence stays outside target."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
QE = dict(zip(('si-20','si-30','si-40','mgo-40','mgo-50','mgo-60'),
              (-16.92463785,-16.92466016,-16.92467579,-152.10176985,-152.10194662,-152.10196248)))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=ROOT/'target/release/dft-cli')
    p.add_argument('--evidence', type=Path, required=True)
    a = p.parse_args()
    nc = os.environ.get('DFT_NC_DIR')
    if not nc or any(not (Path(nc)/f'{s}.upf').is_file() for s in ('Si','Mg','O')):
        print('SKIP: DFT_NC_DIR must contain Si.upf, Mg.upf, O.upf')
        return
    work = a.evidence.resolve()
    if 'target' in work.parts:
        raise ValueError('Evidence must be outside target')
    work.mkdir(parents=True, exist_ok=True)
    binary = a.binary.resolve()
    report = {'binary_sha256': sha(binary), 'qe_source': 'tests/reference/perf/dev08_materials_qe_relayed.md',
              'upf_sha256': {s: sha(Path(nc)/f'{s}.upf') for s in ('Si','Mg','O')}, 'rows': []}
    for case, qe in QE.items():
        source = ROOT/'tests/reference/materials'/f'{case}.toml'
        text = source.read_text().replace('[xc]', '[xc]\npw92_prefactor = "original"')
        for species in ('Si','Mg','O'):
            text = text.replace(f'"{species}.upf"', json.dumps(str((Path(nc)/f'{species}.upf').resolve())))
        for mesh in ('baandi','qe_like'):
            for rule in ('simpson_3_8','simpson_quadratic_end'):
                dest = work/f'{case}-{mesh}-{rule}'
                dest.mkdir(exist_ok=True)
                inp = dest/'input.toml'
                inp.write_text(text + f'\n[pseudo]\nlocal_mesh = "{mesh}"\nradial_even_rule = "{rule}"\n')
                command = [str(binary),'run',str(inp),'--output-dir',str(dest/'out')]
                run = subprocess.run(command, capture_output=True, text=True)
                (dest/'stdout.log').write_text(run.stdout)
                (dest/'stderr.log').write_text(run.stderr)
                if run.returncode:
                    raise RuntimeError(f'{dest.name}: exit {run.returncode}: {run.stderr[-1000:]}')
                r = json.loads((dest/'out/results.json').read_text())
                assert r['convergence']['converged']
                row = dict(case=case, local_mesh=mesh, radial_even_rule=rule,
                           total_energy_ha=r['total_energy']['ha'], total_energy_ry=2*r['total_energy']['ha'],
                           qe_ry=qe, residual_ry=2*r['total_energy']['ha']-qe,
                           iterations=r['convergence']['iterations'], exit_code=run.returncode,
                           input_sha256=sha(inp), density_sha256=sha(dest/'out/density.bin'), command=command)
                report['rows'].append(row)
                (work/'score.json').write_text(json.dumps(report,indent=2)+'\n')
                print(case, mesh, rule, f"{row['total_energy_ry']:.14f}", f"{row['residual_ry']:+.10e}", flush=True)
    # Explicit defaults versus omission on all six original-prefactor cases.
    for case in QE:
        dest = work/f'{case}-omitted'
        dest.mkdir(exist_ok=True)
        base = work/f'{case}-baandi-simpson_3_8'
        inp=dest/'input.toml'
        inp.write_text((base/'input.toml').read_text().split('\n[pseudo]')[0])
        run=subprocess.run([str(binary),'run',str(inp),'--output-dir',str(dest/'out')],capture_output=True,text=True)
        (dest/'stdout.log').write_text(run.stdout); (dest/'stderr.log').write_text(run.stderr)
        assert run.returncode == 0, run.stderr
        before=json.loads((base/'out/results.json').read_text())
        after=json.loads((dest/'out/results.json').read_text())
        for key in ('total_energy','energies','eigenvalues','occupations'):
            assert before.get(key)==after.get(key), (case,key)
        assert (base/'out/density.bin').read_bytes()==(dest/'out/density.bin').read_bytes()
    report['explicit_defaults_equal_omission_all_six']=True
    (work/'score.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__ == '__main__':
    main()
