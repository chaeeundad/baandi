#!/usr/bin/env python3
"""Reference mode SCFs and manual-option identity; no external files are copied."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
REFERENCE = dict(zip(('si-20', 'si-30', 'si-40', 'mgo-40', 'mgo-50', 'mgo-60'),
                    (-16.92463785, -16.92466016, -16.92467579, -152.10176985, -152.10194662, -152.10196248)))
EXPECTED_HASHES = dict(Si='39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282',
                      Mg='15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14',
                      O='2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=ROOT/'target/release/dft-cli')
    p.add_argument('--evidence', type=Path, required=True)
    a = p.parse_args()
    nc = Path(os.environ.get('DFT_NC_DIR', ROOT.parent/'pseudo/pseudodojo/ex'))
    if any(not (nc/f'{s}.upf').is_file() for s in EXPECTED_HASHES):
        print('SKIP: external PseudoDojo Si.upf, Mg.upf, O.upf missing; set DFT_NC_DIR')
        return
    hashes = {s: sha(nc/f'{s}.upf') for s in EXPECTED_HASHES}
    assert hashes == EXPECTED_HASHES, 'present pseudopotentials do not match the evidence identities'
    work = a.evidence.resolve()
    if 'target' in work.parts:
        raise ValueError('Evidence must be outside target')
    work.mkdir(parents=True, exist_ok=True)
    binary = a.binary.resolve()
    report = dict(binary_sha256=sha(binary), upf_sha256=hashes, rows=[])
    for case, reference in REFERENCE.items():
        text = (ROOT/'tests/reference/materials'/f'{case}.toml').read_text()
        for s in EXPECTED_HASHES:
            text = text.replace(f'"{s}.upf"', json.dumps(str((nc/f'{s}.upf').resolve())))
        outputs = {}
        for variant in (['mode', 'manual'] if case in ('si-20', 'mgo-40') else ['mode']):
            dest = work/f'{case}-{variant}'
            dest.mkdir(exist_ok=True)
            inp = dest/'input.toml'
            content = text + '\n[compat]\nmode="reference"\n' if variant == 'mode' else text.replace('[xc]', '[xc]\npw92_prefactor="original"') + '\n[pseudo]\nlocal_mesh="qe_like"\nradial_even_rule="simpson_quadratic_end"\n[symmetry]\nmode="off"\ndensity_group="mesh"\n'
            inp.write_text(content)
            command = [str(binary), 'run', str(inp), '--output-dir', str(dest/'out')]
            run = subprocess.run(command, capture_output=True, text=True)
            (dest/'stdout.log').write_text(run.stdout)
            (dest/'stderr.log').write_text(run.stderr)
            (dest/'command.json').write_text(json.dumps(dict(command=command, exit_code=run.returncode), indent=2))
            assert run.returncode == 0, (case, variant, run.stderr)
            r = json.loads((dest/'out/results.json').read_text())
            assert r['convergence']['converged']
            outputs[variant] = dict(total_energy_ha=r['total_energy']['ha'],
                                    iterations=r['convergence']['iterations'],
                                    density_sha256=sha(dest/'out/density.bin'),
                                    operator_applications=sum(h['operator_applications'] for h in r['convergence']['history']))
            if variant == 'mode':
                assert r['compat_mode'] == 'reference'
                assert r['compat']['pw92_prefactor'] == 'original'
                assert r['compat']['local_mesh'] == 'qe_like'
                assert r['compat']['radial_even_rule'] == 'simpson_quadratic_end'
                assert r['compat']['density_group'] == 'mesh'
        if 'manual' in outputs:
            assert outputs['mode'] == outputs['manual'], (case, outputs)
        row = dict(case=case, **outputs['mode'], total_energy_ry=2*outputs['mode']['total_energy_ha'],
                   reference_ry=reference, residual_ry=2*outputs['mode']['total_energy_ha']-reference,
                   manual_identical=True if 'manual' in outputs else None)
        report['rows'].append(row)
        (work/'score.json').write_text(json.dumps(report, indent=2)+'\n')
        print(case, row, flush=True)


if __name__ == '__main__':
    main()
