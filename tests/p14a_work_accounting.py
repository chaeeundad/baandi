#!/usr/bin/env python3
"""P14A baseline accounting; no numerical policy or certificate changes.

Run with python3 -I tests/p14a_work_accounting.py --output /path/new.json.
Raw inputs/outputs and generated UPFs stay in a temporary directory. Only
summaries and input text (not UPF contents) are written to the evidence file.
Unavailable solver diagnostics are explicitly null, never inferred from totals.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / (name + '.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qe_rows():
    rows = []
    for path in sorted((ROOT / 'tests/reference/kpoints').glob('*/qe.out')):
        iteration = None
        for line in path.read_text().splitlines():
            match = re.search(r'iteration #\s*(\d+)', line)
            if match:
                iteration = int(match[1])
            match = re.search(r'ethr =\s*([\d.E+-]+),\s*avg # of iterations =\s*([\d.]+)', line)
            if match:
                rows.append(dict(source=str(path.relative_to(ROOT)), sha256=sha(path),
                                 iteration=iteration, ethr_ry=float(match[1]),
                                 ethr_ha=float(match[1]) / 2,
                                 average_iterations=float(match[2])))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/dft-cli')
    parser.add_argument('--upf', type=Path, default=ROOT.parent / 'pseudo/Si.pz-vbc.UPF')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cases', default='si8,k444s,syn20-skew,syn4-elongated')
    parser.add_argument('--extra-bands', type=int)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--label', default='default', help='record the experimental build policy')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output must be new')
    scf, suite = module('bench_scf'), module('bench_suite')
    data = dict(command=sys.argv, label=args.label, binary_sha256=sha(args.binary),
                commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                qe=qe_rows(), cases={}, limitations=[
                    'QE ethr is a Ritz-value change threshold, not a residual norm.',
                    'QE rows here use stored k222/k222s outputs; the dev08 relay is summarized separately.',
                    'Root diagnostics require P14A_WORK_ACCOUNTING=1; trace order is only grouped for one thread.'])
    with tempfile.TemporaryDirectory(prefix='p14a-') as td:
        work = Path(td)
        for name in args.cases.split(','):
            directory = work / name
            directory.mkdir()
            if name in ('si8', 'k444s'):
                upf = args.upf.resolve()
                text = scf.case_input(str(upf), name, threads=args.threads)
                occupied = 16 if name == 'si8' else 4
            else:
                if name in ('syn20-skew', 'syn20-cubic', 'syn20-cubic12'):
                    side = 12 if name == 'syn20-cubic12' else 8
                    lattice = ([[side, 0, 0], [0, side, 0], [0, 0, side]] if name.startswith('syn20-cubic')
                               else [[12, 0, 0], [1, 6, 0], [.5, .3, 6]])
                    atoms = [[(i+.13)/5, (j+.21)/2, (k+.32)/2]
                             for i in range(5) for j in range(2) for k in range(2)]
                elif name == 'syn4-elongated':
                    lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 16]]
                    atoms = [[.13, .21, (k+.32)/4] for k in range(4)]
                else:
                    parser.error('unknown case ' + name)
                case = dict(id=name, lattice=lattice, atoms=atoms, cutoff_ry=3,
                            valence=2, channels=[0, 0, 1], shape='soft', max_iterations=100)
                upf = directory / 'synthetic.UPF'
                upf.write_text(suite.synthetic(case))
                text = suite.case_input(case, upf, 1, args.threads)
                occupied = len(atoms)
            if args.extra_bands is not None:
                text += f'nbands = {occupied + args.extra_bands}\n'
            inp = directory / 'input.toml'
            inp.write_text(text)
            output = directory / 'out'
            command = [str(args.binary), 'run', str(inp), '--output-dir', str(output)]
            try:
                proc = subprocess.run(command, capture_output=True, text=True, timeout=600)
                entry = dict(exit_code=proc.returncode, input=text, upf_sha256=sha(upf),
                             stdout=proc.stdout, stderr=proc.stderr)
            except subprocess.TimeoutExpired:
                entry = dict(exit_code=None, input=text, upf_sha256=sha(upf), error='600 s timeout')
            report_path = output / 'results.json'
            if report_path.exists():
                report = json.loads(report_path.read_text())
                if (output / 'density.bin').exists():
                    entry['density_sha256'] = sha(output / 'density.bin')
                entry['measurements'] = suite.measurements(report)
                entry['davidson_work'] = report.get('davidson_work', [])
                bands = report['bands']
                entry['history'] = report['convergence']['history']
                for row in entry['history']:
                    row['h_applications_per_band'] = row['operator_applications'] / bands
                    row['h_applications_per_band_per_k'] = (
                        row['operator_applications'] / bands / len(row['operator_applications_per_k']))
                    phases = row['timings']
                    total = sum(phases.values())
                    row['phase_fractions'] = {k: v / total for k, v in phases.items()} if total else {}
                    row.update(restarts=None, early_converged_bands=None,
                               late_converged_bands=None, buffer_h_fraction=None)
            data['cases'][name] = entry
            args.output.write_text(json.dumps(data, indent=2) + '\n')
            print(name, entry['exit_code'], entry.get('measurements', {}).get('operator_applications'), flush=True)


if __name__ == '__main__':
    main()
