#!/usr/bin/env python3
"""Frozen-H replay runner. Build with P15A_REPLAY=1; outputs contain no UPFs.
python3 -I tests/p15a_replay.py --binary target/release/dft-cli --baseline /tmp/p15a-before-cli
"""
import argparse
import hashlib
import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'tests/reference/perf/p15a'


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / (name + '.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse(stderr):
    solves, metrics, dense = {}, [], {}
    for line in stderr.splitlines():
        if line.startswith('P15A_DENSE,'):
            _, sid, values = line.split(',', 2)
            dense[int(sid)] = json.loads(values)
        elif line.startswith('P15A_SOLVE,'):
            _, sid, variant, tol, strict, buffer, success = line.split(',')
            solves.setdefault(int(sid), {})[int(variant)] = dict(
                tolerance_ha=float(tol), strict_bands=int(strict), buffer_tolerance_ha=float(buffer),
                success=success == 'true', steps=[])
        elif line.startswith('P15A_STEP,'):
            _, sid, variant, step, basis, vectors, blocks, restarts, arrays = line.split(',', 8)
            values, end = json.JSONDecoder().raw_decode(arrays)
            residuals = json.loads(arrays[end + 1:])
            trace = solves[int(sid)][int(variant)]
            previous = trace['steps'][-1]['eigenvalues_ha'] if trace['steps'] else None
            changes = [abs(a-b) for a, b in zip(values, previous)] if previous else None
            unconverged = sum(r > (trace['tolerance_ha'] if j < trace['strict_bands']
                                  else trace['buffer_tolerance_ha']) for j, r in enumerate(residuals))
            trace['steps'].append(dict(ritz_changes_ha=changes, unconverged_residual=unconverged,
                step=int(step), basis=int(basis), h_vectors=int(vectors), h_blocks=int(blocks),
                restarts=int(restarts), eigenvalues_ha=values, residuals_ha=residuals))
        elif line.startswith('P15A_METRIC,'):
            _, iteration, dr2, electrons = line.split(',')
            metrics.append(dict(iteration=int(iteration), dr2_ry=float(dr2), electrons=float(electrons)))
    for sid, values in dense.items():
        solves[sid]['dense_eigenvalues_ha'] = values
    return [solves[i] for i in range(len(solves))], metrics


def inputs(name, directory, upf):
    scf, suite = load_module('bench_scf'), load_module('bench_suite')
    if name in ('si8', 'k444s'):
        return scf.case_input(str(upf), name, threads=1), upf
    if name == 'syn20-skew':
        lattice = [[12, 0, 0], [1, 6, 0], [.5, .3, 6]]
        atoms = [[(i+.13)/5, (j+.21)/2, (k+.32)/2]
                 for i in range(5) for j in range(2) for k in range(2)]
    elif name == 'syn4-elongated':
        lattice = [[4, 0, 0], [0, 4, 0], [0, 0, 16]]
        atoms = [[.13, .21, (k+.32)/4] for k in range(4)]
    elif name == 'syn4-slab':
        lattice = [[8, 0, 0], [0, 8, 0], [0, 0, 20]]
        atoms = [[(i+.13)/2, (j+.21)/2, .45] for i in range(2) for j in range(2)]
    else:
        lattice = [[12, 0, 0], [0, 8, 0], [0, 0, 8]]
        atoms = [[(i+.13)/4, (j+.21)/2, (k+.32)/2]
                 for i in range(4) for j in range(2) for k in range(2)]
    case = dict(id=name, lattice=lattice, atoms=atoms, cutoff_ry=3,
                valence=2, channels=[0, 0, 1], shape='soft', max_iterations=100)
    if name == 'syn16-k211':
        case['kgrid'] = [2, 1, 1]
    pseudo = directory / 'synthetic.UPF'
    pseudo.write_text(suite.synthetic(case))
    return suite.case_input(case, pseudo, 1, 1), pseudo


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k != 'timings'}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def run(binary, inp, out):
    command = [str(binary), 'run', str(inp), '--output-dir', str(out)]
    proc = subprocess.run(command, capture_output=True, text=True, timeout=1800)
    report = json.loads((out / 'results.json').read_text()) if (out / 'results.json').exists() else None
    density = sha(out / 'density.bin') if (out / 'density.bin').exists() else None
    return proc, report, density, command


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', type=Path, default=ROOT / 'target/release/dft-cli')
    ap.add_argument('--baseline', type=Path, default=Path('/tmp/p15a-before-cli'))
    ap.add_argument('--upf', type=Path, default=ROOT.parent / 'pseudo/Si.pz-vbc.UPF')
    ap.add_argument('--verify-repeat', action='store_true')
    ap.add_argument('--cases', default='si8,k444s,syn20-skew,syn4-elongated,syn4-slab,syn16-gamma,syn16-k211')
    args = ap.parse_args()
    for name in args.cases.split(','):
        with tempfile.TemporaryDirectory(prefix='p15a-') as td:
            directory = Path(td)
            text, pseudo = inputs(name, directory, args.upf.resolve())
            inp = directory / 'input.toml'
            inp.write_text(text)
            base, br, bd, bc = run(args.baseline.resolve(), inp, directory / 'before')
            proc, report, density, command = run(args.binary.resolve(), inp, directory / 'replay')
            solves, metrics = parse(proc.stderr)
            result = dict(case=name, commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                          input=text, upf_sha256=sha(pseudo), binary_sha256=sha(args.binary),
                          baseline_binary_sha256=sha(args.baseline), command=command, baseline_command=bc,
                          exit_code=proc.returncode, baseline_exit_code=base.returncode,
                          report=report, baseline_report=br, density_sha256=density, baseline_density_sha256=bd,
                          solves=solves, metrics=metrics,
                          other_stderr='\n'.join(s for s in proc.stderr.splitlines() if not s.startswith('P15A_')))
            if args.verify_repeat:
                repeated, rr, rd, _ = run(args.binary.resolve(), inp, directory / 'repeat')
                rs, rm = parse(repeated.stderr)
                result['determinism'] = dict(exit_equal=repeated.returncode == proc.returncode,
                    trace_identical=rs == solves, metrics_identical=rm == metrics,
                    density_identical=rd == density,
                    energy_identical=rr['total_energy'] == report['total_energy'])
                assert all(result['determinism'].values()), result['determinism']
            if report and br:
                result['equivalence'] = dict(scientific_history_identical=scientific(report['convergence']) == scientific(br['convergence']),
                                            eigenvalues_identical=report['eigenvalues'] == br['eigenvalues'],
                                            energy_delta_ha=report['total_energy']['ha']-br['total_energy']['ha'],
                                            density_identical=density == bd,
                                            history_length_equal=len(report['convergence']['history']) == len(br['convergence']['history']))
            result['trace_file'] = name + '-traces.json.gz'
            (OUT / result['trace_file']).write_bytes(gzip.compress(
                json.dumps(result.pop('solves'), separators=(',', ':'), allow_nan=False).encode(), mtime=0))
            (OUT / (name + '.json')).write_text(json.dumps(result, separators=(',', ':'), allow_nan=False) + '\n')
            print(name, proc.returncode, len(solves), result.get('equivalence'), flush=True)
            if not solves:
                raise RuntimeError('measurement trace missing: ' + proc.stderr[-1000:])


if __name__ == '__main__':
    main()
