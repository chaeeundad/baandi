"""Mac observation only: stress on/off identity and pass cost, no speedup claim."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import tempfile
import time

root = Path(__file__).resolve().parents[3]
bench = runpy.run_path(str(root / 'scripts/bench_scf.py'))
p = argparse.ArgumentParser()
p.add_argument('--binary', required=True)
p.add_argument('--upf', required=True)
p.add_argument('--output', required=True)
a = p.parse_args()
binary = str(Path(a.binary).resolve())
upf = str(Path(a.upf).resolve())
records = {}
for case in ['si8', 'k444s']:
    runs = []
    for enabled in [False, True]:
        with tempfile.TemporaryDirectory(prefix='p22-stress-cost-') as td:
            td = Path(td)
            text = bench['case_input'](upf, case, threads=1)
            text += '\n[output]\nstress=' + ('true' if enabled else 'false') + '\n'
            (td / 'input.toml').write_text(text)
            command = [binary, 'run', str(td / 'input.toml'), '--output-dir', str(td / 'out')]
            start = time.monotonic()
            result = subprocess.run(command, capture_output=True, text=True)
            wall = time.monotonic() - start
            if result.returncode:
                raise RuntimeError(result.stderr)
            report = json.loads((td / 'out/results.json').read_text())
            history = report['convergence']['history']
            runs.append({'enabled': enabled, 'exit_code': result.returncode, 'wall_seconds': wall,
                         'scf_seconds': report['timings']['engine_total_s'],
                         'energy_ha': report['total_energy']['ha'],
                         'density_sha256': hashlib.sha256((td / 'out/density.bin').read_bytes()).hexdigest(),
                         'iterations': len(history),
                         'operator_applications': sum(h['operator_applications'] for h in history),
                         'davidson_iterations': sum(h['davidson_iterations'] for h in history),
                         'stress': report.get('stress')})
    for key in ['energy_ha', 'density_sha256', 'iterations', 'operator_applications', 'davidson_iterations']:
        assert runs[0][key] == runs[1][key], (case, key)
    stress = runs[1]['stress']
    derivative = stress['elapsed_seconds']
    pulay = stress['pulay_assessment']['elapsed_seconds']
    denominator = runs[1]['scf_seconds'] + derivative
    records[case] = {'runs': runs, 'analytic_fraction_of_scf_plus_pass': derivative / denominator,
                     'analytic_fraction_including_pulay': derivative / (denominator + pulay)}
    print(case, 'analytic seconds', derivative, 'Pulay seconds', pulay, 'fraction', derivative / denominator, flush=True)
Path(a.output).write_text(json.dumps({'scope': 'Mac observation, T=1, one pair per case; no speedup acceptance',
                                    'binary_sha256': hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
                                    'upf_sha256': hashlib.sha256(Path(upf).read_bytes()).hexdigest(),
                                    'cases': records}, indent=2) + '\n')
