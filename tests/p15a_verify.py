#!/usr/bin/env python3
"""Run unchanged bench_scf identity gates at the performance-brief thread budgets."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'tests/reference/perf/p15a'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--before', default='/tmp/p15a-before-cli')
    ap.add_argument('--after', default=str(ROOT / 'target/release/dft-cli'))
    ap.add_argument('--upf', default=str(ROOT.parent / 'pseudo/Si.pz-vbc.UPF'))
    args = ap.parse_args()
    commands = []
    for threads in (1, 2, 4, 8, 12):
        before = OUT / f'default-before-t{threads}.json'
        after = OUT / f'default-after-t{threads}.json'
        common = [sys.executable, '-I', str(ROOT / 'scripts/bench_scf.py'), '--upf', args.upf,
                  '--cases', 'gamma,k222s,k222,k444s,si8', '--threads', str(threads), '--repeats', '1']
        for binary, output, extra in [(args.before, before, []), (args.after, after,
                                     ['--baseline', str(before), '--require-identical'])]:
            command = common + ['--binary', binary, '--record', str(output)] + extra
            proc = subprocess.run(command, capture_output=True, text=True)
            commands.append(dict(command=command, exit_code=proc.returncode, stdout=proc.stdout, stderr=proc.stderr))
            (OUT / 'identity-commands.json').write_text(json.dumps(commands, indent=2)+'\n')
            print('threads', threads, output.name, 'exit', proc.returncode, flush=True)
            if proc.returncode:
                raise SystemExit(proc.returncode)
    reference = json.loads((OUT / 'default-after-t1.json').read_text())['cases']
    for threads in (2, 4, 8, 12):
        current = json.loads((OUT / f'default-after-t{threads}.json').read_text())['cases']
        for key, before in reference.items():
            for field in ('total_energy_ha', 'density_sha256', 'iterations', 'davidson_iterations',
                          'operator_applications', 'certificate_applications', 'certificate_iterations'):
                assert before[field] == current[key][field], (threads, key, field)
    print('Cross-thread energy/density/work/certificate identity: PASS', flush=True)


if __name__ == '__main__':
    main()
