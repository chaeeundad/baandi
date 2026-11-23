#!/usr/bin/env python3
"""Run the preselected fixed-potential Mac samples, one job <=590 seconds.
Requires the replay bench built from bench.patch until its path is authorized.
No UPF or binary snapshot bytes are stored in the repository evidence.
"""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', required=True)
    parser.add_argument('--case', required=True)
    parser.add_argument('--phase', choices=['early', 'mid', 'late'])
    args = parser.parse_args()
    lane = Path(__file__).resolve().parent
    root = lane.parents[3]
    index = json.loads((lane / 'capture-index.json').read_text())
    case = next(c for c in index['cases'] if c['case'] == args.case)
    for sample in case['samples']:
        if args.phase and sample['phase'] != args.phase:
            continue
        command = [args.binary, '--input', f"target/scratch/{args.case}.toml",
                   '--snapshot', sample['manifest'], '--repeats', '1', '--filter-replay']
        result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=590)
        name = f"{args.case}-{sample['phase']}"
        (root / f'target/scratch/{name}-replay.log').write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f'{name}: exit {result.returncode}: {result.stderr[-1500:]}')
        output = next(json.loads(line[len('DAVBENCH_JSON '):])
                      for line in result.stdout.splitlines() if line.startswith('DAVBENCH_JSON '))
        replay = output['snapshot']['filter_replay']
        record = dict(command=command, exit_code=result.returncode,
                      exact_engine_replay=output['snapshot'], result=replay)
        # Avoid duplicating the large candidate journal.
        del record['exact_engine_replay']['filter_replay']
        (lane / f'{name}.json').write_text(json.dumps(record, indent=2) + '\n')
        print(name, [(c['degree'], len(c['sweeps']), c['filter_work']['h_columns'], c['fallback'])
                     for c in replay['candidates']], flush=True)


if __name__ == '__main__':
    main()
