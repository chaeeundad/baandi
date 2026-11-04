#!/usr/bin/env python3
"""Paired uninstrumented P10M timing/equivalence experiment (run at repo root).

Each budget/case gets one discarded warmup per binary and >=30 alternating
pairs by default. No compilation or other benchmarks should run concurrently.
JSON retains every elapsed sample, hashes, energies and Davidson counts.
Example:
 python3 -I crates/dft-numerics/src/pool/orthonormalize_bench.py \
   /tmp/before /tmp/after /absolute/Si.pz-vbc.UPF /tmp/p10m-times.json
"""
import argparse
import importlib.util
import json
from pathlib import Path
import statistics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('before')
    ap.add_argument('after')
    ap.add_argument('upf')
    ap.add_argument('output', type=Path)
    ap.add_argument('--threads', default=','.join(map(str, range(1, 13))))
    ap.add_argument('--cases', default='si8,k444s')
    ap.add_argument('--repeats', type=int, default=30)
    args = ap.parse_args()
    spec = importlib.util.spec_from_file_location('bench', 'scripts/bench_scf.py')
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    paths = {k: str(Path(getattr(args, k)).resolve()) for k in ['before', 'after', 'upf']}
    rows = []
    for threads in map(int, args.threads.split(',')):
        for case in args.cases.split(','):
            samples = {'before': [], 'after': []}
            for name in samples:
                bench.run_case(paths[name], paths['upf'], case, 1, threads)
            for repeat in range(args.repeats):
                for name in (['before', 'after'] if repeat % 2 == 0 else ['after', 'before']):
                    samples[name].append(bench.run_case(paths[name], paths['upf'], case, 1, threads))
            reference = samples['before'][0]
            for records in samples.values():
                for r in records:
                    for key in ['total_energy_ha', 'density_sha256', 'davidson_iterations', 'certificate_iterations']:
                        assert r[key] == reference[key], (case, threads, key)
            row = dict(case=case, threads=threads, samples=samples,
                       median_s={k: statistics.median(r['wall_median_s'] for r in v)
                                 for k, v in samples.items()})
            rows.append(row)
            args.output.write_text(json.dumps(dict(repeats=args.repeats, binaries=paths, rows=rows), indent=2) + '\n')
            print(case, threads, row['median_s'], flush=True)


if __name__ == '__main__':
    main()
