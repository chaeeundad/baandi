#!/usr/bin/env python3
"""Build the prepared replay patch in scratch without changing the CLI bench.
Temporary ownership workaround: the ticket names a nonexistent solver example.
The production source and Cargo manifests are never modified by this script.
"""
import argparse
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true', help='also clippy and test the prepared bench')
    args = parser.parse_args()
    lane = Path(__file__).resolve().parent
    root = lane.parents[3]
    scratch = root / 'target/scratch'
    scratch.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TMPDIR=str(scratch))
    subprocess.run(['cargo', 'build', '--release', '-p', 'dft-cli', '--example', 'davidson_bench'],
                   cwd=root, env=env, check=True)
    source = scratch / 'prepared_filter_bench.rs'
    subprocess.run(['patch', '--batch', '-o', str(source),
                    str(root / 'crates/dft-cli/examples/davidson_bench.rs'), str(lane / 'bench.patch')],
                   cwd=root, env=env, check=True)
    libs = ['dft_cli', 'dft_engine', 'dft_geometry', 'dft_mpi', 'dft_physics', 'serde_json',
            'dft_numerics', 'dft_solvers', 'dft_basis', 'dft_units', 'dft_upf']
    common = ['--edition=2024', '-L', f'dependency={root}/target/release/deps', str(source)]
    for lib in libs:
        files = list((root / 'target/release/deps').glob(f'lib{lib}-*.rlib'))
        common += ['--extern', f'{lib}={max(files, key=lambda p: p.stat().st_mtime)}']
    binary = scratch / 'prepared_filter_bench'
    subprocess.run(['rustc'] + common + ['-C', 'opt-level=3', '-o', str(binary)],
                   cwd=root, env=env, check=True)
    if args.check:
        subprocess.run(['clippy-driver'] + common + ['--emit=metadata', '-o', str(scratch / 'prepared-meta'),
                                                    '-D', 'warnings'], cwd=root, env=env, check=True)
        test = scratch / 'prepared_filter_tests'
        subprocess.run(['rustc'] + common + ['--test', '-C', 'opt-level=2', '-o', str(test)],
                       cwd=root, env=env, check=True)
        subprocess.run([str(test)], cwd=root, env=env, check=True)
    print(binary)


if __name__ == '__main__':
    main()
