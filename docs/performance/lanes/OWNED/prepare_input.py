#!/usr/bin/env python3
"""Render one existing benchmark case with the repository's canonical renderer."""
import argparse
import importlib.util
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True, choices=['si16', 'si32', 'si64', 'si128'])
    parser.add_argument('--upf', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[4]
    spec = importlib.util.spec_from_file_location('bench_suite', root / 'scripts/bench_suite.py')
    suite = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(suite)
    case = suite.load_case(root / 'benchmarks/cases' / (args.case + '-gamma.toml'))
    text = suite.case_input(case, args.upf.resolve(strict=True), 601, 1)
    # Same immutable input for both switch arms. CLI default standard path,
    # verification and cutoff ladder are preserved; no hand-written Si128 cell.
    args.output.write_text(text)


if __name__ == '__main__':
    main()
