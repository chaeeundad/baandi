#!/usr/bin/env python3
"""Run the unchanged robustness harness with an explicit locking policy."""
import importlib.util
from pathlib import Path
import sys
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('suite',root/'scripts/bench_suite.py')
suite=importlib.util.module_from_spec(spec);spec.loader.exec_module(suite)
original=suite.case_input
suite.case_input=lambda *args,**kwargs: original(*args,**kwargs).replace('[scf]','[scf]\nlocking="auto"')
if __name__=='__main__':sys.exit(suite.main())
