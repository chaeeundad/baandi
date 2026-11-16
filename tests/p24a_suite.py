#!/usr/bin/env python3
"""Run the unchanged robustness harness with an explicit locking policy."""
import importlib.util
from pathlib import Path
import os
import sys
# ADR 0043 ONEPATH-4: locking="auto" exists only on the strict reference path.
os.environ.setdefault('BAANDI_REFERENCE_STRICT','on')
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('suite',root/'scripts/bench_suite.py')
suite=importlib.util.module_from_spec(spec);spec.loader.exec_module(suite)
original=suite.case_input
suite.case_input=lambda *args,**kwargs: original(*args,**kwargs).replace('[scf]','[scf]\nlocking="auto"')
if __name__=='__main__':sys.exit(suite.main())
