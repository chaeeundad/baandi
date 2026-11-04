#!/usr/bin/env python3
"""Export all scalar SCF histories as standard CSV, with units in column names."""
import csv
import json
from pathlib import Path
import sys

for arg in sys.argv[1:]:
    directory = Path(arg)
    data = json.loads((directory / 'results.json').read_text())
    with (directory / 'traces.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['model', 'policy', 'arm', 'iteration', 'policy_revision',
                         'energy_hartree', 'density_residual_bohr_neg_3_over_2'])
        for r in data['runs']:
            for h in r.get('history', []):
                writer.writerow([r['model'], r['policy'], r['arm'], h['iteration'], h['revision'],
                                 h['energy_hartree'], h['residual']])
