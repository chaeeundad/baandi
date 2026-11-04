#!/usr/bin/env python3
"""Evidence integrity and standard CSV checks; no SCF runs in this test."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

BASE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('validation', BASE / 'validate.py')
validation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validation)


class EvidenceTests(unittest.TestCase):
    def test_predicates_audit_reference_and_csv(self):
        for name in ('measured', 'extended-skewed', 'extended-slab-offset', 'extended-slab'):
            data, _ = validation.validate(BASE / name)
            with (BASE / name / 'traces.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            expected = [(r, h) for r in data['runs'] for h in r['history']]
            self.assertEqual(len(rows), len(expected))
            for row, (r, h) in zip(rows, expected):
                self.assertEqual(len(row), 7)
                self.assertEqual(row['model'], r['model'])
                self.assertEqual(row['policy'], r['policy'])
                self.assertEqual(row['arm'], r['arm'])
                self.assertEqual(int(row['iteration']), h['iteration'])
                self.assertEqual(int(row['policy_revision']), h['revision'])
                self.assertEqual(float(row['energy_hartree']), h['energy_hartree'])
                self.assertEqual(float(row['density_residual_bohr_neg_3_over_2']), h['residual'])

    def test_tampered_recovery_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp) / 'copy'
            shutil.copytree(BASE / 'extended-skewed', directory)
            data = json.loads((directory / 'results.json').read_text())
            r = next(r for r in data['runs'] if r['arm'] == 'intervention')
            # Change both copies to bypass duplicate-file checking; the independent
            # reference comparison must still detect the fabricated agreement.
            r['energy_error_hartree'] = 0.1
            (directory / ('--'.join(validation.key(r)) + '.json')).write_text(json.dumps(r))
            (directory / 'results.json').write_text(json.dumps(data))
            with self.assertRaises(AssertionError):
                validation.validate(directory)


if __name__ == '__main__':
    unittest.main()
