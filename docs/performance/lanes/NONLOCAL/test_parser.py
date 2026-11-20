#!/usr/bin/env python3
"""Contract tests for the NONLOCAL result parser; no timing measurement."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('lane_parser', ROOT / 'scripts/lane_nonlocal_parse.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def fixture(self):
        low, high = {'maps': 2}, {'maps': 2}
        for stage in parser.STAGES:
            key = 'nonlocal.' + stage
            low.update({key: 1., key+'.calls': 2, key+'.flops': 100, key+'.bytes': 80})
            high.update({key: 3., key+'.calls': 4, key+'.flops': 200, key+'.bytes': 160})
        for shape in parser.SHAPES:
            low['nonlocal.' + shape] = 0
            high['nonlocal.' + shape] = 0
        return {'phase_profile': {'ranks': 8, 'stage_timers': 'on', 'buckets': [
            {'bucket': 'fine', 'min': low, 'max': high}]}}

    def parse(self, value, machine='mac'):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'results.json'
            path.write_text(json.dumps(value))
            return parser.parse_result(path, machine)

    def test_mac_omits_times(self):
        rows = self.parse(self.fixture())
        self.assertEqual(len(rows), 6)
        self.assertFalse(any('second' in key for row in rows for key in row))

    def test_rank_extrema_form_bounds(self):
        row = self.parse(self.fixture(), 'dev08')[0]
        self.assertEqual(row['seconds_per_call_lower'], .25)
        self.assertEqual(row['seconds_per_call_upper'], 1.5)
        self.assertEqual(row['seconds_per_map_lower'], .5)
        self.assertEqual(row['seconds_per_map_upper'], 1.5)

    def test_missing_profile_and_nan_reject(self):
        with self.assertRaises(ValueError):
            self.parse({})
        value = self.fixture()
        value['phase_profile']['buckets'][0]['min']['nonlocal.back.calls'] = float('nan')
        with self.assertRaises(ValueError):
            self.parse(value)

    def test_actual_cli_fixture(self):
        path = Path(__file__).with_name('mac-counts-results.json')
        rows = parser.parse_result(path, 'mac')
        self.assertEqual(len(rows), 48)
        fine = next(r for r in rows if r['bucket']=='fine' and r['stage']=='projection')
        self.assertGreater(fine['calls_min'], 0)
        self.assertGreater(fine['flops_min'], 0)
        self.assertEqual(fine['calls_min'], sum(fine[k + '_min'] for k in parser.SHAPES[3:]))
        self.assertTrue(all(r['status']=='unusable_zero_calls' for r in rows if r['stage']=='reduce'))

    def test_zero_calls_are_unusable(self):
        value = self.fixture()
        for bound in ['min', 'max']:
            value['phase_profile']['buckets'][0][bound]['nonlocal.back.calls'] = 0
        row = self.parse(value, 'dev08')[3]
        self.assertEqual(row['status'], 'unusable_zero_calls')
        self.assertIsNone(row['seconds_per_call_lower'])
        self.assertIsNone(row['seconds_per_call_upper'])


if __name__ == '__main__':
    unittest.main()
