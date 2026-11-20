#!/usr/bin/env python3
"""Synthetic schema/statistics tests; fabricated values, never performance data."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


def load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Tools(unittest.TestCase):
    def test_border_prediction_and_overflow(self):
        counter = load('border-counter')
        fields = ['operation', 'n', 'm', 'k', 'detail']
        rows = [[0, 3, 132, 132, 0], [1, 3, 64, 64, 0],
                [1, 3, 128, 64, 0], [1, 3, 132, 4, 0]]
        trace = dict(fields=fields, operations=['upper_border_product', 'upper_border_tile'],
                     per_rank=[dict(rank=0, dropped=0, records=rows)])
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'results.json'
            path.write_text(json.dumps(dict(davidson_profile=trace)))
            c = counter.counters(path)['ranks'][0]
            self.assertEqual(c['full_rectangle_flop'], 104544)
            self.assertEqual(c['executed_flop'], 76896)
            self.assertEqual(c['executed_flop'], c['predicted_tiled_flop'])
            trace['per_rank'][0]['dropped'] = 1
            path.write_text(json.dumps(dict(davidson_profile=trace)))
            with self.assertRaises(ValueError):
                counter.counters(path)

    def test_crossed_pair_statistics(self):
        summary = load('paired-summary')
        with tempfile.TemporaryDirectory() as temp:
            job = Path(temp)
            for phase in ('aa', 'ab'):
                for pair in range(1, 4):
                    for arm in ('a', 'b'):
                        seconds = 9 if (phase, arm) == ('ab', 'b') else 10
                        stem = job / f'{phase}-{pair}-{arm}'
                        stem.mkdir()
                        (stem / 'results.json').write_text(json.dumps(dict(run=dict(completed=True))))
                        Path(str(stem) + '.external-s').write_text(str(seconds))
            result = summary.summarize(job, 'screen')
            self.assertTrue(result['statistical_gain'])
            self.assertTrue(result['exceeds_three_times_aa_floor'])
            self.assertAlmostEqual(result['phases']['ab']['mean'], -.1)
            self.assertEqual(result['phases']['ab']['se'], 0)


if __name__ == '__main__':
    unittest.main()
