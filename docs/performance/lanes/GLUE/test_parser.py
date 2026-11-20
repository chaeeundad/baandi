"""Parser contracts, including incomplete/invalid reports and count-only data."""
import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('glue_parse', ROOT/'scripts/lane_glue_parse.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.data = dict(run=dict(completed=True, rank=1, ranks=8, counters={
            'glue_lane.serial.mix.calls': 2, 'glue_lane.serial.mix.ns': 2000000,
            'glue_lane.grid.potential.calls': 0, 'glue_lane.grid.potential.ns': 0,
            'glue_lane.dist.dot_elements': 17}))

    def parse(self, data, **kwargs):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'results.json'
            p.write_text(json.dumps(data))
            return parser.parse(p, **kwargs)

    def test_units_rank_and_work(self):
        r = self.parse(self.data)
        self.assertEqual(r['rank'], 1)
        self.assertEqual(r['spans'][0]['seconds_per_call'], 0.001)
        self.assertEqual(r['logical_work']['dist_dot_flops'], 34)
        self.assertIsNone(r['spans'][1]['seconds_per_call'])

    def test_count_only_and_mac_suppression(self):
        self.assertIsNone(self.parse(self.data, correctness_only=True)['spans'][0]['seconds'])
        self.data['run']['counters']['glue_lane.serial.mix.ns'] = 0
        r = self.parse(self.data)
        self.assertEqual(r['state'], 'counts_only')
        self.assertIsNone(r['spans'][0]['seconds'])

    def test_invalid(self):
        for kind in ('incomplete', 'missing', 'rank', 'negative', 'nan', 'float', 'orphan', 'no_calls'):
            d = copy.deepcopy(self.data)
            c = d['run']['counters']
            if kind == 'incomplete': d['run']['completed'] = False
            if kind == 'missing': c.clear()
            if kind == 'rank': d['run']['rank'] = 8
            if kind in ('negative', 'nan', 'float'):
                c['glue_lane.serial.mix.ns'] = {'negative': -1, 'nan': float('nan'), 'float': 1.5}[kind]
            if kind == 'orphan': del c['glue_lane.serial.mix.ns']
            if kind == 'no_calls': c['glue_lane.serial.mix.calls'] = 0
            with self.subTest(kind=kind), self.assertRaises(ValueError): self.parse(d)


if __name__ == '__main__': unittest.main()
