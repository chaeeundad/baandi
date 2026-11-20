#!/usr/bin/env python3
"""Parser transport-contract tests; no timing/performance inference."""
import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('lane', Path(__file__).resolve().parents[1] / 'lane_davidson_parse.py')
lane = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lane)


class ParseTests(unittest.TestCase):
    def setUp(self):
        self.sample = dict(run=dict(completed=True), davidson_profile=dict(
            schema='davidson-profile-1', fields=lane.FIELDS,
            operations=['solve', 'reduced_gram', 'tree_sum'], stages=['fine'],
            kinds=['ordinary'], reduction_sites=['reduced_gram'],
            per_rank=[dict(rank=0, dropped=0, records=[
                [1, 1, 1, 0, 0, 2, 0, 100, 8, 6, .1, 0, 1e-8],
                [2, 1, 1, 0, 0, 2, 1, 100, 16, 8, .02, 0, 0],
                [3, 1, 1, 0, 0, 2, 2, 128, 1, 0, .001, 0, 0]])]))

    def parse(self, sample):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'results.json'
            p.write_text(json.dumps(sample))
            return lane.parse(p)

    def test_tables_and_floor_inputs(self):
        r = self.parse(self.sample)
        self.assertEqual(len(r['calls']), 3)
        self.assertEqual(r['calls'][1]['gemm_flop_model'], 25600)
        self.assertEqual(r['calls'][2]['reduction_payload_bytes'], 1024)
        self.assertEqual(len(r['map_rank_table']), 3)
        self.assertEqual(len(r['map_table']), 3)
        self.assertEqual(r['map_table'][0]['participating_ranks'], [0])

    def test_rejects_incomplete_corrupt_and_truncated_profiles(self):
        mutations = [
            lambda x: x['run'].update(completed=False),
            lambda x: x.pop('davidson_profile'),
            lambda x: x['davidson_profile']['per_rank'][0].update(dropped=1),
            lambda x: x['davidson_profile']['per_rank'][0].update(rank=1),
            lambda x: x['davidson_profile']['per_rank'][0]['records'].append(
                x['davidson_profile']['per_rank'][0]['records'][0]),
            lambda x: x['davidson_profile']['per_rank'][0]['records'][1].__setitem__(10, float('nan')),
            lambda x: x['davidson_profile']['per_rank'][0]['records'][1].__setitem__(7, -1),
            lambda x: x['davidson_profile']['per_rank'][0]['records'][1].__setitem__(6, 99),
            lambda x: x['davidson_profile']['per_rank'][0]['records'].pop(0),
        ]
        for mutate in mutations:
            sample = copy.deepcopy(self.sample)
            mutate(sample)
            with self.assertRaises(ValueError):
                self.parse(sample)


if __name__ == '__main__':
    unittest.main()
