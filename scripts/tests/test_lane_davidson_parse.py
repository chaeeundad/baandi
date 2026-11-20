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

    def test_eigen_counts_use_map_and_solve_ids(self):
        p = self.sample['davidson_profile']
        p['operations'] += ['reduced_eigen_step', 'reduced_eigen', 'bordered_update']
        rows = p['per_rank'][0]['records']
        rows += [
            [4, 1, 1, 0, 0, 2, 3, 8, 0, 8, .01, 0, 0],
            [5, 1, 1, 0, 0, 2, 4, 8, 8, 8, .009, 0, 0],
            [6, 1, 1, 0, 0, 2, 3, 9, 8, 9, .01, 0, 0],
            [7, 1, 1, 0, 0, 2, 5, 9, 8, 9, .009, 0, 0],
            [8, 2, 2, 0, 0, 3, 0, 100, 8, 6, .1, 0, 1e-8],
            [9, 2, 2, 0, 0, 3, 3, 12, 0, 12, .01, 0, 0],
            [10, 2, 2, 0, 0, 3, 4, 12, 12, 12, .009, 0, 0],
        ]
        a, b = self.parse(self.sample)['eigen_map_work']
        self.assertEqual(a['solve_ids'], [1])
        self.assertEqual(a['reduced_steps'], 2)
        self.assertEqual(a['full_evd_calls'], 1)
        self.assertEqual(a['bordered_attempts'], 1)
        self.assertEqual(a['full_evd_sum_m3'], 512)
        self.assertEqual(b['map'], 2)
        self.assertEqual(b['solve_ids'], [2])
        self.assertEqual(b['full_evd_order_counts'], {12: 1})
        self.assertEqual(b['full_evd_per_solve'], 1)

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
