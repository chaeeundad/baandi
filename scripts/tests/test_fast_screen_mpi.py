"""Decision logic of scripts/fast_screen_mpi.py on synthetic records (no engine run)."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from fast_screen_mpi import verdict


def record(wall, tree_calls, energy=-126.7, bits='0x1', fnv='0xa', panel_norm=100, barrier=0):
    return {
        'wall_s': wall, 'energy_ha': energy, 'energy_bits': bits, 'density_fnv': fnv,
        'scf': 13, 'davidson': 74, 'h_columns': 2884,
        'reductions': {'panel_norm': [panel_norm, panel_norm], 'setup': [0, 0]},
        'collectives': {'tree_sum': [tree_calls, 8 * tree_calls, 0.1], 'barrier': [barrier, 0, 0.0]},
        'stages': {'orthonormalize': 0.5},
    }


class Verdicts(unittest.TestCase):
    def test_counter_drop_with_wall_in_spread_is_candidate(self):
        arms = {'A': [record(4.2, 1000), record(4.4, 1000), record(4.3, 1000)],
                'B': [record(4.25, 600, bits='0x2', energy=-126.7 + 1e-12), record(4.3, 600, bits='0x2', energy=-126.7 + 1e-12),
                      record(4.35, 600, bits='0x2', energy=-126.7 + 1e-12)]}
        res = verdict(arms, [('collectives.tree_sum.calls', 'down')], 5e-9, 0.10)
        self.assertEqual(res['verdict'], 'CANDIDATE')
        self.assertEqual(res['correctness']['identity'], 'rounding class')
        self.assertIn('collectives.barrier.calls', res['unusable'])
        self.assertIn('reductions.setup.calls', res['unusable'])

    def test_counter_not_moving_is_reject(self):
        arms = {'A': [record(4.2, 1000)] * 3, 'B': [record(4.0, 1000)] * 3}
        res = verdict(arms, [('collectives.tree_sum.calls', 'down')], 5e-9, 0.10)
        self.assertEqual(res['verdict'], 'REJECT')

    def test_energy_outside_class_is_reject(self):
        arms = {'A': [record(4.2, 1000)] * 3, 'B': [record(3.0, 500, energy=-126.7 + 1e-8)] * 3}
        res = verdict(arms, [], 5e-9, 0.10)
        self.assertEqual(res['verdict'], 'REJECT')

    def test_all_zero_expected_counter_is_unusable(self):
        arms = {'A': [record(4.2, 1000)] * 3, 'B': [record(4.2, 1000)] * 3}
        res = verdict(arms, [('collectives.barrier.calls', 'down')], 5e-9, 0.10)
        self.assertEqual(res['verdict'], 'UNUSABLE')
        self.assertEqual(res['expect_unusable'], ['collectives.barrier.calls'])

    def test_identical_arms_are_neutral_and_bit_identical(self):
        arms = {'A': [record(4.2, 1000), record(4.3, 1000)], 'B': [record(4.25, 1000), record(4.22, 1000)]}
        res = verdict(arms, [], 5e-9, 0.10)
        self.assertEqual(res['verdict'], 'NEUTRAL')
        self.assertEqual(res['correctness']['identity'], 'bit-identical')

    def test_wall_beyond_spread(self):
        arms = {'A': [record(4.2, 1000), record(4.25, 1000)], 'B': [record(3.5, 1000), record(3.52, 1000)]}
        self.assertEqual(verdict(arms, [], 5e-9, 0.10)['verdict'], 'CANDIDATE')
        arms = {'A': [record(3.5, 1000), record(3.52, 1000)], 'B': [record(4.2, 1000), record(4.25, 1000)]}
        self.assertEqual(verdict(arms, [], 5e-9, 0.10)['verdict'], 'REJECT')

    def test_nonreproducible_density_is_reject(self):
        arms = {'A': [record(4.2, 1000)] * 2, 'B': [record(4.2, 900, fnv='0xb'), record(4.2, 900, fnv='0xc')]}
        self.assertEqual(verdict(arms, [], 5e-9, 0.10)['verdict'], 'REJECT')


if __name__ == '__main__':
    unittest.main()
