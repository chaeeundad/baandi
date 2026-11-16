"""Decision logic of scripts/fast_screen_dav.py on synthetic records (no engine run)."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from fast_screen_dav import verdict


def kind(name, ms, h=190, dav=3, gemm=4.7e9, fnv='0x1', eig=(-0.5, -0.4, 0.1), tol=5e-3,
         converged=3, extra=None, identical=True):
    counters = {'gemm_flop': gemm, 'gemm_calls': 110, 'coll.tree_sum.calls': 0, 'dist.rotate.calls': 0,
                'lean.proj_cache_columns': 0}
    counters.update(extra or {})
    return {'kind': name, 'tolerance_ha': tol, 'solve_ms': [ms, ms * 1.01, ms * 0.99],
            'solve_ms_median': ms, 'davidson_iterations': dav, 'h_columns': h, 'converged': converged,
            'eigen_fnv': fnv, 'repeats_identical': identical, 'counters_deterministic': True,
            'eigenvalues': list(eig), 'band_energy_occupied': sum(eig[:2]), 'residual_max': 1e-3,
            'counters': counters}


def rec(*kinds):
    return {'occupied': 2, 'kinds': list(kinds)}


class Verdicts(unittest.TestCase):
    def test_identical_arms_are_neutral(self):
        arms = {'A': [rec(kind('cold', 300), kind('warm', 450, tol=1e-10))] * 3,
                'B': [rec(kind('cold', 301), kind('warm', 449, tol=1e-10))] * 3}
        res = verdict(arms, [], 0.05, 0.01)
        self.assertEqual(res['verdict'], 'NEUTRAL')
        self.assertTrue(res['kinds']['cold']['identity']['bit_identical'])
        self.assertIn('coll.tree_sum.calls', res['kinds']['cold']['unusable'])

    def test_counter_drop_inside_spread_is_candidate(self):
        arms = {'A': [rec(kind('cold', 300))] * 3,
                'B': [rec(kind('cold', 299, h=150, fnv='0x2', eig=(-0.5, -0.4 + 1e-6, 0.1)))] * 3}
        res = verdict(arms, [], 0.05, 0.01)
        self.assertEqual(res['verdict'], 'CANDIDATE')
        self.assertFalse(res['kinds']['cold']['identity']['bit_identical'])

    def test_slower_is_reject(self):
        arms = {'A': [rec(kind('cold', 300))] * 3, 'B': [rec(kind('cold', 330))] * 3}
        self.assertEqual(verdict(arms, [], 0.05, 0.01)['verdict'], 'REJECT')

    def test_eigenvalue_beyond_gate_is_reject(self):
        arms = {'A': [rec(kind('warm', 300, tol=1e-10))] * 3,
                'B': [rec(kind('warm', 200, tol=1e-10, fnv='0x3', eig=(-0.5, -0.4 + 1e-8, 0.1)))] * 3}
        self.assertEqual(verdict(arms, [], 0.05, 0.01)['verdict'], 'REJECT')

    def test_not_reproducible_is_reject(self):
        arms = {'A': [rec(kind('cold', 300))] * 3, 'B': [rec(kind('cold', 300, identical=False))] * 3}
        self.assertEqual(verdict(arms, [], 0.05, 0.01)['verdict'], 'REJECT')

    def test_expected_counter_zero_on_both_arms_is_unusable(self):
        arms = {'A': [rec(kind('cold', 300))] * 3, 'B': [rec(kind('cold', 300))] * 3}
        res = verdict(arms, [('lean.proj_cache_columns', 'up')], 0.05, 0.01)
        self.assertEqual(res['verdict'], 'UNUSABLE')

    def test_expected_counter_wrong_direction_is_reject(self):
        arms = {'A': [rec(kind('cold', 300, extra={'coll.tree_sum.calls': 100}))] * 3,
                'B': [rec(kind('cold', 300, extra={'coll.tree_sum.calls': 120}))] * 3}
        res = verdict(arms, [('coll.tree_sum.calls', 'down')], 0.05, 0.01)
        self.assertEqual(res['verdict'], 'REJECT')

    def test_fewer_than_three_pairs_does_not_judge_time(self):
        arms = {'A': [rec(kind('cold', 300))] * 2, 'B': [rec(kind('cold', 200))] * 2}
        self.assertEqual(verdict(arms, [], 0.05, 0.01)['verdict'], 'NEUTRAL')


if __name__ == '__main__':
    unittest.main()
