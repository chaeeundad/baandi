#!/usr/bin/env python3
"""Independent specification tests for offline decisions, units and censoring."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('analysis', Path(__file__).with_name('p15a_analyze.py'))
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class Decisions(unittest.TestCase):
    def trace(self):
        return dict(steps=[dict(step=i+1, h_vectors=2*(i+1), h_blocks=i+1,
                               eigenvalues_ha=[v, v+1], residuals_ha=[r, 2*r])
                           for i, (v, r) in enumerate([(1., .1), (.9, .01), (.901, .02), (.9011, .001)])])

    def test_no_first_step_ritz_and_strict_comparison(self):
        trace = self.trace()
        self.assertEqual(a.stop(trace, 'ritz', 1.0)['step'], 2)
        self.assertEqual(a.stop(trace, 'ritz', .01)['step'], 3)
        self.assertIsNone(a.stop(trace, 'ritz', 1e-8))
        self.assertEqual(a.stop(trace, 'residual', .002)['step'], 4)

    def test_recrossing_is_not_locking(self):
        rates = a.rates(self.trace())['0.01']
        self.assertEqual(rates[0], dict(first=2, sustained=4, recrossed=True))
        self.assertEqual(rates[1], dict(first=4, sustained=4, recrossed=False))

    def test_allowance_and_cap_are_distinct(self):
        trace = self.trace()
        self.assertEqual(a.stop(trace, 'ritz', 1e-9, 2)['step'], 2)
        self.assertIsNone(a.stop(trace, 'ritz', 1e-9, 2, cap_only=True))
        # A native initial projection is not a QE expansion step.
        trace['steps'][1]['step'] = 21
        self.assertEqual(a.stop(trace, 'ritz', 1e-9, 2, cap_only=True)['expansion_steps'], 20)

    def test_ratchet_units_floor_and_redo(self):
        metrics = [dict(iteration=i+1, dr2_ry=d, electrons=10.)
                   for i, d in enumerate([.01, .1, 0., 1.])]
        s = a.schedule(metrics)
        self.assertTrue(s[0]['first_redo_trigger'])
        self.assertEqual(s[0]['ethr_ha'], .005)
        self.assertEqual(s[0]['redo_ethr_ha'], .00005)
        self.assertEqual(s[1]['ethr_ha'], .00005)
        self.assertEqual(s[2]['ethr_ha'], .00005)
        self.assertEqual(s[3]['ethr_ha'], 5e-14)


if __name__ == '__main__':
    unittest.main()
