#!/usr/bin/env python3
"""Synthetic paired-statistics and analytic counter model tests, no DFT timing."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('fusion_stats', Path(__file__).with_name('fusion_stats.py'))
stats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stats)


class FusionTools(unittest.TestCase):
    def data(self, savings):
        def record(wall):
            return dict(exit=0, converged=True, verified=True, density_sha256='same', energy_ha=-1., wall_s=wall)
        return dict(records=dict(A=[record(10.) for _ in savings], B=[record(10. - x) for x in savings]))

    def test_resolved_gain_and_noise(self):
        self.assertTrue(stats.summarize(self.data([.8, .9, .85, .88]))['resolved_positive'])
        self.assertFalse(stats.summarize(self.data([-.1, .1, -.1, .1]))['resolved_positive'])
        self.assertFalse(stats.summarize(self.data([.8, .8, .8]))['resolved_positive'])
        self.assertFalse(stats.summarize(self.data([-.8, -.9, -.85]))['resolved_positive'])

    def test_scientific_gate(self):
        data = self.data([.8, .9, .85])
        data['records']['B'][0]['energy_ha'] += 1e-12
        with self.assertRaises(AssertionError):
            stats.summarize(data)

    def test_model_bytes(self):
        points = 72**3 * 128 * 15
        self.assertEqual(points, 716636160)
        self.assertEqual(32 * points - 16 * points, 11466178560)


if __name__ == '__main__':
    unittest.main()
