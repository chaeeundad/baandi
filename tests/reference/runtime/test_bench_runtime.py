"""Acceptance-statistics regressions, independent of expensive SCF runs."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('bench_runtime', Path(__file__).resolve().parents[3] / 'scripts/bench_runtime.py')
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)


class StatisticsTests(unittest.TestCase):
    def test_paired_not_ratio_of_unpaired_medians(self):
        result = bench.compare([1, 10, 100], [2, 10, 50], 1, 7, True)
        self.assertEqual(result['paired_overhead_percent'], [100, 0, -50])
        self.assertEqual(result['status'], 'inconclusive')

    def test_strict_threshold_and_insufficient_samples(self):
        self.assertEqual(bench.compare([100]*20, [101]*20, 1, 7, True)['status'], 'not_met')
        self.assertEqual(bench.compare([100]*20, [100]*20, 1, 7, True)['status'], 'met')
        self.assertEqual(bench.compare([100], [1], 1, 7, False)['status'], 'insufficient_samples')

    def test_repeatable_bootstrap(self):
        values = [98, 101, 99, 100, 102] * 4
        self.assertEqual(bench.compare([100]*20, values, 1, 2026, True),
                         bench.compare([100]*20, values, 1, 2026, True))


if __name__ == '__main__':
    unittest.main()
