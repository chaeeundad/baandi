#!/usr/bin/env python3
"""Validate saved measurement evidence independently of the summary generator."""
import gzip
import json
import math
from pathlib import Path
import unittest

OUT = Path(__file__).resolve().parent / 'reference/perf/p15a'
CASES = ('si8', 'k444s', 'syn20-skew', 'syn4-elongated', 'syn4-slab', 'syn16-gamma', 'syn16-k211')


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items() if k != 'timings'}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


class Evidence(unittest.TestCase):
    def test_all_frozen_inputs_and_native_results(self):
        for name in CASES:
            with self.subTest(case=name):
                data = json.loads((OUT / (name + '.json')).read_text())
                self.assertTrue(all(data['determinism'].values()))
                self.assertEqual(data['exit_code'], data['baseline_exit_code'])
                self.assertEqual(data['report']['total_energy'], data['baseline_report']['total_energy'])
                self.assertEqual(data['density_sha256'], data['baseline_density_sha256'])
                self.assertEqual(scientific(data['report']['convergence']),
                                 scientific(data['baseline_report']['convergence']))
                self.assertEqual(data['report']['eigenvalues'], data['baseline_report']['eigenvalues'])
                traces = json.loads(gzip.decompress((OUT / data['trace_file']).read_bytes()))
                for solve in traces:
                    first = solve['0']['steps'][0]['eigenvalues_ha']
                    for variant in map(str, range(7)):
                        trace = solve[variant]
                        self.assertEqual(first, trace['steps'][0]['eigenvalues_ha'])
                        previous = None
                        for step in trace['steps']:
                            self.assertTrue(all(math.isfinite(v) and v >= 0 for v in step['residuals_ha']))
                            expected = ([abs(a-b) for a,b in zip(step['eigenvalues_ha'], previous)]
                                        if previous is not None else None)
                            self.assertEqual(step['ritz_changes_ha'], expected)
                            self.assertEqual(step['unconverged_residual'], sum(
                                r > (trace['tolerance_ha'] if j < trace['strict_bands'] else trace['buffer_tolerance_ha'])
                                for j, r in enumerate(step['residuals_ha'])))
                            previous = step['eigenvalues_ha']
                    self.assertTrue(solve['6']['success'])
                    self.assertLessEqual(max(solve['6']['steps'][-1]['residuals_ha']), 1e-10)
                offset = 0
                for h in data['report']['convergence']['history']:
                    n = len(h['operator_applications_per_k'])
                    self.assertEqual([s['0']['steps'][-1]['h_vectors'] for s in traces[offset:offset+n]],
                                     h['operator_applications_per_k'])
                    offset += n
                    if h['fixed_point_certificate']:
                        c = h['fixed_point_certificate']
                        self.assertEqual([s['0']['steps'][-1]['h_vectors'] for s in traces[offset:offset+n]],
                                         c['operator_applications_per_k'])
                        offset += n
                self.assertEqual(offset, len(traces))


if __name__ == '__main__':
    unittest.main()
