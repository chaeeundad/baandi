"""Synthetic contract tests only, not dev08 or production instrumentation evidence."""
import copy
import runpy
import unittest
from pathlib import Path
parse = runpy.run_path(str(Path(__file__).with_name('wait_parse.py')))['parse']


def fixture():
    ranks = []
    for rank, entry in enumerate([100, 160, 120]):
        ranks.append({'rank': rank, 'maps': [{'map_id': 1, 'bucket': 'fine',
            'start_ns': 0, 'layout': {'local_sticks': rank+2}, 'events': [
                {'seq': 0, 'op': 'alltoall', 'site': 'fft.inverse', 'payload_bytes': 512,
                 'enter_ns': entry, 'exit_ns': 200,
                 'preceding_stages': [{'stage': 'fft0', 'start_ns': 0, 'end_ns': entry,
                     'points': 20, 'logical_bytes': 640, 'flops': 400}]}]}]})
    return {'schema': 'arrival-skew-1', 'mode': 'timed', 'complete': True,
            'rank_count': 3, 'clock_valid': True, 'clock_error_ns': 1, 'ranks': ranks}


class WaitParse(unittest.TestCase):
    def test_known_arrival_and_exclusive_stage(self):
        r = parse(fixture())[0]
        self.assertEqual(r['last_rank'], 1)
        self.assertEqual(r['ranks'][0]['arrival_gap_bounds_ns'], [58, 62])
        self.assertEqual(r['last_rank_stage_excess_vs_other_median_ns']['fft0'], 50)
        self.assertEqual(r['last_vs_earliest_decomposition_ns']['local_interval'], 60)

    def test_clock_tie_and_early_completion_are_not_wait(self):
        d = fixture()
        d['clock_error_ns'] = 25
        d['ranks'][0]['maps'][0]['events'][0]['exit_ns'] = 110
        r = parse(d)[0]
        self.assertIsNone(r['last_rank'])
        self.assertEqual(r['ranks'][0]['pre_last_entry_residence_ns'], 10)
        self.assertNotIn('wait_ns', r['ranks'][0])

    def test_counts_and_rejections(self):
        d = fixture()
        d['mode'] = 'counts'
        self.assertIsNone(parse(d)[0]['last_rank'])
        for mutate in [lambda x: x.update(complete=False),
                       lambda x: x.update(clock_valid=False),
                       lambda x: x['ranks'].pop(),
                       lambda x: x['ranks'][0]['maps'][0]['events'][0].update(site='mismatch'),
                       lambda x: x['ranks'][0]['maps'][0]['events'][0].update(enter_ns=-1),
                       lambda x: x['ranks'][0]['maps'][0]['events'][0]['preceding_stages'].append(
                           copy.deepcopy(x['ranks'][0]['maps'][0]['events'][0]['preceding_stages'][0]))]:
            d = fixture()
            mutate(d)
            with self.assertRaises(ValueError): parse(d)


if __name__ == '__main__':
    unittest.main()
