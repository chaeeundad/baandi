"""Rejection rules of scripts/parse_results.py on synthetic results documents
(no engine run; the CLI test crates/dft-cli/tests/run_config.rs runs the parser
on real dft-cli output)."""
from pathlib import Path
import copy
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from parse_results import check, main

SHA = 'ab' * 32


def bucket(label, wall, maps, wall_min=None, other=0.0):
    return {'bucket': label,
            'max': {'wall': wall, 'solve.other': other, 'maps': maps},
            'min': {'wall': wall if wall_min is None else wall_min, 'solve.other': other, 'maps': maps}}


def rank(r, wall, covered):
    return {'rank': r, 'wall_s': wall, 'covered_s': covered, 'residual_s': wall - covered}


def doc(profile=None, switches=None):
    d = {'run': {'schema': 'run-report-1', 'run_id': '0123456789abcdef', 'exit_reason': 'converged',
                 'binary': {'path': '/x/dft-cli', 'sha256': SHA}, 'completed': True},
         'scf': {'switches': {k: {'value': v} for k, v in (switches or {}).items()}},
         'convergence': {'work': {'davidson_iterations_recorded': 12,
                                  'davidson_stages': {'stage_s': {'h_apply': 0.0, 'rotate': 0.0}}}}}
    if profile is not None:
        d['phase_profile'] = profile
    return d


def profile(per_rank, buckets=None):
    covered = [r['covered_s'] for r in per_rank]
    return {'schema': 'phase-profile-2', 'stage_timers': 'on',
            'buckets': buckets or [bucket('init', 0.2, 0), bucket('fine', 1.0, 9)],
            'per_rank': per_rank,
            'closure': {'covered_max_s': max(covered), 'covered_min_s': min(covered),
                        'residual_max_s': max(r['residual_s'] for r in per_rank),
                        'residual_min_s': min(r['residual_s'] for r in per_rank)}}


class ParseResults(unittest.TestCase):
    def test_accepts_a_complete_run(self):
        self.assertEqual(check(doc(), SHA), [])
        self.assertEqual(check(doc(profile([rank(0, 1.3, 1.2)]),
                                   {'BAANDI_PHASE_PROFILE': 'on'})), [])

    def test_interrupted_run_is_rejected(self):
        d = doc()
        del d['run']['completed']
        self.assertIn('no completion marker', check(d)[0])
        del d['run']
        self.assertIn('no completion marker', check(d)[0])
        with tempfile.TemporaryDirectory() as t:
            path = Path(t) / 'results.json'
            text = json.dumps(doc(), indent=1)
            path.write_text(text[: len(text) // 2])  # killed mid-write
            out = io.StringIO()
            with redirect_stdout(out):
                code = main([str(path)])
            self.assertEqual(code, 2)
            self.assertIn('REJECT: truncated: not valid JSON', out.getvalue())
            path.write_text(text)
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main([str(path), '--expect-sha', SHA]), 0)

    def test_stale_binary_is_rejected(self):
        self.assertIn('stale binary', check(doc(), 'cd' * 32)[0])
        d = doc()
        del d['run']['binary']['sha256']
        self.assertIn('no binary SHA-256', check(d, SHA)[0])
        with tempfile.TemporaryDirectory() as t:
            m = Path(t) / 'manifest.json'
            m.write_text(json.dumps({'binary_sha256': 'cd' * 32}))
            r = Path(t) / 'results.json'
            r.write_text(json.dumps(doc()))
            with redirect_stdout(io.StringIO()):
                self.assertEqual(main([str(r), '--manifest', str(m)]), 2)

    def test_missing_or_zero_clocks_of_an_active_path(self):
        d = doc(switches={'BAANDI_PHASE_PROFILE': 'on'})
        self.assertIn('no phase_profile', check(d)[0])
        p = profile([rank(0, 1.3, 1.2)], [bucket('fine', 0.0, 9)])
        self.assertIn('zero wall', ' '.join(check(doc(p))))
        p = profile([rank(0, 1.3, 0.0), rank(1, 1.3, 0.0)], [bucket('fine', 0.0, 0)])
        self.assertIn("every rank's covered time is zero", ' '.join(check(doc(p))))
        # Count-only mode: zero clocks are expected, not an error.
        p['stage_timers'] = 'off'
        self.assertEqual(check(doc(p)), [])
        d = doc(switches={'BAANDI_DAVIDSON_STAGES': 'on'})
        self.assertIn('every stage second zero', check(d)[0])
        d['scf']['switches']['BAANDI_STAGE_TIMERS'] = {'value': 'off'}
        self.assertEqual(check(d), [])
        d = doc()
        d['parallel'] = {'stage_timers': 'on', 'collectives': [{'op': 'tree_sum', 'calls': 40, 'seconds': 0.0}]}
        self.assertIn('collective tree_sum', check(d)[0])

    def test_impossible_rows(self):
        p = profile([rank(0, 1.3, 1.2)], [bucket('fine', 1.0, 9, wall_min=1.5)])
        self.assertIn('min 1.5 > max 1.0', check(doc(p))[0])
        p = profile([rank(0, 1.3, 1.2)], [bucket('fine', 1.0, 9, other=-0.2)])
        self.assertIn('solve.other negative', check(doc(p))[0])

    def test_rank_imbalance_closes_per_rank(self):
        # Simulated imbalance: the legacy sum of per-bucket maxima (1.0 + 3.0 +
        # 2.0 ...) can exceed every wall, the per-rank closure cannot unless a
        # rank double counts.
        ok = profile([rank(0, 4.5, 4.0), rank(1, 4.0, 3.5), rank(2, 4.4, 2.0)])
        self.assertEqual(check(doc(ok)), [])
        bad = profile([rank(0, 4.5, 4.0), rank(1, 4.0, 4.6)])
        problems = check(doc(bad))
        self.assertTrue(any('rank 1 covered 4.6 s exceeds its wall 4.0 s' in p for p in problems), problems)
        # Within clock granularity: accepted.
        fine = profile([rank(0, 4.0, 4.0005)])
        self.assertEqual(check(doc(fine)), [])

    def test_json_verdict(self):
        with tempfile.TemporaryDirectory() as t:
            r = Path(t) / 'results.json'
            d = copy.deepcopy(doc())
            d['run']['completed'] = False
            r.write_text(json.dumps(d))
            out = io.StringIO()
            with redirect_stdout(out):
                self.assertEqual(main([str(r), '--json']), 2)
            v = json.loads(out.getvalue())
            self.assertFalse(v['accepted'])
            self.assertEqual(v['run_id'], '0123456789abcdef')


if __name__ == '__main__':
    unittest.main()
