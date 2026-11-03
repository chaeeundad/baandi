"""Synthetic k-point interchange checks; no QE execution or real UPF."""
import sys
from pathlib import Path
import unittest
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from make_qe_input import kpoint_input, compare_kpoints, generate
from qe_common import Error, read_json, write_json


class KPoints(unittest.TestCase):
    def block(self, grid=None, shift=None, tr=False):
        return {'kpoints': {'grid': grid, 'shift': shift, 'time_reversal': tr,
            'count': 2, 'points': [dict(fractional=[0.25, 0, 0], weight=0.5),
                                 dict(fractional=[-0.25, 0, 0], weight=0.5)]}}

    def test_legacy(self):
        card, noinv, points = kpoint_input({})
        self.assertEqual(card, ['K_POINTS automatic', '1 1 1 0 0 0'])
        self.assertTrue(noinv)
        self.assertEqual(points[0]['weight'], 1)

    def test_grid(self):
        for tr in (False, True):
            card, noinv, _ = kpoint_input(self.block([2, 3, 4], [0.5, 0, 0.5], tr))
            self.assertEqual(card, ['K_POINTS automatic', '2 3 4 1 0 1'])
            self.assertEqual(noinv, not tr)

    def test_explicit(self):
        for tr in (False, True):
            card, noinv, _ = kpoint_input(self.block(tr=tr))
            self.assertEqual(card, ['K_POINTS crystal', '2', '0.25 0 0 0.5', '-0.25 0 0 0.5'])
            self.assertTrue(noinv)

    def test_errors(self):
        for grid, shift in [([0,1,1], [0,0,0]), ([2,1,1], [0.1,0,0])]:
            with self.assertRaises(Error):
                kpoint_input(self.block(grid, shift))
        r = self.block()
        r['kpoints']['points'][0]['weight'] = 0.4
        with self.assertRaises(Error):
            kpoint_input(r)

    def test_comparison(self):
        a = [dict(fractional=[0.25,0,0], weight=1)]
        b = [dict(fractional=[0.75,0,0], weight=1)]
        self.assertFalse(compare_kpoints(a,b))
        self.assertTrue(compare_kpoints(a,b,True))
        b[0]['weight'] = 0.5
        self.assertFalse(compare_kpoints(a,b,True))

    def test_generated_flags_and_metadata(self):
        fixture = Path(__file__).resolve().parents[2]/'tests/fixtures_format/results.json'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for grid, tr in (([2,1,1], True), ([2,1,1], False), (None, True), (None, False)):
                r = read_json(fixture)
                r.update(self.block(grid, [0.5,0,0] if grid else None, tr))
                if grid and tr:
                    r['kpoints']['points'] = [dict(fractional=[0.25,0,0], weight=1)]
                    r['kpoints']['count'] = 1
                write_json(root/'results.json', r)
                meta = generate(root/'results.json', '/synthetic', 'synthetic.UPF', root/'qe.in')
                text = (root/'qe.in').read_text()
                self.assertIn('nosym = .true.', text)
                self.assertEqual('noinv = .true.' in text, not (grid and tr))
                self.assertIn('K_POINTS automatic' if grid else 'K_POINTS crystal', text)
                self.assertEqual(meta['expected_kpoints'], r['kpoints']['points'])
