"""Independent analytic fixtures and repository-only QE density regression."""
import cmath
import hashlib
import itertools
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from qe_rhog import FormatError, dft, read_qe, validate_grid, write_qe
from compare_rhog import DEFAULTS, compare


class ReaderTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'charge-density.dat'
        self.values = {(0, 0, 0): 1+0j, (1, -1, 0): .2-.3j}
        write_qe(self.path, self.values)
        self.raw = self.path.read_bytes()

    def test_round_trip(self):
        values, meta = read_qe(self.path)
        self.assertEqual(values, self.values)
        self.assertEqual(meta['ngm_g'], 2)
        self.assertEqual(meta['nspin'], 1)

    def reject(self, raw, message):
        self.path.write_bytes(raw)
        with self.assertRaisesRegex(FormatError, message):
            read_qe(self.path)

    def test_markers_and_sizes_each_record(self):
        for offset in [0, 20, 100, 132]:
            size, = struct.unpack_from('<i', self.raw, offset)
            with self.subTest(record_offset=offset, problem='leading'):
                raw = bytearray(self.raw)
                struct.pack_into('<i', raw, offset, size+1)
                self.reject(raw, 'wrong size')
            with self.subTest(record_offset=offset, problem='trailing'):
                raw = bytearray(self.raw)
                struct.pack_into('<i', raw, offset+4+size, size+1)
                self.reject(raw, 'marker mismatch')

    def test_truncation(self):
        for end in range(len(self.raw)):
            with self.subTest(end=end):
                self.reject(self.raw[:end], 'missing leading marker|truncated payload or marker')

    def test_trailing_bytes(self):
        self.reject(self.raw+b'x', 'trailing bytes')

    def test_gamma(self):
        raw = bytearray(self.raw)
        struct.pack_into('<i', raw, 4, 1)
        self.reject(raw, 'gamma_only')

    def test_spin(self):
        raw = bytearray(self.raw)
        struct.pack_into('<i', raw, 12, 2)
        self.reject(raw, 'nspin')

    def test_count(self):
        raw = bytearray(self.raw)
        struct.pack_into('<i', raw, 8, -1)
        self.reject(raw, 'ngm_g must be positive')

    def test_duplicate(self):
        raw = bytearray(self.raw)
        raw[116:128] = raw[104:116]
        self.reject(raw, 'duplicate Miller')

    def test_nonfinite(self):
        raw = bytearray(self.raw)
        struct.pack_into('<d', raw, 136, float('nan'))
        self.reject(raw, 'finite number')


class DftTests(unittest.TestCase):
    def test_analytic_modes_and_parseval(self):
        grid = (4, 5, 6)
        points = list(itertools.product(*(range(n) for n in grid)))
        requested = list(itertools.product(*(range(-(n//2), (n-1)//2+1) for n in grid)))
        for modes in [{(1, -1, 2): .3+.2j},
                      {(1, -1, 2): .3+.2j, (-2, 1, -3): -.7+.1j}]:
            with self.subTest(modes=modes):
                values = [sum(c*cmath.exp(2j*math.pi*sum(m*x/n for m, x, n in zip(mode, p, grid)))
                              for mode, c in modes.items()) for p in points]
                actual = dft(values, grid, requested)
                for m in requested:
                    self.assertLess(abs(actual[m]-modes.get(m, 0)), 1e-12)
                self.assertAlmostEqual(sum(abs(z)**2 for z in actual.values()),
                                       sum(abs(z)**2 for z in values)/len(values), delta=1e-12)
                selected = dft(values, grid, modes)
                for m in modes:
                    self.assertLess(abs(selected[m]-modes[m]), 1e-12)

    def test_alias_and_positive_nyquist_rejected(self):
        for indices in [[(2, 0, 0)], [(0, 0, 0), (0, 0, 0)], [(-3, 0, 0)]]:
            with self.subTest(indices=indices), self.assertRaises(FormatError):
                validate_grid((4, 5, 6), indices)


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.dir = Path(self.temp.name)
        self.ours = self.dir / 'ours'
        self.ours.mkdir()
        self.grid = [5, 5, 5]
        self.values = {(0, 0, 0): 1+0j, (1, 0, 0): .05+0j, (-1, 0, 0): .05+0j,
                       (0, 1, 0): -.03j, (0, -1, 0): .03j, (0, 0, 1): 0j, (0, 0, -1): 0j}
        density = [1+.1*math.cos(2*math.pi*i/5)+.06*math.sin(2*math.pi*j/5)
                   for i, j, k in itertools.product(range(5), repeat=3)]
        raw = struct.pack('<125d', *density)
        (self.ours / 'density.bin').write_bytes(raw)
        volume = (2*math.pi)**3
        lattice = [[2*math.pi if i == j else 0 for j in range(3)] for i in range(3)]
        self.meta = dict(schema_version=1, dtype='f64', endianness='little', fft_grid=self.grid,
                         ordering='C: (i*n2+j)*n3+k; k fastest', units='electrons per Bohr^3',
                         normalization=dict(expected_electrons=volume, integral_electrons=volume),
                         sha256=hashlib.sha256(raw).hexdigest(), volume_bohr3=volume,
                         lattice_columns_bohr=lattice)
        self.save_meta()
        (self.ours / 'results.json').write_text(json.dumps(dict(
            cutoff=dict(ecutwfc_ha=.2), basis=dict(fft_grid=self.grid, n_density_g=7),
            cell=dict(lattice_columns_bohr=lattice))))
        self.qe = self.dir / 'qe.dat'
        write_qe(self.qe, self.values)
        self.counter = 0

    def save_meta(self):
        (self.ours / 'density.json').write_text(json.dumps(self.meta))

    def run_tool(self, values=None, extra=(), ours=None, summary=None):
        if values is not None:
            self.qe.unlink()
            write_qe(self.qe, values)
        self.counter += 1
        report = self.dir / f'report{self.counter}.json'
        md = report.with_suffix('.md')
        args = [sys.executable, '-I', str(ROOT/'scripts/compare_rhog.py'),
                '--ours-dir', str(ours or self.ours),
                '--qe-summary' if summary else '--qe-file', str(summary or self.qe),
                '--json-report', str(report), '--markdown-report', str(md), *extra]
        proc = subprocess.run(args, capture_output=True, text=True)
        result = json.loads(report.read_text()) if report.exists() else None
        if result:
            self.assertIn('Fourier density comparison', md.read_text())
        return proc, result

    def test_identical_passes(self):
        proc, report = self.run_tool()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertLess(report['max_coefficient_difference']['abs_difference'], 1e-12)
        self.assertEqual(report['set_comparison']['only_ours_count'], 0)
        self.assertEqual(report['set_comparison']['only_qe_count'], 0)
        self.assertLess(report['relative_l2_difference'], 1e-12)

    def test_above_tolerance_names_triplet(self):
        values = dict(self.values)
        values[(1, 0, 0)] += 2e-7j
        proc, report = self.run_tool(values)
        self.assertEqual(proc.returncode, 5, proc.stderr)
        self.assertEqual(report['max_coefficient_difference']['miller'], [1, 0, 0])
        self.assertIn('(1, 0, 0)', ' '.join(report['failures']))

    def test_below_tolerance(self):
        values = dict(self.values)
        values[(1, 0, 0)] += 1e-9j
        proc, _ = self.run_tool(values)
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_individual_aggregate_gates(self):
        for gate in ['rho0', 'parseval', 'relative_l2']:
            with self.subTest(gate=gate):
                tolerances = {key: 1.0 for key in DEFAULTS}
                tolerances[gate] = DEFAULTS[gate]
                values = dict(self.values)
                values[(0, 0, 0)] += 1e-5
                report = compare(self.values, values, tolerances)
                self.assertEqual(report['exit_code'], 5)
                self.assertEqual(len(report['failures']), 1)

    def test_summary_maximum_magnitude_gate(self):
        summary = self.summary()
        summary['max_abs_nonzero']['value'] += 2e-7
        report = compare(self.values, self.values, DEFAULTS, summary)
        self.assertEqual(report['exit_code'], 5)
        self.assertEqual(report['failures'], ['maximum nonzero magnitude tolerance exceeded'])

    def test_set_mismatch(self):
        values = dict(self.values)
        del values[(1, 0, 0)]
        values[(1, 1, 0)] = 0j
        proc, report = self.run_tool(values)
        self.assertEqual(proc.returncode, 3, proc.stderr)
        self.assertEqual(report['set_comparison']['only_ours'], [[1, 0, 0]])
        self.assertEqual(report['set_comparison']['only_qe'], [[1, 1, 0]])

    def test_bad_hash(self):
        (self.ours/'density.bin').write_bytes(b'bad')
        proc, _ = self.run_tool()
        self.assertEqual(proc.returncode, 3)
        self.assertIn('SHA-256', proc.stderr)

    def test_bad_integral(self):
        self.meta['normalization']['expected_electrons'] += 1e-8
        self.save_meta()
        proc, _ = self.run_tool()
        self.assertEqual(proc.returncode, 3)
        self.assertIn('integral', proc.stderr)

    def test_bad_metadata(self):
        for key, value in [('dtype', 'f32'), ('endianness', 'big'), ('ordering', 'F'),
                           ('units', 'electrons'), ('fft_grid', [0, 5, 5])]:
            original = self.meta[key]
            with self.subTest(key=key):
                self.meta[key] = value
                self.save_meta()
                proc, _ = self.run_tool()
                self.assertEqual(proc.returncode, 3, proc.stderr)
            self.meta[key] = original

    def test_unrepresentable_qe(self):
        values = dict(self.values)
        values[(3, 0, 0)] = 0j
        proc, _ = self.run_tool(values)
        self.assertEqual(proc.returncode, 3)
        self.assertIn('aliasing', proc.stderr)

    def test_overwrite_refused(self):
        existing = self.dir/'existing.json'
        existing.write_text('evidence')
        proc, _ = self.run_tool(extra=['--json-report', str(existing)])
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(existing.read_text(), 'evidence')
        self.assertFalse((self.dir/'report1.md').exists())

    def test_override_flagged(self):
        proc, report = self.run_tool(extra=['--coefficient-tol', '2e-7'])
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(report['non_default_tolerances'], {'coefficient': 2e-7})
        self.assertEqual(report['override_flags_supplied'], ['--coefficient-tol'])

    def test_invalid_tolerances(self):
        for value in ['nan', '-1', 'inf']:
            proc, _ = self.run_tool(extra=['--coefficient-tol', value])
            self.assertEqual(proc.returncode, 3)

    def summary(self):
        return dict(gamma_only=0, nspin=1, ngm_g=7, rho_000=dict(re=1., im=0.),
                    sum_abs2=sum(abs(z)**2 for z in self.values.values()),
                    max_abs_nonzero=dict(value=.05, miller=[-1, 0, 0]),
                    samples=[dict(miller=m, re=z.real, im=z.imag) for m, z in self.values.items()])

    def test_summary_symmetry_max_does_not_gate_triplet(self):
        summary = self.summary()
        report = compare(self.values, self.values, DEFAULTS, summary)
        self.assertTrue(report['passed'])
        self.assertNotEqual(report['max_abs_nonzero']['ours_miller'],
                            report['max_abs_nonzero']['qe_miller'])
        self.assertAlmostEqual(report['symmetry_noise']['qe_spread'], math.hypot(.05, .03))

    def test_empty_summary_rejected(self):
        summary = self.summary()
        summary['samples'] = []
        path = self.dir/'summary.json'
        path.write_text(json.dumps(summary))
        proc, _ = self.run_tool(summary=path)
        self.assertEqual(proc.returncode, 3)
        self.assertIn('non-empty', proc.stderr)

    def test_usage_and_io(self):
        proc = subprocess.run([sys.executable, '-I', str(ROOT/'scripts/compare_rhog.py')],
                              capture_output=True)
        self.assertEqual(proc.returncode, 1)
        self.qe.unlink()
        proc, _ = self.run_tool()
        self.assertEqual(proc.returncode, 1)

    def test_real_evidence_regression(self):
        evidence = ROOT/'tests/reference/si30_ours'
        proc, report = self.run_tool(ours=evidence/'ours', summary=evidence/'density/qe_rhoG_summary.json')
        self.assertEqual(proc.returncode, 0, proc.stderr+proc.stdout)
        self.assertLess(report['rho0']['abs_difference'], 1e-12)
        for source in ['ours', 'qe']:
            self.assertAlmostEqual(report['rho0'][source]['re'], .0296178744919464, delta=1e-12)
        self.assertLess(report['parseval']['relative_difference'], 1e-6)
        self.assertLess(report['max_coefficient_difference']['abs_difference'], 1e-7)
        self.assertEqual(report['non_default_tolerances'], {})
        self.assertFalse(report['set_comparison']['verified'])
        print('\nREAL QE regression: rho0 difference = %.17g e/Bohr^3; '
              'sum_abs2 relative difference = %.17g; largest sample difference = %.17g e/Bohr^3; '
              'QE symmetry spread = %.17g e/Bohr^3; ours symmetry spread = %.17g e/Bohr^3'
              % (report['rho0']['abs_difference'], report['parseval']['relative_difference'],
                 report['max_coefficient_difference']['abs_difference'],
                 report['symmetry_noise']['qe_spread'], report['symmetry_noise']['ours_spread']))


if __name__ == '__main__':
    unittest.main()
