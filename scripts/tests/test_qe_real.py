"""Regression against immutable stdout returned by the remote QE 7.6 runner."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
BUNDLE = ROOT / 'tests/reference/si30_ours'
LEVEL_LABEL = 'highest occupied, lowest unoccupied level (ev):'


class RealBundleTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.temp = Path(temporary.name)
        self.bundle = self.temp / 'bundle'
        shutil.copytree(BUNDLE, self.bundle)
        self.output = (BUNDLE / 'qe.out').read_text()
        self.report_number = 0

    def invoke(self, bundle=None):
        bundle = self.bundle if bundle is None else bundle
        self.report_number += 1
        report = self.temp / f'report-{self.report_number}.json'
        markdown = report.with_suffix('.md')
        process = subprocess.run(
            [sys.executable, '-I', str(ROOT / 'scripts/verify_qe.py'),
             '--bundle', str(bundle), '--ours', str(bundle / 'ours/results.json'),
             '--json', str(report), '--markdown', str(markdown)],
            cwd=self.temp, capture_output=True, text=True)
        return process, report, markdown

    def assert_refused(self, code, message):
        process, report, markdown = self.invoke()
        self.assertEqual(process.returncode, code, process.stderr)
        self.assertIn(message, process.stderr)
        self.assertFalse(report.exists())
        self.assertFalse(markdown.exists())
        return process.stderr

    def edit_json(self, relative, edit):
        path = self.bundle / relative
        data = json.loads(path.read_text())
        edit(data)
        path.write_text(json.dumps(data))

    def test_real_bundle_full_comparison(self):
        # Exercise the repository evidence directly; destinations and cwd are temporary.
        process, report_path, markdown_path = self.invoke(BUNDLE)
        self.assertEqual(process.returncode, 0, process.stderr)
        report = json.loads(report_path.read_text())
        q = report['qe']
        lines = self.output.splitlines()

        def record(prefix):
            found = [line.strip() for line in lines if line.strip().startswith(prefix)]
            self.assertEqual(len(found), 1, prefix)
            return found[0]

        # Independent oracle: literal line selection and token splitting, no imports
        # from the verifier or its regex/number/unit helpers.
        for key, prefix, headline in [
            ('total', '!    total energy', -14.59106606),
            ('one_electron_ha', 'one-electron contribution', 5.58356185),
            ('hartree_ha', 'hartree contribution', 1.67209213),
            ('xc_ha', 'xc contribution', -5.04777627),
            ('ewald_ha', 'ewald contribution', -16.79894377),
        ]:
            printed = float(record(prefix).split('=')[1].split()[0])
            self.assertEqual(printed, headline)
            self.assertEqual(q[key], printed / 2)
        for key, prefix, headline in [
            ('electrons', 'number of electrons', 8),
            ('bands', 'number of Kohn-Sham states', 8),
        ]:
            printed = float(record(prefix).split('=')[1])
            self.assertEqual(printed, headline)
            self.assertEqual(q[key], printed)
        iterations = int(record('convergence has been achieved').split()[-2])
        self.assertEqual(iterations, 11)
        self.assertEqual(q['iterations'], iterations)
        dense = record('Dense  grid:')
        count = int(dense.split(':')[1].split()[0])
        grid = [int(x) for x in dense.split('(')[1].split(')')[0].split(',')]
        self.assertEqual(count, 5985)
        self.assertEqual(q['dense_g_vectors'], count)
        self.assertEqual(grid, [25, 25, 25])
        self.assertEqual(q['grid'], grid)
        final = self.output.split('End of self-consistent calculation')[1]
        band_header, band_body = final.split('bands (ev):', 1)
        npw = int(band_header.split('(')[1].split()[0])
        self.assertEqual(npw, 749)
        self.assertEqual(q['npw'], npw)
        eigenvalues = [float(x) for x in band_body.strip().split('\n\n')[0].split()]
        self.assertEqual(eigenvalues, [-5.1733, 7.0378, 7.0378, 7.0378,
                                      9.1742, 9.1742, 9.1742, 10.2981])
        self.assertEqual(q['eigenvalues_ha'], [x / 27.211386245988 for x in eigenvalues])
        levels = [float(x) for x in record(LEVEL_LABEL).split(':')[1].split()]
        self.assertEqual(levels, [7.0378, 9.1742])
        self.assertEqual(q['levels_ev'], levels)
        functional = record('Exchange-correlation=').split('=')[1].strip()
        self.assertEqual(functional, 'SLA  PZ   NOGX NOGC')
        self.assertEqual(q['xc'], functional)
        self.assertTrue(report['passed'])
        self.assertEqual(len(report['checks']), 17)
        for check in report['checks']:
            self.assertIs(check['passed'], True, check['item'])
            self.assertIs(check['informational'], False, check['item'])
        self.assertEqual(report['nondefault_tolerances'], {})
        warning = 'Pseudopotential provenance is unverified'
        self.assertIn(warning, report['header'])
        markdown = markdown_path.read_text()
        self.assertIn(warning, markdown.split('## Provenance')[0])
        rows = [line for line in markdown.splitlines() if line.endswith('| PASS |')]
        self.assertEqual(len(rows), len(report['checks']))

    def test_input_hash_mismatch(self):
        path = self.bundle / 'qe.in'
        original = path.read_text()
        path.write_text(original.replace('&CONTROL', '&CONTROX', 1))
        self.assertNotEqual(path.read_text(), original)
        self.assert_refused(2, 'provenance mismatch: input_sha256')

    def test_upf_hash_mismatch(self):
        self.edit_json('provenance.json', lambda p: p.update(upf_sha256='0' * 64))
        self.assert_refused(2, 'provenance mismatch: upf_sha256')

    def test_fixture_refusal(self):
        self.edit_json('provenance.json', lambda p: p.update(is_format_fixture=True))
        self.assert_refused(4, 'FORMAT FIXTURE: not a reference')

    def test_missing_total_energy(self):
        lines = self.output.splitlines(keepends=True)
        kept = [line for line in lines if not line.startswith('!    total energy')]
        self.assertEqual(len(lines) - len(kept), 1)
        (self.bundle / 'qe.out').write_text(''.join(kept))
        message = self.assert_refused(3, 'expected ! total energy: /')
        self.assertIn('(found 0)', message)

    def test_ours_cutoff_mismatch(self):
        self.edit_json('ours/results.json',
                       lambda r: r['cutoff'].update(ecutwfc_ha=16.0, ecutwfc_ry=32.0))
        message = self.assert_refused(3, 'conditions mismatch:')
        self.assertIn('ecutwfc', message)

    def test_level_variants(self):
        original = next(line for line in self.output.splitlines() if LEVEL_LABEL in line)
        variants = [
            (LEVEL_LABEL + ' 7.0378 9.1742', [7.0378, 9.1742]),
            (LEVEL_LABEL.replace('level (', 'levels (') + ' 7.0378 9.1742', [7.0378, 9.1742]),
            ('highest occupied level (ev): 7.0378', [7.0378]),
        ]
        for line, expected in variants:
            with self.subTest(line=line):
                (self.bundle / 'qe.out').write_text(self.output.replace(original, line))
                process, path, _ = self.invoke()
                self.assertEqual(process.returncode, 0, process.stderr)
                self.assertEqual(json.loads(path.read_text())['qe']['levels_ev'], expected)

    def test_malformed_levels(self):
        original = next(line for line in self.output.splitlines() if LEVEL_LABEL in line)
        for line in [
            LEVEL_LABEL.replace('level (', 'levelz (') + ' 7.0378 9.1742',
            LEVEL_LABEL + ' 7.0378',
            LEVEL_LABEL + ' 7.0378 9.1742 10.2981',
            LEVEL_LABEL.replace('(ev)', '(eV)') + ' 7.0378 9.1742',
            'highest occupied levels (ev): 7.0378',
            'highest occupied level (ev): 7.0378 9.1742',
        ]:
            with self.subTest(line=line):
                (self.bundle / 'qe.out').write_text(self.output.replace(original, line))
                message = self.assert_refused(3, 'expected highest occupied')
                self.assertIn('(found 0)', message)


if __name__ == '__main__':
    unittest.main()
