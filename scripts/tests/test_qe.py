"""FORMAT FIXTURES and synthetic self-tests only; no external reference data."""
import contextlib
import io
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from qe_common import Error, EV_PER_HA, UNVERIFIED_SHA, read_json, sha256, write_json
from make_qe_input import generate, main as generator_main
from verify_qe import parse_output, parse_input, compare, main as verify_main, markdown
from convergence_scan import edit_table, minimum_grid

FIXTURES = ROOT/'tests/fixtures_format'


def synthetic_output(r, perturbation=0):
    """Construct text from our own invented numbers, never from a pw.x run."""
    terms = r['energies']
    ev = ' '.join(format(x*EV_PER_HA, '.17g') for x in r['eigenvalues']['ha'])
    grid = ', '.join(map(str, r['basis']['fft_grid']))
    return f'''FORMAT FIXTURE: synthetic self-test, not a reference
Program PWSCF v.TEST starts on invented date
number of electrons = {r['electrons']}
number of Kohn-Sham states = {r['bands']}
kinetic-energy cutoff = {r['cutoff']['ecutwfc_ry']} Ry
charge density cutoff = {4*r['cutoff']['ecutwfc_ry']} Ry
number of k points = 1
Dense grid: 1111 G-vectors FFT dimensions: ( {grid} )
Exchange-correlation = SLA PZ NOGX NOGC ( 1 1 0 0 )
convergence has been achieved in 11 iterations
End of self-consistent calculation
k = 0.0000 0.0000 0.0000 ( {r['basis']['n_plane_waves']} PWs) bands (ev):

{ev}

highest occupied, lowest unoccupied levels (ev): {r['eigenvalues']['ha'][0]*EV_PER_HA:.17g} {r['eigenvalues']['ha'][1]*EV_PER_HA:.17g}
! total energy = {2*(terms['total']['ha']+perturbation):.17g} Ry
one-electron contribution = {2*sum(terms[k]['ha'] for k in ('kinetic','local','nonlocal')):.17g} Ry
hartree contribution = {2*terms['hartree']['ha']:.17g} Ry
xc contribution = {2*terms['xc']['ha']:.17g} Ry
ewald contribution = {2*terms['ewald']['ha']:.17g} Ry
PWSCF : 1.11s CPU 2.22s WALL
JOB DONE.
'''


def provenance_for(path, r):
    return dict(command='FORMAT FIXTURE ONLY; no command was run',
                qe_version='Program PWSCF v.TEST starts on invented date', host='invented',
                mpi_processes=1, omp_threads=1, date_utc='2026-01-01T00:00:00Z',
                input_sha256=sha256(path), upf_sha256=r['species'][0]['sha256'],
                wall_time_s=1111.1111, returned_by='test generator',
                notes='synthetic self-test, not a reference', is_format_fixture=True)


class ToolTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ours = self.root/'results.json'
        self.r = read_json(FIXTURES/'results.json')
        write_json(self.ours, self.r)
        self.bundle = self.root/'bundle'
        self.bundle.mkdir()
        self.inp = self.bundle/'qe.in'
        generate(self.ours, '/synthetic/pseudo', 'synthetic.UPF', self.inp)
        self.output = synthetic_output(self.r)
        (self.bundle/'qe.out').write_text(self.output)
        self.p = provenance_for(self.inp, self.r)
        self.save_provenance()

    def save_provenance(self):
        write_json(self.bundle/'provenance.json', self.p)

    def invoke(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
            code = verify_main(['--ours', str(self.ours), '--bundle', str(self.bundle),
                                '--json', str(self.root/'report.json'), '--markdown', str(self.root/'report.md')])
        return code, stderr.getvalue()

    def test_golden_input(self):
        self.assertEqual(self.inp.read_text(), (FIXTURES/'golden.in').read_text())
        self.assertNotIn('input_dft', self.inp.read_text())
        meta = read_json(str(self.inp)+'.meta.json')
        self.assertEqual(meta['source_results_sha256'], sha256(self.ours))
        self.assertEqual(meta['input_sha256'], sha256(self.inp))

    def test_rows_columns_and_bit_exact_round_trip(self):
        a = [[math.nextafter(4.0, 5.0), 0.12345678901234567, -0.31],
             [0.5, 5.987654321098765, 0.9876543210987654], [0.875, 1.25, 6.0]]
        self.r['cell']['lattice_columns_bohr'] = a
        write_json(self.ours, self.r)
        generate(self.ours, '/unused', 'synthetic.UPF', self.inp)
        rows = parse_input(self.inp.read_text())[2]
        self.assertEqual(rows, [list(x) for x in zip(*a)])
        self.assertNotEqual(rows, a)

    def test_auto_grid_and_override(self):
        generate(self.ours, '/unused', 'synthetic.UPF', self.inp, fft_grid='auto', ecut_ry=25)
        text = self.inp.read_text()
        self.assertNotIn('nr1', text)
        self.assertIn('ecutwfc = 25,', text)
        self.assertIn('ecutrho = 100,', text)

    def test_generator_isolated_cli(self):
        p = subprocess.run([sys.executable, '-I', str(ROOT/'scripts/make_qe_input.py'),
                            '--results', str(self.ours), '--pseudo-dir', '/no/upf/opened',
                            '--upf-name', 'synthetic.UPF', '--out', str(self.inp)], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)

    def test_generator_rejects_injected_text(self):
        with self.assertRaisesRegex(Error, 'QE text'):
            generate(self.ours, '/unused', 'synthetic.UPF', self.inp, prefix="x'\nnspin=2")

    def test_format_fixture_refusal(self):
        self.assertEqual(self.invoke(), (4, 'FORMAT FIXTURE: not a reference\n'))
        self.assertFalse((self.root/'report.json').exists())

    def test_static_fixture_refusal(self):
        with self.assertRaisesRegex(Error, 'FORMAT FIXTURE'):
            compare(self.ours, FIXTURES/'bundle')

    def test_provenance_modified_input(self):
        # False markers only exist transiently within synthetic self-tests.
        self.p['is_format_fixture'] = False
        self.save_provenance()
        self.inp.write_text(self.inp.read_text()+'! changed\n')
        code, message = self.invoke()
        self.assertEqual(code, 2)
        self.assertIn('input_sha256', message)

    def test_provenance_wrong_upf(self):
        self.p.update(is_format_fixture=False, upf_sha256='b'*64)
        self.save_provenance()
        self.assertEqual(self.invoke()[0], 2)

    def test_provenance_missing_key(self):
        self.p['is_format_fixture'] = False
        del self.p['command']
        self.save_provenance()
        self.assertEqual(self.invoke()[0], 2)

    def test_provenance_wrong_version(self):
        self.p.update(is_format_fixture=False, qe_version='invented different version')
        self.save_provenance()
        self.assertEqual(self.invoke()[0], 2)

    def test_synthetic_self_test_not_a_reference_pass_and_fail(self):
        self.p['is_format_fixture'] = False
        self.save_provenance()
        self.assertEqual(self.invoke()[0], 0)
        report = read_json(self.root/'report.json')
        self.assertIn('synthetic self-test, not a reference', report['header'])
        self.assertIn('synthetic self-test, not a reference', (self.root/'report.md').read_text())
        (self.bundle/'qe.out').write_text(synthetic_output(self.r, 0.01))
        self.assertEqual(self.invoke()[0], 5)
        report = read_json(self.root/'report.json')
        self.assertFalse(report['checks'][0]['passed'])
        self.assertAlmostEqual(report['checks'][0]['absolute_difference'], 0.01, places=9)

    def test_synthetic_self_test_not_a_reference_auto_grid_warning(self):
        generate(self.ours, '/unused', 'synthetic.UPF', self.inp, fft_grid='auto')
        self.p.update(is_format_fixture=False, input_sha256=sha256(self.inp))
        self.save_provenance()
        (self.bundle/'qe.out').write_text(self.output.replace('12, 15, 18', '16, 18, 20'))
        report = compare(self.ours, self.bundle)
        self.assertTrue(report['passed'])
        self.assertTrue(report['warnings'])
        self.assertEqual(report['qe_grid'], [16,18,20])

    def test_synthetic_self_test_not_a_reference_tolerances_and_note(self):
        self.p.update(is_format_fixture=False, upf_sha256=UNVERIFIED_SHA)
        self.r['species'][0]['sha256'] = UNVERIFIED_SHA
        write_json(self.ours, self.r)
        self.save_provenance()
        report = compare(self.ours, self.bundle, {'energy': 1e-5}, 'synthetic convention demonstration')
        self.assertIn('Pseudopotential provenance is unverified', markdown(report))
        self.assertEqual(report['nondefault_tolerances'], {'energy':1e-5})
        self.assertTrue(report['checks'][1]['informational'])

    def test_synthetic_self_test_not_a_reference_eigenvalue_failure(self):
        self.p['is_format_fixture'] = False
        self.save_provenance()
        perturbed = read_json(self.ours)
        perturbed['eigenvalues']['ha'][2] += 0.01
        (self.bundle/'qe.out').write_text(synthetic_output(perturbed))
        self.assertEqual(self.invoke()[0], 5)
        report = read_json(self.root/'report.json')
        self.assertEqual(report['maximum_differences']['eigenvalues']['item'], 'eigenvalue band 3')
        self.assertAlmostEqual(report['maximum_differences']['eigenvalues']['absolute_difference'], 0.01, places=9)

    def test_condition_mismatch_geometry_and_input_cutoff(self):
        self.inp.write_text(self.inp.read_text().replace('4 0.5 0.875', '4.001 0.5 0.875')
                           .replace('Xx 0.125 0.25 0.5', 'Xx 0.125 0.25 0.51')
                           .replace('ecutwfc = 20,', 'ecutwfc = 21,'))
        self.p.update(is_format_fixture=False, input_sha256=sha256(self.inp))
        self.save_provenance()
        code, message = self.invoke()
        self.assertEqual(code, 3)
        for expected in ('lattice', 'atomic position', 'input ecutwfc'):
            self.assertIn(expected, message)

    def test_nonconverged_ours_refused(self):
        self.p['is_format_fixture'] = False
        self.save_provenance()
        self.r['convergence']['converged'] = False
        write_json(self.ours, self.r)
        self.assertEqual(self.invoke()[0], 3)

    def test_provenance_types_and_all_required_keys(self):
        for key in list(self.p):
            with self.subTest(missing=key):
                p = dict(self.p, is_format_fixture=False)
                del p[key]
                write_json(self.bundle/'provenance.json', p)
                self.assertEqual(self.invoke()[0], 2)
        for key, value in [('is_format_fixture', 'false'), ('mpi_processes', 0),
                           ('omp_threads', 2), ('wall_time_s', -1),
                           ('date_utc', 'not a date'), ('date_utc', '2026-01-01T00:00:00+01:00')]:
            with self.subTest(invalid=key, value=value):
                p = dict(self.p, is_format_fixture=False)
                p[key] = value
                write_json(self.bundle/'provenance.json', p)
                self.assertEqual(self.invoke()[0], 2)

    def test_optional_malformed_levels_refused(self):
        text = self.output.replace('highest occupied, lowest unoccupied levels (ev):',
                                   'highest occupied, lowest unoccupied levels (unknown):')
        with self.assertRaisesRegex(Error, 'expected highest occupied'):
            parse_output(text)

    def test_optional_malformed_wall_refused(self):
        with self.assertRaisesRegex(Error, 'expected PWSCF'):
            parse_output(self.output.replace('2.22s WALL', 'unknown WALL'))
        self.assertEqual(parse_output(self.output)['wall_time_s'], 2.22)

    def test_units(self):
        q = parse_output(self.output)
        self.assertEqual(q['total'], self.r['energies']['total']['ha'])
        for ours, parsed in zip(self.r['eigenvalues']['ha'], q['eigenvalues_ha']):
            self.assertAlmostEqual(ours, parsed, places=11)
        self.assertEqual(EV_PER_HA, 27.211386245988)

    def test_optional_lines_absent(self):
        output = '\n'.join(x for x in self.output.splitlines() if not x.startswith(('highest occupied', 'PWSCF :')))
        q = parse_output(output)
        self.assertIsNone(q['levels_ev'])
        self.assertIsNone(q['wall_time_text'])

    def test_bad_eigenvalue_token(self):
        with self.assertRaisesRegex(Error, 'eigenvalues'):
            parse_output(self.output.replace(format(self.r['eigenvalues']['ha'][0]*EV_PER_HA, '.17g'), 'bad', 1))

    def test_multiple_runs_rejected(self):
        with self.assertRaisesRegex(Error, 'version'):
            parse_output(self.output+self.output)

    def test_unexpected_input_parameter(self):
        with self.assertRaisesRegex(Error, 'supported SYSTEM'):
            parse_input(self.inp.read_text().replace('  ibrav = 0,', '  ibrav = 0,\n  nspin = 2,'))

    def test_usage_exit_one(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(verify_main([]), 1)
            self.assertEqual(generator_main([]), 1)


# One discoverable test per required output record; no broad catch-all success.
MISSING = {
    'version': 'Program PWSCF', 'electrons': 'number of electrons',
    'bands': 'number of Kohn-Sham', 'ecutwfc': 'kinetic-energy cutoff',
    'ecutrho': 'charge density cutoff', 'nk': 'number of k points',
    'grid': 'Dense grid', 'xc': 'Exchange-correlation',
    'convergence': 'convergence has been achieved', 'final_block': 'End of self-consistent',
    'plane_waves': 'k = ', 'total': '! total energy', 'one_electron': 'one-electron contribution',
    'hartree': 'hartree contribution', 'xc_energy': 'xc contribution', 'ewald': 'ewald contribution',
    'done': 'JOB DONE.',
}
for name, prefix in MISSING.items():
    def test(self, prefix=prefix):
        stripped = '\n'.join(line for line in self.output.splitlines() if not line.startswith(prefix))
        with self.assertRaisesRegex(Error, 'expected'):
            parse_output(stripped)
    setattr(ToolTests, 'test_missing_pattern_' + name, test)


def missing_eigenvalues(self):
    line = ' '.join(format(x*EV_PER_HA, '.17g') for x in self.r['eigenvalues']['ha'])
    with self.assertRaisesRegex(Error, 'eigenvalues'):
        parse_output(self.output.replace(line, ''))
ToolTests.test_missing_pattern_eigenvalues = missing_eigenvalues

for name, before, after in [
    ('cutoff', 'kinetic-energy cutoff = 20.0', 'kinetic-energy cutoff = 25.0'),
    ('grid', '12, 15, 18', '12, 16, 18'),
    ('electrons', 'number of electrons = 2.0', 'number of electrons = 4.0'),
    ('functional', 'SLA PZ NOGX NOGC ( 1 1 0 0 )', 'SLA PW NOGX NOGC'),
    ('npw', '( 1111 PWs)', '( 1112 PWs)'),
]:
    def test(self, before=before, after=after):
        self.p['is_format_fixture'] = False
        self.save_provenance()
        (self.bundle/'qe.out').write_text(self.output.replace(before, after))
        code, message = self.invoke()
        self.assertEqual(code, 3, message)
        self.assertIn('conditions mismatch', message)
    setattr(ToolTests, 'test_condition_mismatch_' + name, test)


if __name__ == '__main__':
    unittest.main()
