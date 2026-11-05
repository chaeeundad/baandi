"""Independent analytic norms and fail-closed equal-accuracy evidence tests."""
import copy
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import equal_accuracy as ea
from qe_common import Error
from verify_qe import require_strict_result
import verify_qe
import compare_materials_qe


class EqualAccuracyTests(unittest.TestCase):
    def test_validation_rejects_standard_at_each_provenance_location(self):
        for key in (None, 'scf', 'convergence'):
            for label in ('standard', 'typo', None):
                result = {} if key is None else {key: {}}
                section = result if key is None else result[key]
                section['convergence_profile'] = label
                with self.assertRaises(Error):
                    require_strict_result(result)
            result = {'certified': False} if key is None else {key: {'certified': False}}
            with self.assertRaises(Error):
                require_strict_result(result)
        require_strict_result({})
        require_strict_result({'convergence_profile': 'strict', 'certified': True})

    def test_validation_entrypoints_reject_before_reading_qe(self):
        for mod, loader in ((verify_qe, 'results'), (compare_materials_qe, 'materials_results')):
            with patch.object(mod, loader, return_value=(Path('unused'), {'convergence_profile': 'standard'})):
                with self.assertRaisesRegex(Error, 'strict convergence_profile'):
                    mod.compare('unused', '/nonexistent/qe')

    def test_analytic_density_metric_and_zero_mode(self):
        side, amplitude = 8., .03
        grid = [8, 1, 1]
        meta = dict(fft_grid=grid, lattice_columns_bohr=[[side,0,0],[0,side,0],[0,0,side]], volume_bohr3=side**3)
        coeff = {(0,0,0): .2, (1,0,0): amplitude/2, (-1,0,0): amplitude/2}
        zero = {k: 0j for k in coeff}
        real = [.2 + amplitude*math.cos(2*math.pi*i/8) for i in range(8)]
        measured = ea.density_accuracy((coeff, meta, real), (zero, meta, [0.]*8))
        expected = math.pi*side**3*amplitude**2/(2*math.pi/side)**2
        self.assertAlmostEqual(measured['density_hartree_energy_ha'], expected, places=13)
        self.assertAlmostEqual(measured['density_l2_bohr_m3_over_2']**2, side**3*(.2**2+amplitude**2/2), places=12)
        self.assertEqual(ea.density_accuracy((coeff,meta,real), (coeff,meta,real))['density_hartree_energy_ha'], 0.)
        wrong = copy.deepcopy(meta)
        wrong['fft_grid'][0] = 9
        with self.assertRaises(ValueError):
            ea.density_accuracy((coeff,meta,real), (coeff,wrong,real))

    def test_accuracy_includes_first_empty_and_rejects_k_order(self):
        ref = dict(occupied=2, kpoints=[dict(fractional=[0,0,0],weight=1)],
                   eigenvalues_ha=[[1.,2.,3.]], total_energy_ha=-5.)
        row = copy.deepcopy(ref)
        row['eigenvalues_ha'][0][2] += .5
        row['total_energy_ha'] += .125
        got = ea.accuracy(row, ref)
        self.assertEqual(got['max_occupied_eigenvalue_difference_ha'], 0.)
        self.assertEqual(got['max_first_empty_eigenvalue_difference_ha'], .5)
        self.assertEqual(got['energy_difference_ha'], .125)
        row['kpoints'][0]['fractional'][0] = .5
        with self.assertRaises(ValueError):
            ea.accuracy(row, ref)

    def test_scf_settings_does_not_change_other_tables(self):
        text = '[scf]\nthreads = 1\n[other]\nthreads = 4\n'
        got = ea.scf_settings(text, {'threads': 8, 'convergence_profile': 'standard'})
        self.assertIn('threads = 8', got)
        self.assertTrue(got.endswith('[other]\nthreads = 4\n'))
        self.assertEqual(got.count('threads ='), 2)

    def test_qe_units_rounding_timing_and_missing_work(self):
        fixture = ea.module(ROOT / 'scripts/tests/test_materials_qe.py')
        r = json.loads((ROOT / 'scripts/tests/fixtures/materials/si-20.json').read_text())
        text = fixture.output(r)
        import re
        text = re.sub(r'! total energy = .* Ry', '! total energy = -10.12345678 Ry', text)
        text += ' PWSCF : 1.2s CPU 1m 2.50s WALL\n'
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'qe.out'
            path.write_text(text)
            q = ea.qe_row(path)
            self.assertEqual(q['total_energy_ha'], -10.12345678/2)
            self.assertEqual(q['energy_rounding_bound_ha'], 2.5e-9)
            self.assertEqual(q['wall_s'], 62.5)
            self.assertIsNone(q['h_vectors'])
            self.assertIsNone(q['davidson_steps'])
            path.write_text(text.replace('JOB DONE.', ''))
            with self.assertRaises(Error):
                ea.qe_row(path)

    def test_failed_reference_has_no_accuracy(self):
        with tempfile.TemporaryDirectory() as temp:
            row, density = ea.run_one(Path('/nonexistent/cli'), '[scf]\n', Path(temp)/'run', 'tight', 'strict', 1.)
            self.assertEqual(row['status'], 'failed')
            self.assertIsNone(row['accuracy'])
            self.assertIsNone(density)

    def test_qe_pair_provenance_threshold_and_own_reference(self):
        fixture = ea.module(ROOT / 'scripts/tests/test_materials_qe.py')
        r = json.loads((ROOT / 'scripts/tests/fixtures/materials/si-20.json').read_text())
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            entry = dict(case='synthetic', upf_sha256={'Si': 'a'*64})
            for key, threshold, energy in [('standard', '1d-8', -5.), ('tight', '1d-14', -5.125)]:
                (directory / (key+'.in')).write_text('&electrons\nconv_thr='+threshold+'\n/\n')
                r['total_energy']['ha'] = energy
                (directory / (key+'.out')).write_text(fixture.output(r))
                entry[key] = dict(input=key+'.in', output=key+'.out', command=['pw.x', '-in', key+'.in'])
            manifest = directory / 'pairs.json'
            manifest.write_text(json.dumps([entry]))
            rows = ea.qe_pairs(manifest)
            self.assertEqual(rows[0]['accuracy']['energy_difference_ha'], .125)
            self.assertFalse(rows[0]['reference_check']['verified'])
            self.assertEqual(rows[1]['accuracy']['energy_difference_ha'], 0.)
            (directory / 'tight.in').write_text('&electrons\nconv_thr=1d-12\n/\n')
            with self.assertRaisesRegex(ValueError, 'conv_thr mismatch'):
                ea.qe_pairs(manifest)

    def test_standard_without_false_certification_rejected(self):
        base = dict(convergence_profile='standard', convergence=dict(converged=True, history=[{}]))
        with self.assertRaisesRegex(ValueError, 'explicitly be uncertified'):
            ea.baandi_measurement(base, 'standard')

    def test_stored_qe_outputs_are_read_without_new_remote_runs(self):
        paths = list((ROOT/'tests/reference/kpoints').glob('*/qe.out')) + list((ROOT/'tests/reference').glob('si30_*/qe.out'))
        self.assertEqual(len(paths), 4)
        for path in paths:
            row = ea.qe_row(path)
            self.assertTrue(row['final_estimated_scf_error_ha'] > 0)
            self.assertTrue(row['davidson_steps'] > 0)
            self.assertTrue(row['h_blocks'] > 0)
            self.assertIsNone(row['h_vectors'])

    def test_optional_qe_density_metric_from_stdout_geometry(self):
        from qe_rhog import write_qe
        side, amplitude = 8., .03
        text = 'lattice parameter (alat) = 8.0 a.u.\n'
        text += 'b(1) = (1 0 0)\nb(2) = (0 1 0)\nb(3) = (0 0 1)\n'
        with tempfile.TemporaryDirectory() as temp:
            a,b = Path(temp)/'a.dat',Path(temp)/'b.dat'
            coeff = {(0,0,0):0.2, (1,0,0):amplitude/2, (-1,0,0):amplitude/2}
            write_qe(a,coeff);write_qe(b,{k:0j for k in coeff})
            got = ea.qe_density_accuracy(a,b,text)
            expected = math.pi*side**3*amplitude**2/(2*math.pi/side)**2
            self.assertAlmostEqual(got['density_hartree_energy_ha'],expected,places=12)


if __name__ == '__main__':
    unittest.main()
