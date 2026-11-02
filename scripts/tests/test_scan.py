"""Real rust-dft runs on a temporary synthetic UPF; never QE references."""
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from convergence_scan import edit_table, minimum_grid, smooth_at_least
from qe_common import Error, read_json, sha256


def synthetic_upf():
    n = 401
    r = [i*0.025 for i in range(n)]
    def array(values):
        return ' '.join(format(x, '.17g') for x in values)
    return f'''<UPF version="2.0.1"><PP_INFO>SYNTHETIC FORMAT TEST ONLY; NOT A REFERENCE</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx"
pseudo_type="NC" relativistic="nonrelativistic" is_ultrasoft="F" is_paw="F" is_coulomb="F"
has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0"
l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="1"/>
<PP_MESH dx="0.025" mesh="{n}" xmin="0" rmax="10" zmesh="1"><PP_R>{array(r)}</PP_R><PP_RAB>{array([0.025]*n)}</PP_RAB></PP_MESH>
<PP_LOCAL>{array([-8/math.sqrt(x*x+1) for x in r])}</PP_LOCAL>
<PP_NONLOCAL><PP_BETA.1 index="1" angular_momentum="0" cutoff_radius="10" cutoff_radius_index="401">{array([x*math.exp(-x*x) for x in r])}</PP_BETA.1><PP_DIJ>2</PP_DIJ></PP_NONLOCAL>
<PP_RHOATOM>{array([0]*n)}</PP_RHOATOM></UPF>'''


INPUT = '''[system]
name = 'synthetic scan, not a reference'
lattice_unit = 'bohr'
lattice_vectors = [[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]
[[system.atoms]]
species = 'Xx'
fractional = [0.13,0.21,0.32]
[pseudopotentials]
Xx = 'synthetic.xml'
[basis]
ecutwfc_ry = 3.0
kpoints = [[0.0,0.0,0.0]]
[xc]
functional = 'lda-pz'
'''


class ScanTests(unittest.TestCase):
    def test_text_edit_preserves_unrelated_settings(self):
        changed = edit_table(INPUT, 'basis', {'ecutwfc_ry':'4.0', 'fft_grid':'[5,5,5]'})
        self.assertIn('ecutwfc_ry = 4.0', changed)
        self.assertIn('fft_grid = [5,5,5]', changed)
        self.assertIn('kpoints = [[0.0,0.0,0.0]]', changed)
        self.assertEqual(changed.split('[basis]')[0], INPUT.split('[basis]')[0])
        self.assertNotIn('fft_grid', edit_table(changed, 'basis', {'fft_grid':None}))
        self.assertIn('[scf]\nenergy_tolerance_ha = 1e-8', edit_table(INPUT, 'scf', {'energy_tolerance_ha':'1e-8'}))

    def test_multiline_edit_refused(self):
        with self.assertRaisesRegex(Error, 'single-line'):
            edit_table('[basis]\nfft_grid = [\n3,3,3]\n', 'basis', {'fft_grid':None})

    def test_grid_rules(self):
        self.assertEqual(smooth_at_least(7), 8)
        self.assertEqual(smooth_at_least(17), 18)
        r = {'basis':{'fft_grid':[5,5,5]}, 'cutoff':{'ecutwfc_ry':3},
             'cell':{'reciprocal_columns_inverse_bohr':[[math.pi/2,0,0],[0,math.pi/2,0],[0,0,math.pi/2]]}}
        self.assertEqual(minimum_grid(r), [5,5,5])

    def test_real_binary_synthetic_scan_all(self):
        binary = Path(os.environ.get('DFT_CLI', str(ROOT/'target/debug/dft-cli'))).resolve()
        if not binary.is_file():
            self.skipTest('dft-cli not built: run scripts/check.sh or set DFT_CLI; driver receives --dft-cli')
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder/'input.toml'
            inp.write_text(INPUT)
            (folder/'synthetic.xml').write_text(synthetic_upf())
            original = sha256(inp)
            out = folder/'scan'
            command = [sys.executable, '-I', str(ROOT/'scripts/convergence_scan.py'),
                       '--input', str(inp), '--dft-cli', str(binary), '--out-dir', str(out),
                       '--scan', 'all', '--ecuts', '3,4', '--timeout', '30',
                       '--emit-qe-inputs', '--pseudo-dir', '/remote/synthetic', '--upf-name', 'synthetic.xml']
            p = subprocess.run(command, capture_output=True, text=True, timeout=120)
            self.assertEqual(p.returncode, 0, p.stdout+p.stderr)
            self.assertEqual(sha256(inp), original)
            manifest = read_json(out/'manifest.json')
            rows = read_json(out/'scan.json')
            self.assertEqual(manifest['status'], 'done')
            self.assertEqual(len(rows), 11)  # two cutoffs, four grids, five tolerances
            self.assertEqual(len(manifest['qe_inputs']), 22)
            for record in manifest['points']:
                self.assertEqual(record['exit_code'], 0)
                self.assertEqual(record['input_sha256'], sha256(record['input']))
                self.assertTrue(record['commit_hash'])
                self.assertTrue(read_json(record['results'])['convergence']['converged'])
            for record in manifest['qe_inputs']:
                self.assertEqual(record['sha256'], sha256(record['input']))
                self.assertEqual(record['our_results_sha256'], sha256(record['our_results']))
            for series in ('ecut','grid','tolerance'):
                finest = max((r for r in rows if r['scan']==series), key=lambda r:r['rank'])
                self.assertEqual(finest['difference_to_finest_ha'], 0)
            self.assertTrue((out/'scan.csv').is_file())
            self.assertTrue((out/'scan.md').is_file())
            print('\nSynthetic rust-dft scan (not a reference): '
                  f"{len(rows)} certified points; {len(manifest['qe_inputs'])} QE inputs; "
                  f"energies {min(r['total_energy_ha'] for r in rows):.17g} to "
                  f"{max(r['total_energy_ha'] for r in rows):.17g} Ha")
            # Existing destination is rejected instead of reusing stale results.
            p = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(p.returncode, 1)

    def test_failed_run_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            inp = folder/'input.toml'
            inp.write_text(INPUT)
            command = [sys.executable, '-I', str(ROOT/'scripts/convergence_scan.py'),
                       '--input', str(inp), '--dft-cli', '/usr/bin/false', '--out-dir', str(folder/'scan'),
                       '--scan', 'ecut', '--ecuts', '3']
            p = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(p.returncode, 1)
            manifest = read_json(folder/'scan/manifest.json')
            self.assertEqual(manifest['status'], 'failed')
            self.assertEqual(manifest['points'][0]['exit_code'], 1)
            self.assertFalse((folder/'scan/scan.json').exists())


if __name__ == '__main__':
    unittest.main()
