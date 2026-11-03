"""Independent XML checks and process-level suite regression tests (no Rust build)."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('bench_suite', ROOT / 'scripts/bench_suite.py')
suite = importlib.util.module_from_spec(spec)
spec.loader.exec_module(suite)


def fixture_report(converged=True):
    return {'convergence': {'converged': converged, 'history': [
        {'operator_applications': 10, 'fixed_point_certificate': {'operator_applications': 3}},
        {'operator_applications': 7, 'fixed_point_certificate': {'operator_applications': 5}}]},
        'total_energy': {'ha': -2.0}, 'timings': {'phases': {'density_s': .1}, 'engine_total_s': .2},
        'basis': {'fft_grid': [7,11,13]}, 'bands': 3, 'occupations': [2,0,0],
        'eigenvalues': {'ha': [-1,.2,.2]}}


class SuiteTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)
        self.case = suite.load_case(ROOT / 'benchmarks/cases/syn-s1.toml')

    def fake(self, body):
        p = self.path / 'fake-cli'
        p.write_text('#!' + sys.executable + '\nimport sys,json,time\nfrom pathlib import Path\n' + body)
        p.chmod(0o755)
        return p

    def run_case(self, binary, **kwargs):
        return suite.run_case(binary, self.case, None, self.path / 'work', 1, 1,
                              kwargs.get('repeats', 1), kwargs.get('timeout', 2))

    def passing_binary(self, converged=True):
        return self.fake("out=Path(sys.argv[-1]);out.mkdir()\n"
                         f"(out/'results.json').write_text({json.dumps(json.dumps(fixture_report(converged)))})\n"
                         "(out/'density.bin').write_bytes(b'density')\n" +
                         ('' if converged else "sys.exit(3)\n"))

    def test_catalog_and_render(self):
        cases = [suite.load_case(p) for p in (ROOT / 'benchmarks/cases').glob('*.toml')]
        self.assertGreaterEqual(len(cases), 24)
        self.assertEqual(len(cases), len({c['id'] for c in cases}))
        for c in cases:
            text = suite.case_input(c, '/some/path with space.UPF', 7, 4)
            self.assertIn('threads = 4', text)
            self.assertIn('seed = 7', text)
            if c['potential'] == 'synthetic':
                ET.fromstring(suite.synthetic(c))
        for invalid in ('id = "x"\nid = "y"', 'unknown = 1', 'id = NaN'):
            f = self.path / 'bad.toml'; f.write_text(invalid)
            with self.assertRaises(ValueError):
                suite.load_case(f)

    def test_upf_independent_xml_and_units(self):
        for z, channels, shape in [(1,[0],'gaussian'), (4,[0,0,0],'soft'), (14,[0,1,2],'gaussian'), (4,[2,2,2],'gaussian')]:
            root = ET.fromstring(suite._upf.generate(z, channels, .8, shape))
            header = root.find('PP_HEADER')
            self.assertEqual(header.attrib['number_of_proj'], str(len(channels)))
            self.assertEqual(header.attrib['pseudo_type'], 'NC')
            self.assertIn('NON-PHYSICAL', root.find('PP_INFO').text)
            r = list(map(float, root.find('PP_MESH/PP_R').text.split()))
            rab = list(map(float, root.find('PP_MESH/PP_RAB').text.split()))
            local = list(map(float, root.find('PP_LOCAL').text.split()))
            self.assertEqual(len(r), 1201)
            self.assertEqual(len(local), len(r))
            self.assertAlmostEqual(r[-1], 12)
            self.assertTrue(all(y>x for x,y in zip(r,r[1:])))
            self.assertTrue(all(x>0 for x in rab))
            expected = -4*z/(.8*math.sqrt(math.pi)) if shape == 'gaussian' else -2*z/.8
            self.assertAlmostEqual(local[0], expected)
            self.assertAlmostEqual(local[-1], -2*z/12, delta=.003*z)
            dij = list(map(float, root.find('PP_NONLOCAL/PP_DIJ').text.split()))
            n = len(channels)
            self.assertEqual(len(dij), n*n)
            for i,l in enumerate(channels):
                beta = root.find(f'PP_NONLOCAL/PP_BETA.{i+1}')
                values = list(map(float, beta.text.split()))
                self.assertEqual(int(beta.attrib['angular_momentum']), l)
                self.assertEqual(len(values), len(r))
                self.assertAlmostEqual(values[101], r[101]**(l+1)*math.exp(-r[101]**2/(i+1)))
                for j in range(n):
                    self.assertEqual(dij[i*n+j], dij[j*n+i])
                    if channels[i] != channels[j]:
                        self.assertEqual(dij[i*n+j], 0)
        for options in ({'valence':0}, {'valence':15}, {'channels':[0]*4}, {'channels':[3]}, {'radius':float('nan')}, {'points':100}):
            with self.assertRaises(ValueError):
                suite._upf.generate(**options)

    def test_variants_are_reproducible_independent_and_bounded(self):
        original = copy.deepcopy(self.case)
        a = suite.variant(self.case, 7, 0, .01)
        self.assertEqual(a, suite.variant(self.case, 7, 0, .01))
        self.assertNotEqual(a['atoms'], suite.variant(self.case, 8, 0, .01)['atoms'])
        self.assertEqual(original, self.case)
        for site, base in zip(a['atoms'], self.case['atoms']):
            for x,y in zip(site,base):
                self.assertLessEqual(abs((x-y+.5)%1-.5), .01)

    def test_success_counters_hash_and_repeats(self):
        result = self.run_case(self.passing_binary(), repeats=2)
        self.assertEqual(result['status'], 'pass')
        self.assertTrue(result['repeat_identical'])
        self.assertEqual(result['operator_applications'], 25)
        self.assertEqual(result['runs'][0]['minimum_direct_gap_ha'], 1.2)
        self.assertEqual(result['density_sha256'], suite.digest(b'density'))
        self.assertIn('synthetic.UPF', result['input_toml'])
        self.assertNotIn(str(self.path), result['input_toml'])

    def test_error_preserves_message(self):
        self.case['expected_error'] = 'not converged after 400 iterations'
        r = self.run_case(self.fake("sys.stderr.write('not converged after 400 iterations');sys.exit(4)\n"))
        self.assertEqual(r['status'], 'error')
        self.assertTrue(r['expected_failure'])
        self.assertIsNone(r['iterations'])
        self.assertIn(r['input_toml'].strip(), suite.markdown({'cases': {'x':r}, 'comparison':[]}))

    def test_timeout_is_not_expected_solver_failure(self):
        self.case['expected_error'] = 'not converged'
        r = self.run_case(self.fake('time.sleep(10)\n'), timeout=.05)
        self.assertEqual(r['status'], 'timeout')
        self.assertFalse(r['expected_failure'])
        self.assertLess(r['wall_median_s'], 2)

    def test_budget_exhaustion_retains_energy(self):
        r = self.run_case(self.passing_binary(False))
        self.assertEqual(r['status'], 'nonconverged')
        self.assertFalse(r['converged'])
        self.assertEqual(r['total_energy_ha'], -2)

    def test_exit_zero_missing_results_fails(self):
        r = self.run_case(self.fake('pass\n'))
        self.assertEqual(r['status'], 'invalid_output')

    def test_malformed_or_nonfinite_output_fails(self):
        report = fixture_report()
        report['total_energy']['ha'] = float('nan')
        binary = self.fake("out=Path(sys.argv[-1]);out.mkdir()\n"
                           f"(out/'results.json').write_text({json.dumps(json.dumps(report))})\n")
        r = self.run_case(binary)
        self.assertEqual(r['status'], 'invalid_output')
        self.assertIn('nonfinite', r['runs'][0]['error'])

    def test_repeat_nondeterminism_is_visible(self):
        report = json.dumps(json.dumps(fixture_report()))
        binary = self.fake("out=Path(sys.argv[-1]);out.mkdir()\n"
                           f"(out/'results.json').write_text({report})\n"
                           "(out/'density.bin').write_text(out.name)\n")
        r = self.run_case(binary, repeats=2)
        self.assertEqual(r['status'], 'nondeterministic')
        self.assertFalse(r['repeat_identical'])

    def test_fingerprint_is_portable_but_captures_seed(self):
        binary = self.passing_binary()
        first = self.run_case(binary)
        second = suite.run_case(binary, self.case, None, self.path/'elsewhere', 1, 1, 1, 2)
        changed = suite.run_case(binary, self.case, None, self.path/'seed7', 7, 1, 1, 2)
        self.assertEqual(first['fingerprint'], second['fingerprint'])
        self.assertNotEqual(first['input_sha256'], second['input_sha256'])
        self.assertNotEqual(first['fingerprint'], changed['fingerprint'])

    def test_comparison(self):
        r = self.run_case(self.passing_binary())
        base = {'x':r}
        self.assertFalse(suite.compare(base,base)[0]['regressions'])
        for field,value in [('fingerprint','changed'), ('iterations',100), ('total_energy_ha', -2+2e-12), ('converged',False), ('status','timeout')]:
            new = copy.deepcopy(base);new['x'][field]=value
            self.assertTrue(suite.compare(new,base)[0]['regressions'], field)
        new = copy.deepcopy(base);new['x']['density_sha256']='different'
        self.assertFalse(suite.compare(new,base)[0]['regressions'])
        self.assertTrue(suite.compare(new,base,True)[0]['regressions'])
        self.assertTrue(suite.compare(base,{})[0]['regressions'])
        old = copy.deepcopy(base);old['x'].update(status='error',converged=False)
        comp=suite.compare(base,old)[0]
        self.assertFalse(comp['regressions'])
        self.assertTrue(comp['improvements'])

    def test_robustness_cli_continues_and_missing_baseline_errors(self):
        binary = self.fake("sys.stderr.write('solver failed');sys.exit(4)\n")
        output = self.path / 'report.json'
        args = [sys.executable, '-I', str(ROOT / 'scripts/bench_suite.py'), '--binary', str(binary),
                '--cases', 'syn-s1,syn-s2', '--output', str(output), '--work-dir', str(self.path/'runs')]
        p = subprocess.run(args+['--robustness'], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        report = json.loads(output.read_text())
        self.assertEqual(len(report['cases']), 2)
        self.assertTrue(output.with_suffix('.md').exists())
        args[args.index(str(output))] = str(self.path/'missing.json')
        p = subprocess.run(args+['--baseline',str(self.path/'absent.json')], capture_output=True)
        self.assertEqual(p.returncode, 1)


if __name__ == '__main__':
    unittest.main()
