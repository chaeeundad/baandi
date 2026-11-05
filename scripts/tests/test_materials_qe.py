"""Metadata goldens plus synthetic stdout/density tests; no new QE reference claims."""
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import shlex
import struct
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from make_qe_input import generate, MANIFEST
from qe_common import Error, EV_PER_HA, read_json, write_json
from compare_materials_qe import compare, parse_output, main
from qe_rhog import write_qe

FIX = ROOT/'scripts/tests/fixtures'


def output(r):
    """Invented stdout for exact perturbation tests, not a pw.x measurement."""
    points = r['kpoints']['points']
    lines = ['Program PWSCF v.TEST starts synthetic fixture',
             f"number of electrons = {r['electrons']}",
             f"number of Kohn-Sham states = {r['bands']}",
             f"number of k points = {len(points)}",
             f"kinetic-energy cutoff = {r['cutoff']['ecutwfc_ry']} Ry",
             f"charge density cutoff = {4*r['cutoff']['ecutwfc_ry']} Ry",
             'Exchange-correlation = SLA PW PBX PBC',
             'Dense grid: 1 G-vectors FFT dimensions: ( ' + ', '.join(map(str,r['basis']['fft_grid'])) + ' )']
    for title in ('cart. coord. in units 2pi/alat', 'cryst. coord.'):
        lines.append(title)
        for i,p in enumerate(points):
            lines.append(f"k( {i+1}) = ( " + ' '.join(map(str,p['fractional'])) + f"), wk = {2*p['weight']}")
        lines.append('')
    lines += ['convergence has been achieved in 7 iterations', 'End of self-consistent calculation']
    for p in points:
        lines += ['k = ' + ' '.join(map(str,p['fractional'])) + ' ( 1 PWs) bands (ev):', '',
                  ' '.join(format(x*EV_PER_HA,'.17g') for x in p['eigenvalues']['ha']), '',
                  'occupation numbers', ' '.join(map(str,r['occupations'])), '']
    lines += [f"! total energy = {2*r['total_energy']['ha']} Ry", 'JOB DONE.']
    return '\n'.join(lines)+'\n'


class MaterialsTests(unittest.TestCase):
    def setUp(self):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        self.tmp = Path(t.name)
        self.source = self.tmp/'results.json'
        self.out = self.tmp/'qe.in'
        self.r = read_json(FIX/'materials/si-20.json')

    def generate(self, **kw):
        write_json(self.source,self.r)
        return generate(self.source,'/never-opened',None,self.out,command='synthetic test invocation',**kw)

    def test_legacy_byte_goldens(self):
        folder = FIX/'qe_legacy'
        for row in read_json(folder/'index.json'):
            with self.subTest(source=row['source']):
                generate(ROOT/row['source'],'/legacy/pseudo','Si.pz-vbc.UPF',self.out)
                self.assertEqual(self.out.read_bytes(),(folder/row['golden']).read_bytes())

    def test_material_goldens(self):
        for name in ('si-20','mgo-40'):
            args = ['python3','-I','scripts/make_qe_input.py','--results',f'scripts/tests/fixtures/materials/{name}.json',
                    '--pseudo-dir','/scr/mino/rustdft/pseudo','--out',f'tests/reference/materials/qe_inputs/{name}.in']
            generate(ROOT/args[4],args[6],None,self.out,command=shlex.join(args))
            self.assertEqual(self.out.read_bytes(),(ROOT/args[8]).read_bytes())

    def test_isolated_cli_and_options(self):
        write_json(self.source,self.r)
        args = [sys.executable,'-I',str(ROOT/'scripts/make_qe_input.py'),'--results',str(self.source),
                '--pseudo-dir','/not/present','--out',str(self.out),'--disk-io','low','--outdir','./work with spaces',
                '--conv-thr','1e-10']
        p = subprocess.run(args,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        text = self.out.read_text()
        self.assertIn("disk_io = 'low'",text)
        self.assertIn("outdir = './work with spaces'",text)
        self.assertIn('! Exact command: '+shlex.join(['python3',*args[1:]]),text)

    def test_rejections(self):
        original = copy.deepcopy(self.r)
        changes = [dict(spin={'mode':'collinear'}),dict(smearing='fermi-dirac'),
                   dict(occupation_policy='fermi_dirac'),dict(schema_version=2),
                   dict(occupations=[1]*self.r['bands']),dict(input_sha256='missing')]
        for change in changes:
            self.r = dict(copy.deepcopy(original),**change)
            with self.subTest(change=change), self.assertRaises(Error):
                self.generate()
        self.r = original
        for kw in (dict(disk_io='none'),dict(fft_grid='auto'),dict(ecut_ry=80)):
            with self.assertRaises(Error):
                self.generate(**kw)
        m = read_json(MANIFEST)
        for d in m['datasets']:
            d.pop('sha256')
        path = self.tmp/'manifest.json'; write_json(path,m)
        with self.assertRaisesRegex(Error,'manifest hash'):
            self.generate(manifest=path)
        with self.assertRaisesRegex(Error,'summary'):
            generate(ROOT/'tests/reference/materials/results.json','/unused',None,self.out,command='test')

    def test_multi_species_hash_filename_and_overwrite(self):
        self.r = read_json(FIX/'materials/mgo-40.json')
        self.generate()
        text = self.out.read_text()
        self.assertIn('Mg 1.0 Mg.upf\nO 1.0 O.upf',text)
        self.r['species'][1]['sha256']='a'*64
        with self.assertRaisesRegex(Error,'manifest hash'):
            self.generate()
        self.r = read_json(FIX/'materials/mgo-40.json')
        self.r['species'][1]['path']='/unused/wrong.upf'
        with self.assertRaisesRegex(Error,'filename'):
            self.generate()
        self.r = read_json(FIX/'materials/mgo-40.json'); write_json(self.source,self.r)
        before = self.source.read_bytes()
        with self.assertRaises(Error):
            generate(self.source,'/unused',None,self.source,command='test')
        self.assertEqual(self.source.read_bytes(),before)

    def prepare_comparison(self):
        self.r['kpoints'] = dict(count=2,grid=None,shift=None,time_reversal=False,points=[
            dict(fractional=[0,0,0],weight=.25,eigenvalues={'ha':list(range(self.r['bands']))}),
            dict(fractional=[.25,0,0],weight=.75,eigenvalues={'ha':[x+.1 for x in range(self.r['bands'])]})])
        write_json(self.source,self.r)
        self.qe = self.tmp/'qe.out'; self.qe.write_text(output(self.r))

    def test_comparison_values_reordering_and_no_caps(self):
        self.prepare_comparison()
        q = copy.deepcopy(self.r)
        q['kpoints']['points'].reverse()
        q['total_energy']['ha'] += .125
        q['kpoints']['points'][0]['eigenvalues']['ha'][4] += .0625
        self.qe.write_text(output(q))
        report = compare(self.source,self.qe)
        self.assertEqual(report['total_energy_ha']['difference'],-.125)
        self.assertEqual(report['total_energy_ry']['difference'],-.25)
        self.assertEqual(len(report['eigenvalues']),10)
        self.assertEqual(report['eigenvalues'][0]['qe_k_index'],2)
        self.assertAlmostEqual(report['eigenvalues'][-1]['difference'],-.0625)
        self.assertEqual(report['iterations']['qe'],7)
        self.assertIsNone(report['tolerances'])
        self.assertNotIn('passed',report)
        self.assertIsNone(report['density'])

    def test_gap_differences_and_missing_empty_band(self):
        self.prepare_comparison()
        q = copy.deepcopy(self.r)
        for p in q['kpoints']['points']:
            p['eigenvalues']['ha'][4] += .25
        self.qe.write_text(output(q))
        report = compare(self.source,self.qe)
        for gap in report['band_gap'].values():
            self.assertAlmostEqual(gap['difference'],-.25)
        self.r['bands'] = 4
        self.r['occupations'] = [2]*4
        for p in self.r['kpoints']['points']:
            p['eigenvalues']['ha'] = p['eigenvalues']['ha'][:4]
        write_json(self.source,self.r)
        self.qe.write_text(output(self.r))
        self.assertIsNone(compare(self.source,self.qe)['band_gap'])

    def test_bad_stdout_and_conditions(self):
        self.prepare_comparison()
        good = self.qe.read_text()
        for bad in (good+good,good.replace('JOB DONE.',''),good.replace('bands (ev):','bands broken:'),
                    good.replace('wk = 0.5','wk = 0.7'),good.replace('SLA PW PBX PBC','SLA PZ'),
                    good.replace('30, 30, 30','31, 30, 30'),good.replace('k = 0.25 0 0','k = 0.5 0 0')):
            self.qe.write_text(bad)
            with self.assertRaises(Error):
                compare(self.source,self.qe)
        self.qe.write_text(good)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(['--ours',str(self.source),'--qe-out',str(self.qe),'--json',str(self.source)]),1)

    def test_density_reuses_independent_transform(self):
        self.prepare_comparison()
        self.r.update(electrons=2,bands=3,occupations=[2,0,0])
        self.r['cutoff']=dict(ecutwfc_ha=.1,ecutwfc_ry=.2)
        lattice=[[2,0,0],[0,2,0],[0,0,2]]
        self.r['cell']['lattice_columns_bohr']=lattice
        self.r['basis'].update(fft_grid=[3,3,3],n_density_g=1)
        for p in self.r['kpoints']['points']:
            p['eigenvalues']['ha']=[0,1,2]
        write_json(self.source,self.r); self.qe.write_text(output(self.r))
        raw=struct.pack('<27d',*([.25]*27)); (self.tmp/'density.bin').write_bytes(raw)
        write_json(self.tmp/'density.json',dict(schema_version=1,dtype='f64',endianness='little',
            ordering='C: (i*n2+j)*n3+k; k fastest',units='electrons per Bohr^3',fft_grid=[3,3,3],
            sha256=hashlib.sha256(raw).hexdigest(),volume_bohr3=8,lattice_columns_bohr=lattice,
            normalization=dict(expected_electrons=2,integral_electrons=2)))
        density=self.tmp/'charge-density.dat'; write_qe(density,{(0,0,0):.25+0j})
        report=compare(self.source,self.qe,density)['density']
        self.assertEqual(report['relative_l2_difference'],0)
        self.assertEqual(report['qe_integral_electrons'],2)
        self.assertNotIn('passed',report)
        density.write_bytes(b'broken')
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main(['--ours',str(self.source),'--qe-out',str(self.qe),'--qe-density',str(density)]),3)

    def test_existing_real_stdout(self):
        # Existing independently archived reference only; no new material QE run.
        for name in ('k222_c','k222s_c'):
            folder=ROOT/'tests/reference/kpoints'/name
            report=compare(folder/'ours_results.json',folder/'qe.out')
            self.assertEqual(report['status'],'REPORTED')
            self.assertLess(report['total_energy_ha']['absolute_difference'],1e-4)


if __name__ == '__main__':
    unittest.main()
