#!/usr/bin/env python3
"""Runner contract tests with fake launchers; no MPI or scientific calculation."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

LANE=Path(__file__).resolve().parent
ROOT=LANE.parents[3]


class RunnerContracts(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(dir=ROOT/'target/scratch')
        self.addCleanup(self.tmp.cleanup)
        self.work=Path(self.tmp.name)
        self.env=dict(os.environ, PATH=str(self.work)+os.pathsep+os.environ['PATH'])
        for name,shifts in [('taskset',2),('mpirun',6)]:
            p=self.work/name
            p.write_text('#!/bin/sh\nshift '+str(shifts)+'\nexec "$@"\n')
            p.chmod(0o755)
        cli=self.work/'cli'
        cli.write_text('#!'+sys.executable+'\n'+'''import json,sys,pathlib
p=pathlib.Path(sys.argv[-1]); p.mkdir(parents=True)
d={'convergence':{'converged':True,'history':[{'standard_error_ha':1e-12,'density_residual':1e-8,'pulay':None,'mixer_reset':False}]},'total_energy':{'ha':-1.0},'verified':True}
(p/'results.json').write_text(json.dumps(d)); (p/'density.bin').write_bytes(b'fake')
''')
        cli.chmod(0o755)
        (self.work/'input').write_text('fake')
        (self.work/'pseudo').write_text('fake, not a UPF')
        (self.work/'strict.json').write_text(json.dumps({'single_path':{'effective_profile':'strict'},
            'convergence':{'converged':True,'history':[{'fixed_point_certificate':{}}]},
            'total_energy':{'ha':-1.0}}))

    def command(self,budget=600,name='test-aa'):
        return [sys.executable,'-I',str(LANE/'dev08_pairs.py'),'--binary',str(self.work/'cli'),
            '--input',str(self.work/'input'),'--upf',str(self.work/'pseudo'),
            '--strict-result',str(self.work/'strict.json'),'--out',str(self.work/'out'),
            '--job-name',name,'--sequence','AA','--budget',str(budget)]

    def test_aa_records_identity_and_acceptance(self):
        p=subprocess.run(self.command(),env=self.env,capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr)
        d=json.loads((self.work/'out/job.json').read_text())
        self.assertTrue(d['complete_and_correct'])
        self.assertEqual(len(d['runs']),2)
        self.assertEqual(d['runs'][0]['density_sha256'],d['runs'][1]['density_sha256'])
        self.assertTrue(all(r['environment']['BAANDI_MIX_CELL_KERKER']=='off' for r in d['runs']))
        self.assertEqual(d['runs'][0]['predicted_ratio_trace'],[None])

    def test_budget_caps_reject_before_launch(self):
        for name,budget in [('test-aa',601),('final-test',7201)]:
            p=subprocess.run(self.command(budget,name),env=self.env,capture_output=True,text=True)
            self.assertEqual(p.returncode,2)
            self.assertFalse((self.work/'out').exists())

    def test_uncertified_reference_rejects(self):
        (self.work/'strict.json').write_text(json.dumps({'convergence':{'converged':False}}))
        p=subprocess.run(self.command(),env=self.env,capture_output=True,text=True)
        self.assertEqual(p.returncode,2)
        self.assertFalse((self.work/'out/job.json').exists())


if __name__=='__main__': unittest.main()
