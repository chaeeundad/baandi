"""Independent ASE optimizer and format readers against the CLI relaxation.

Run with the current native wheel and BAANDI_CLI pointing at the release CLI.
The synthetic periodic H-labelled model is not a physical hydrogen prediction.
"""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_binding import potential, Options, Session


@unittest.skipIf(Session is None, "native extension unavailable")
@unittest.skipUnless(os.environ.get("BAANDI_CLI"), "set BAANDI_CLI")
class RelaxationTests(unittest.TestCase):
    def test_ase_bfgs_matches_cli_and_reads_formats(self):
        from ase import Atoms
        from ase.constraints import FixAtoms
        from ase.optimize import BFGS
        from ase.io import read
        from baandi.ase import Baandi
        from baandi._baandi import angstrom_to_bohr, hartree_to_ev
        import numpy as np
        ang = 1.0 / angstrom_to_bohr(1.0)
        ev = hartree_to_ev(1.0)
        text = potential("H", 1, (0,))
        opts = Options()
        opts.cutoff_hartree = 1.5
        opts.bands = 4
        opts.max_iterations = 300
        opts.density_tolerance = 1e-10
        opts.energy_tolerance = 1e-12
        opts.eigen_tolerance = 1e-11
        atoms = Atoms("H2", cell=np.eye(3)*8*ang,
                      scaled_positions=[[0, 0, 0], [.4, .4, .4]], pbc=True)
        atoms.set_constraint(FixAtoms(indices=[0]))
        atoms.calc = Baandi(pseudopotentials={"H": text}, options=opts)
        optimizer = BFGS(atoms, logfile=None)
        self.assertTrue(optimizer.run(fmax=1e-8*ev/ang, steps=200))
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "H.upf").write_text(text)
            (root / "input.toml").write_text('''[system]
name='synthetic periodic relaxation'
lattice_unit='bohr'
lattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]
[[system.atoms]]
species='H'
fractional=[0.0,0.0,0.0]
[[system.atoms]]
species='H'
fractional=[0.4,0.4,0.4]
[pseudopotentials]
H='H.upf'
[basis]
ecutwfc_ry=3.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional='lda-pz'
[scf]
nbands=4
max_iterations=300
density_residual_tolerance=1e-10
energy_tolerance_ha=1e-12
eigen_tolerance=1e-11
[relax]
max_steps=200
force_tol_ha_bohr=1e-8
energy_tol_ha=1e-12
fixed=[{atom=0}]
''')
            run = subprocess.run([os.environ["BAANDI_CLI"], "relax", str(root / "input.toml"),
                                  "--output-dir", str(root / "out")], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            result = json.loads((root / "out/results.json").read_text())["relaxation"]
            final = result["frames"][-1]
            delta = np.max(np.abs(atoms.positions/ang - np.array(final["positions"]).reshape(-1, 3)))
            self.assertLess(delta, 1e-4)
            de = abs(atoms.get_potential_energy()/ev - final["energy_ha"])
            self.assertLess(de, 1e-10)
            trajectory = read(root / "out/trajectory.extxyz", index=":")
            self.assertEqual(len(trajectory), len(result["frames"]))
            np.testing.assert_allclose(trajectory[-1].positions/ang, np.array(final["positions"]).reshape(-1, 3), atol=1e-13, rtol=0)
            self.assertAlmostEqual(trajectory[-1].get_potential_energy()/ev, final["energy_ha"], places=13)
            np.testing.assert_allclose(trajectory[-1].get_forces()*ang/ev, final["forces_ha_bohr"], atol=1e-14, rtol=0)
            poscar = read(root / "out/POSCAR", format="vasp")
            np.testing.assert_allclose(poscar.positions, trajectory[-1].positions, atol=1e-14, rtol=0)
            print(f"ASE BFGS vs CLI: max coordinate difference={delta:.17e} Bohr; energy difference={de:.17e} Ha; ASE steps={optimizer.nsteps}; CLI steps={result['steps']}; extxyz/POSCAR independent readers passed")
