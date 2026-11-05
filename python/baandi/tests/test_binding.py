"""Run with unittest discovery; extension and CLI-dependent tests skip explicitly.

Synthetic UPFs are generated in a temporary directory, never committed.
BAANDI_CLI points to a CLI built with the same profile/toolchain as the wheel.
BAANDI_EVIDENCE optionally records measured comparisons as JSON.
"""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

try:
    from baandi import Options, Session
except ImportError as exc:
    Options = Session = None
    EXTENSION_ERROR = str(exc)

ROOT = Path(__file__).resolve().parents[3]


def potential(symbol, valence, channels):
    spec = importlib.util.spec_from_file_location("synthetic", ROOT / "scripts/gen_synthetic_upf.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.generate(valence, channels, radius=1.4, points=401).replace('element="Xx"', f'element="{symbol}"')


CASES = [
    ("cubic_h2", "H", 1, (0,), [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]], [[0., 0., 0.], [.5, .5, .5]]),
    ("skew_he", "He", 2, (0, 1), [[4., 0., 0.], [.5, 4., 0.], [0., .5, 5.]], [[0., 0., 0.]]),
]


def settings():
    o = Options()
    o.cutoff_hartree = 4.
    o.bands = 4
    return o


@unittest.skipIf(Session is None, "native extension unavailable: " + globals().get("EXTENSION_ERROR", ""))
class BindingTests(unittest.TestCase):
    def test_two_systems_stepwise_and_owned_results(self):
        for name, symbol, z, channels, cell, positions in CASES:
            with self.subTest(system=name):
                pp = {symbol: potential(symbol, z, channels)}
                args = (cell, [symbol] * len(positions), positions, pp, settings())
                a = Session(*args)
                with self.assertRaises(RuntimeError):
                    a.result()
                with self.assertRaises(RuntimeError):
                    a.density_summary()
                a.initialize()
                self.assertEqual(a.iteration, 0)
                self.assertEqual(a.status, "Ready")
                a.step()
                self.assertEqual(a.iteration, 1)
                while a.status == "Ready":
                    a.step()
                ra = a.result()
                self.assertTrue(ra.converged)
                b = Session(*args)
                b.run()
                rb = b.result()
                self.assertEqual(ra.energy_hartree.hex(), rb.energy_hartree.hex())
                self.assertEqual(ra.eigenvalues_hartree, rb.eigenvalues_hartree)
                shape, low, high, mean, charge = a.density_summary()
                self.assertEqual(len(shape), 3)
                self.assertLessEqual(low, mean)
                self.assertLessEqual(mean, high)
                self.assertAlmostEqual(charge, z * len(positions), places=10)
                del a, b
                self.assertTrue(ra.eigenvalues_hartree)

    def test_invalid_inputs_and_exhaustion(self):
        _, symbol, z, channels, cell, positions = CASES[0]
        pp = {symbol: potential(symbol, z, channels)}
        with self.assertRaises(ValueError):
            Session(cell, [symbol], positions, pp)
        o = settings()
        o.max_iterations = 0
        with self.assertRaises(ValueError):
            Session(cell, [symbol] * 2, positions, pp, o)
        o.max_iterations = 1
        s = Session(cell, [symbol] * 2, positions, pp, o)
        s.run()
        self.assertFalse(s.result().converged)
        self.assertEqual(s.status, "NotConverged")
        with self.assertRaises(RuntimeError):
            s.forces()

    def test_ase_forces_and_units(self):
        try:
            from ase import Atoms
            from baandi.ase import Baandi
            from baandi._baandi import hartree_to_ev, angstrom_to_bohr
            import numpy as np
        except ImportError as exc:
            self.skipTest(str(exc))
        atoms = Atoms("H2", cell=[4., 4., 4.], scaled_positions=[[0., 0., 0.], [.27, .23, .19]], pbc=True)
        atoms.calc = Baandi(pseudopotentials={"H": potential("H", 1, (0,))}, options=settings())
        forces = atoms.get_forces()
        native = np.asarray(atoms.calc.session.forces())
        np.testing.assert_allclose(forces, native * hartree_to_ev(1.) * angstrom_to_bohr(1.), rtol=1e-14, atol=1e-14)
        self.assertEqual(forces.shape, (2, 3))
        self.assertTrue(np.isfinite(forces).all())
        self.assertGreater(np.max(np.abs(forces)), 1e-5)

    def test_ase_unsupported_properties(self):
        try:
            from baandi.ase import Baandi
        except ImportError as exc:
            self.skipTest(str(exc))
        with self.assertRaises(NotImplementedError):
            Baandi(pseudopotentials={}, hubbard={"H": 1.0})
        calc = Baandi(pseudopotentials={})
        for prop in ("stress", "free_energy", "dipole", "charges"):
            with self.subTest(property=prop), self.assertRaises(NotImplementedError):
                calc.get_property(prop)
        with self.assertRaises(NotImplementedError):
            calc.calculate(properties=["stress"])

    @unittest.skipUnless(os.environ.get("BAANDI_CLI"), "set BAANDI_CLI for exact CLI comparison")
    def test_native_cli_bitwise(self):
        evidence = []
        for name, symbol, z, channels, cell, positions in CASES:
            with self.subTest(system=name), tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                pp = potential(symbol, z, channels)
                (folder / "synthetic.upf").write_text(pp)
                text = f'[system]\nname="{name}"\nlattice_unit="bohr"\nlattice_vectors={json.dumps(cell)}\n'
                for pos in positions:
                    text += f'[[system.atoms]]\nspecies="{symbol}"\nfractional={json.dumps(pos)}\n'
                text += '[pseudopotentials]\n' + symbol + '= "synthetic.upf"\n[basis]\necutwfc_ry=8.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional="lda-pz"\n[scf]\nnbands=4\n'
                (folder / "input.toml").write_text(text)
                completed = subprocess.run([os.environ["BAANDI_CLI"], "run", str(folder / "input.toml"), "--output-dir", str(folder / "out")], capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                result = json.loads((folder / "out/results.json").read_text())
                session = Session(cell, [symbol] * len(positions), positions, {symbol: pp}, settings())
                session.run()
                native = session.result()
                self.assertTrue(native.converged)
                self.assertEqual(native.energy_hartree.hex(), result["total_energy"]["ha"].hex())
                self.assertEqual(native.eigenvalues_hartree, result["eigenvalues"]["ha"])
                evidence.append({"system": name, "energy_ha": native.energy_hartree, "energy_hex_ha": native.energy_hartree.hex(), "cli_energy_hex_ha": result["total_energy"]["ha"].hex(), "iterations": native.iterations, "eigenvalues_ha": native.eigenvalues_hartree, "charge_electrons": native.charge_electrons, "reference": "same-build BAANDI CLI; synthetic NC-format fixture, not QE"})
        if os.environ.get("BAANDI_NATIVE_EVIDENCE"):
            Path(os.environ["BAANDI_NATIVE_EVIDENCE"]).write_text(json.dumps(evidence, indent=2) + "\n")

    @unittest.skipUnless(os.environ.get("BAANDI_CLI"), "set BAANDI_CLI for exact CLI comparison")
    def test_cli_and_ase_bitwise(self):
        try:
            from ase import Atoms
            from baandi.ase import Baandi, ConvergenceError
        except ImportError as exc:
            self.skipTest(str(exc))
        from baandi._baandi import angstrom_to_bohr
        evidence = []
        for name, symbol, z, channels, _, positions in CASES:
            with self.subTest(system=name), tempfile.TemporaryDirectory() as folder:
                folder = Path(folder)
                pp = potential(symbol, z, channels)
                (folder / "synthetic.upf").write_text(pp)
                # Exactly representable ASE cell and scaled positions. Serialize
                # ASE's actual fractional coordinates to avoid inverse-rounding.
                cell = [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]] if name == "cubic_h2" else [[2., 0., 0.], [.25, 2., 0.], [0., .25, 2.5]]
                atoms = Atoms([symbol] * len(positions), scaled_positions=positions, cell=cell, pbc=True)
                fractions = atoms.get_scaled_positions(wrap=False).tolist()
                text = f'[system]\nname="{name}"\nlattice_unit="angstrom"\nlattice_vectors={json.dumps(cell)}\n'
                for pos in fractions:
                    text += f'[[system.atoms]]\nspecies="{symbol}"\nfractional={json.dumps(pos)}\n'
                text += '[pseudopotentials]\n' + symbol + '="synthetic.upf"\n[basis]\necutwfc_ry=8.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional="lda-pz"\n[scf]\nnbands=4\n'
                (folder / "input.toml").write_text(text)
                completed = subprocess.run([os.environ["BAANDI_CLI"], "run", str(folder / "input.toml"), "--output-dir", str(folder / "out")], capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                result = json.loads((folder / "out/results.json").read_text())
                calc = Baandi(pseudopotentials={symbol: pp}, options=settings())
                atoms.calc = calc
                energy = atoms.get_potential_energy()
                native = calc.session.result()
                self.assertEqual(energy.hex(), result["total_energy"]["ev"].hex())
                self.assertEqual(native.energy_hartree.hex(), result["total_energy"]["ha"].hex())
                self.assertEqual(native.eigenvalues_hartree, result["eigenvalues"]["ha"])
                evidence.append({"system": name, "energy_ha": native.energy_hartree, "energy_ev": energy, "energy_hex_ev": energy.hex(), "cli_energy_hex_ev": result["total_energy"]["ev"].hex(), "iterations": native.iterations, "reference": "same-build BAANDI CLI; synthetic NC-format fixture, not QE"})
                # ASE and native geometry conversion must agree for skew cells.
                direct = Session([[angstrom_to_bohr(x) for x in row] for row in cell], [symbol] * len(positions), fractions, {symbol: pp}, settings())
                direct.run()
                self.assertEqual(direct.result().energy_hartree.hex(), native.energy_hartree.hex())
                o = settings()
                o.max_iterations = 1
                atoms.calc = Baandi(pseudopotentials={symbol: pp}, options=o)
                with self.assertRaises(ConvergenceError):
                    atoms.get_potential_energy()
                self.assertEqual(atoms.calc.results, {})
        if os.environ.get("BAANDI_EVIDENCE"):
            Path(os.environ["BAANDI_EVIDENCE"]).write_text(json.dumps(evidence, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
