"""R5/R8d native acceptance (pytest or unittest); skips without a built extension.

BAANDI_CLI: CLI from the same checkout/build profile as the extension.
BAANDI_SI_UPF: external Si.pz-vbc.UPF (never committed).
BAANDI_R5_EVIDENCE: optional NEW R5 evidence JSON destination.
BAANDI_R8D_EVIDENCE: optional NEW synthetic step/batch evidence JSON destination.
Run the Si2 tests explicitly: stress performs additional cutoff scans.
"""
import hashlib
import importlib.util
import json
import math
import re
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import warnings

try:
    from baandi import Options, Session
except ImportError:
    Options = Session = None

ROOT = Path(__file__).resolve().parents[3]
CELL = [[0., 5.13, 5.13], [5.13, 0., 5.13], [5.13, 5.13, 0.]]
# Displace the second atom so force equality also covers nonzero derivatives.
POSITIONS = [[0., 0., 0.], [.251, .249, .252]]


def synthetic():
    spec = importlib.util.spec_from_file_location("synthetic", ROOT / "scripts/gen_synthetic_upf.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.generate(2, (0,), radius=1.4, points=101)


def small_session(options):
    return Session([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
                   ["Xx"], [[0., 0., 0.]], {"Xx": synthetic()}, options)


def bits(values):
    if isinstance(values, (tuple, list)):
        return [bits(v) for v in values]
    return float(values).hex()


@unittest.skipIf(Session is None, "native extension not built/importable")
class ReferenceBindingTests(unittest.TestCase):
    def test_defaults_and_overrides(self):
        o = Options()
        self.assertEqual(o.driver, "reference")
        self.assertEqual(list(o.kgrid), [1, 1, 1])
        self.assertEqual((o.mixing_beta, o.mixing_history), (0.7, 8))
        self.assertTrue(o.verify)
        self.assertEqual(o.standard_threshold_ha, 5e-9)
        o.driver = "legacy"
        o.reference_strict = False
        self.assertEqual((o.mixing_beta, o.mixing_history), (0.5, 6))
        o.reference_strict = True
        self.assertEqual(o.mixing_beta, 0.3)
        o.mixing_beta, o.mixing_history = 0.4, 5
        o.driver = "reference"
        self.assertEqual((o.mixing_beta, o.mixing_history), (0.4, 5))

    def test_reference_invalid_start_fails_without_fallback(self):
        o = Options()
        o.reference_strict = False
        s = small_session(o)
        self.assertEqual(s.status, "Created")
        self.assertEqual(s.iteration, 0)
        for operation in (s.result, s.forces, s.density_summary, s.step):
            with self.assertRaises(RuntimeError):
                operation()
        with self.assertRaises(ValueError):
            s.stress()
        with self.assertRaisesRegex(RuntimeError, "requires compatible NC atomic charge"):
            s.initialize()
        self.assertEqual(s.status, "Failed")
        with self.assertRaises(RuntimeError):
            s.run()

    def test_reference_stepwise_batch_bit_identity(self):
        # Nonphysical atomic-start data, generated in memory; no UPF committed.
        text = synthetic()
        radii = [float(v) for v in re.search(r"<PP_R>(.*?)</PP_R>", text).group(1).split()]
        values = lambda xs: " ".join(format(v, ".17e") for v in xs)
        rho = values([r*r*math.exp(-r*r) for r in radii])
        text = re.sub(r"<PP_RHOATOM>.*?</PP_RHOATOM>", f"<PP_RHOATOM>{rho}</PP_RHOATOM>", text)
        text = text.replace('number_of_wfc="0"', 'number_of_wfc="2"')
        chi = "".join(f'<PP_CHI.{l+1} l="{l}" label="orbital" occupation="2">{values([r**(l+1)*math.exp(-0.6*r*r) for r in radii])}</PP_CHI.{l+1}>' for l in (0, 1))
        text = text.replace("</UPF>", f"<PP_PSWFC>{chi}</PP_PSWFC></UPF>")
        records = []
        for mesh in ([1, 1, 1], [2, 1, 1]):
            with self.subTest(kgrid=mesh):
                o = Options()
                o.reference_strict = False
                o.cutoff_hartree, o.bands = 2.5, 2
                o.kgrid = mesh
                def session():
                    return Session([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
                                   ["Xx"], [[0., 0., 0.]], {"Xx": text}, o)
                a, b = session(), session()
                self.assertEqual(a.status, "Created")
                a.initialize()
                self.assertEqual(a.status, "Initialized")
                self.assertEqual(a.iteration, 0)
                with self.assertRaises(RuntimeError):
                    a.initialize()
                reports = []
                while a.status == "Initialized":
                    report = a.step()
                    self.assertEqual(a.iteration, report.iteration)
                    self.assertIn(report.kind, ("Ordinary", "Verification", "Retry"))
                    self.assertEqual(report.stage, "fine")
                    self.assertEqual(report.estimator_ha, report.dr2_ry / 2.)
                    self.assertGreaterEqual(report.h_columns, 0)
                    self.assertGreater(report.ethr_ry, 0.)
                    self.assertIsInstance(report.safe_point, bool)
                    self.assertIsInstance(report.reused, bool)
                    reports.append(report)
                self.assertEqual(a.status, "Converged")
                self.assertIn("Verification", [r.kind for r in reports])
                with self.assertRaises(RuntimeError):
                    a.step()
                # Resume a partially stepped batch session, without replaying maps.
                b.initialize()
                b.step()
                b.run()
                ra, rb = a.result(), b.result()
                self.assertTrue(ra.converged and rb.converged)
                self.assertEqual(bits(ra.energy_hartree), bits(rb.energy_hartree))
                self.assertEqual(bits(ra.eigenvalues_hartree), bits(rb.eigenvalues_hartree))
                self.assertEqual(bits(a.density_summary()), bits(b.density_summary()))
                self.assertEqual(bits(a.density()), bits(b.density()))
                self.assertEqual(a.iteration, ra.iterations)
                self.assertEqual(a.iteration, b.iteration)
                c = session()
                c.run()
                self.assertEqual(bits(ra.energy_hartree), bits(c.result().energy_hartree))
                self.assertEqual(bits(ra.eigenvalues_hartree), bits(c.result().eigenvalues_hartree))
                self.assertEqual(bits(a.density_summary()), bits(c.density_summary()))
                self.assertEqual(bits(a.density()), bits(c.density()))
                self.assertEqual(bits(a.forces()), bits(c.forces()))
                a.run()  # retains delivery, does not run SCF again
                self.assertEqual(bits(a.result().energy_hartree), bits(ra.energy_hartree))
                records.append({"kgrid": mesh, "energy_ha": ra.energy_hartree,
                                "energy_hex": bits(ra.energy_hartree), "iterations": a.iteration,
                                "maps": len(reports), "max_energy_difference_ha": 0.,
                                "eigenvalues_bit_identical": True, "density_bit_identical": True, "forces_bit_identical": True,
                                "density_hex_sha256": hashlib.sha256(json.dumps(bits(a.density())).encode()).hexdigest()})
        if os.environ.get("BAANDI_R8D_EVIDENCE"):
            Path(os.environ["BAANDI_R8D_EVIDENCE"]).write_text(json.dumps({
                "task_id": "R8D-PY-MCP-STEP",
                "command": "PYTHONPATH=target/tmp/r8d-native python3 -m unittest discover -s crates/baandi-py/tests -p test_reference.py -v",
                "reference": "same-build reference batch run(), not QE",
                "scope": "synthetic cubic one-atom model; Gamma and k211; energy/eigenvalues/full density/forces",
                "cases": records}, indent=2) + "\n")

    def test_invalid_options(self):
        for name, value in (("driver", "typo"), ("mixing_beta", 0.),
                            ("mixing_history", 0), ("max_iterations", 0),
                            ("standard_threshold_ha", float("nan")),
                            ("reference_strict", True), ("kgrid", [0, 2, 2])):
            o = Options()
            o.reference_strict = False
            setattr(o, name, value)
            with self.subTest(option=name), self.assertRaises(ValueError):
                small_session(o)

    def test_legacy_rejects_non_gamma_mesh(self):
        o = Options()
        o.driver = "legacy"
        o.kgrid = [2, 2, 2]
        with self.assertRaisesRegex(ValueError, "non-Gamma kgrid requires"):
            small_session(o)

    def test_legacy_stepwise_batch_and_owned_result(self):
        o = Options()
        o.driver = "legacy"
        o.reference_strict = False
        o.cutoff_hartree, o.bands = 2., 2
        a, b = small_session(o), small_session(o)
        a.initialize()
        self.assertEqual(a.iteration, 0)
        while a.status == "Ready":
            a.step()
        b.run()
        ra, rb = a.result(), b.result()
        self.assertTrue(ra.converged and rb.converged)
        self.assertEqual(bits(ra.energy_hartree), bits(rb.energy_hartree))
        self.assertEqual(bits(ra.eigenvalues_hartree), bits(rb.eigenvalues_hartree))
        self.assertEqual(a.iteration, ra.iterations)
        self.assertAlmostEqual(a.density_summary()[-1], 2., places=10)
        del a, b
        self.assertTrue(ra.eigenvalues_hartree)

    @unittest.skipUnless(os.environ.get("BAANDI_CLI") and os.environ.get("BAANDI_SI_UPF"),
                         "set BAANDI_CLI and BAANDI_SI_UPF for real Si2 parity")
    def test_si2_gamma_and_k222_cli_bit_identity(self):
        pp_path = Path(os.environ["BAANDI_SI_UPF"]).resolve()
        pp = pp_path.read_text()
        records = []
        for mesh in ([1, 1, 1], [2, 2, 2]):
            with self.subTest(kgrid=mesh), tempfile.TemporaryDirectory(prefix="r5-python-") as folder:
                folder = Path(folder)
                o = Options()
                o.reference_strict = False
                o.cutoff_hartree, o.bands = 4., 8
                o.standard_threshold_ha = 5e-11
                o.kgrid = mesh
                text = ('[system]\nname="R5 Si2"\nlattice_unit="bohr"\n'
                        f'lattice_vectors={json.dumps(CELL)}\n')
                for pos in POSITIONS:
                    text += f'[[system.atoms]]\nspecies="Si"\nfractional={json.dumps(pos)}\n'
                text += (f'[pseudopotentials]\nSi={json.dumps(str(pp_path))}\n'
                         f'[basis]\necutwfc_ry=8.0\nkgrid={json.dumps(mesh)}\n'
                         '[xc]\nfunctional="lda-pz"\n[scf]\ndriver="reference"\n'
                         'nbands=8\nstandard_threshold_ha=5e-11\n'
                         '[output]\nforces=true\nstress=true\n')
                path = folder / "input.toml"
                path.write_text(text)
                command = [os.environ["BAANDI_CLI"], "run", str(path),
                           "--output-dir", str(folder / "out")]
                cli = subprocess.run(command, capture_output=True, text=True, timeout=600)
                self.assertEqual(cli.returncode, 0, cli.stdout + cli.stderr)
                expected = json.loads((folder / "out/results.json").read_text())
                self.assertEqual(expected["scf"]["driver"], "reference")
                self.assertTrue(expected["verified"])
                s = Session(CELL, ["Si", "Si"], POSITIONS, {"Si": pp}, o)
                s.run()
                r = s.result()
                self.assertTrue(r.converged)
                # Stress first must perform the facade's required force delivery.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    stress = s.stress()
                forces = s.forces()
                self.assertEqual(bits(r.energy_hartree), bits(expected["total_energy"]["ha"]))
                self.assertEqual(bits(forces), bits(expected["forces"]["values"]))
                self.assertEqual(bits(stress), bits(expected["stress"]["values"]))
                self.assertGreater(max(abs(v) for row in forces for v in row), 1e-6)
                self.assertEqual(r.iterations, expected["convergence"]["iterations"])
                self.assertAlmostEqual(s.density_summary()[-1], 8., places=9)
                s.run()  # repeated calls retain the same delivered state
                self.assertEqual(bits(s.result().energy_hartree), bits(r.energy_hartree))
                self.assertEqual(bits(s.forces()), bits(forces))
                records.append({"kgrid": mesh, "input_toml": text,
                                "cli_command": command, "driver": expected["scf"]["driver"],
                                "cli_tool": expected["tool"],
                                "cli_energy_hex": bits(expected["total_energy"]["ha"]),
                                "cli_forces_hex": bits(expected["forces"]["values"]),
                                "cli_stress_hex": bits(expected["stress"]["values"]),
                                "upf_sha256": hashlib.sha256(pp_path.read_bytes()).hexdigest(),
                                "energy_ha": r.energy_hartree, "energy_hex": bits(r.energy_hartree),
                                "forces_ha_bohr": forces, "forces_hex": bits(forces),
                                "stress_ha_bohr3": stress, "stress_hex": bits(stress),
                                "iterations": r.iterations, "max_energy_difference_ha": 0.,
                                "max_force_difference_ha_bohr": 0.,
                                "max_stress_difference_ha_bohr3": 0.})
        self.assertEqual(len(records), 2, "do not publish incomplete parity evidence")
        if os.environ.get("BAANDI_R5_EVIDENCE"):
            Path(os.environ["BAANDI_R5_EVIDENCE"]).write_text(
                json.dumps({"task": "R5-PYTHON-RUN", "reference": "same-build BAANDI CLI, not QE",
                            "cases": records}, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
