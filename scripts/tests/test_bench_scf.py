"""The optional historical-anchor bypass must not bypass baseline gates."""
import contextlib
import copy
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "bench_scf", Path(__file__).resolve().parents[1] / "bench_scf.py")
BENCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCH)


class AnchorSelection(unittest.TestCase):
    def run_gate(self, skip, delta=0.0, density="same", historical=False,
                 baseline_reference="off", run_driver="strict_oracle"):
        result = dict(wall_median_s=1., engine_total_s=1.,
                      phases_s={"diagonalization_s": 0.5}, iterations=10,
                      total_energy_ha=-30., density_sha256="same",
                      scf={"driver": run_driver, "driver_source": "env",
                           "driver_reason": "strict request"}, effective_profile="strict")
        baseline = dict(result, total_energy_ha=-30. + delta,
                        density_sha256=density)
        if historical:
            baseline["scf"] = {
                "driver": "legacy",
                "reference_switches": {
                    "BAANDI_REFERENCE_SCF": {"value": baseline_reference},
                    "BAANDI_REFERENCE_STRICT": {"value": "on"}}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            binary, upf, stored = [root / n for n in ("binary", "upf", "base.json")]
            binary.touch()
            upf.touch()
            stored.write_text(json.dumps({"cases": {"si8": baseline}}))
            args = ["--binary", str(binary), "--upf", str(upf), "--cases", "si8",
                    "--baseline", str(stored), "--require-identical"]
            if skip:
                args.append("--skip-legacy-si8-anchor")
            with patch.object(BENCH, "run_case", return_value=result), \
                    contextlib.redirect_stdout(io.StringIO()):
                return BENCH.main(args)

    def test_legacy_anchor_remains_default(self):
        self.assertEqual(self.run_gate(False), 5)

    def test_explicit_bypass_accepts_identical_profile_baseline(self):
        self.assertEqual(self.run_gate(True), 0)

    def test_explicit_bypass_keeps_energy_and_density_gates(self):
        self.assertEqual(self.run_gate(True, delta=1e-8), 5)
        self.assertEqual(self.run_gate(True, density="different"), 5)

    def test_historical_oracle_label_keeps_all_numerical_gates(self):
        self.assertEqual(self.run_gate(True, historical=True), 0)
        self.assertEqual(self.run_gate(True, historical=True, delta=1e-8), 5)
        self.assertEqual(self.run_gate(True, historical=True, density="different"), 5)

    def test_ambiguous_historical_switches_do_not_authorize_oracle_identity(self):
        self.assertEqual(self.run_gate(True, historical=True, baseline_reference="on"), 3)

    def test_reference_run_cannot_pass_as_historical_oracle(self):
        self.assertEqual(self.run_gate(True, historical=True, run_driver="reference"), 3)


class RoutingMetadata(unittest.TestCase):
    def test_default_baseline_checks_reference_without_selecting_switches(self):
        root = Path(__file__).resolve().parents[2]
        stored = json.loads((root / "docs/performance/baseline-reference.json").read_text())
        result = copy.deepcopy(stored["cases"]["gamma"])
        result["scf"]["driver_source"] = "default"
        result["scf"]["driver_reason"] = None
        with tempfile.TemporaryDirectory() as tmp:
            binary, upf = [Path(tmp) / n for n in ("binary", "upf")]
            binary.touch()
            upf.touch()
            args = ["--binary", str(binary), "--upf", str(upf), "--cases", "gamma",
                    "--require-identical"]
            with patch.object(BENCH, "run_case", return_value=result) as run, \
                    patch.dict(BENCH.os.environ, {}, clear=True), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(BENCH.main(args), 0)
                run.assert_called_once_with(str(binary), str(upf), "gamma", 3, None, [])
                self.assertFalse(any(k.startswith("BAANDI_") for k in BENCH.os.environ))
                result["density_sha256"] = "changed"
                self.assertEqual(BENCH.main(args), 5)
                self.assertEqual(BENCH.main(args + ["--baseline", str(Path(tmp) / "missing")]), 1)
                self.assertEqual(BENCH.main(args + ["--cases", "k444s"]), 1)

    def test_template_does_not_force_legacy_but_explicit_controls_remain(self):
        for case in BENCH.CASES:
            text = BENCH.case_input("synthetic.upf", case)
            self.assertNotIn("energy_tolerance_ha", text)
            self.assertNotIn("density_residual_tolerance", text)
        text = BENCH.case_input("synthetic.upf", "gamma", extra_toml=["driver = 'reference'"])
        self.assertIn("driver = 'reference'", text)
        text = BENCH.case_input("synthetic.upf", "gamma", extra_toml=["energy_tolerance_ha = 1e-8"])
        self.assertIn("energy_tolerance_ha = 1e-8", text)

    def test_selected_driver_provenance_is_copied_from_result(self):
        scf = {"driver": "legacy", "driver_reason": "L11: missing atomic data",
               "driver_source": "input", "switches": {
                   "BAANDI_REFERENCE_SCF": {"value": "on"},
                   "BAANDI_REFERENCE_STRICT": {"value": "off"}}}
        metadata = BENCH.driver_metadata({"scf": scf})
        for key in ("driver", "driver_reason", "driver_source"):
            self.assertEqual(metadata[key], scf[key])
        for key in ("driver", "driver_source", "driver_reason", "switches"):
            incomplete = dict(scf)
            del incomplete[key]
            with self.assertRaises(RuntimeError):
                BENCH.driver_metadata({"scf": incomplete})

    def test_baseline_normalization_is_narrow_and_does_not_mutate_evidence(self):
        root = Path(__file__).resolve().parents[2]
        for name, expected in [("reference", "reference"), ("oracle", "strict_oracle")]:
            stored = json.loads((root / f"docs/performance/baseline-{name}.json").read_text())
            before = copy.deepcopy(stored)
            for case in stored["cases"].values():
                self.assertEqual(BENCH.baseline_driver(case), expected)
                if name == "oracle":
                    for field in ("effective_profile", "scf"):
                        incomplete = copy.deepcopy(case)
                        del incomplete[field]
                        self.assertNotEqual(BENCH.baseline_driver(incomplete), "strict_oracle")
                    for switch in BENCH.DRIVER_SWITCHES:
                        incomplete = copy.deepcopy(case)
                        del incomplete["scf"]["reference_switches"][switch]
                        self.assertEqual(BENCH.baseline_driver(incomplete), "legacy")
            self.assertEqual(stored, before)
