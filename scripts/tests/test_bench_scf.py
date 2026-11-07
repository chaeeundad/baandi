"""The optional historical-anchor bypass must not bypass baseline gates."""
import contextlib
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
    def run_gate(self, skip, delta=0.0, density="same"):
        result = dict(wall_median_s=1., engine_total_s=1.,
                      phases_s={"diagonalization_s": 0.5}, iterations=10,
                      total_energy_ha=-30., density_sha256="same")
        baseline = dict(result, total_energy_ha=-30. + delta,
                        density_sha256=density)
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
