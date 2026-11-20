#!/usr/bin/env python3
"""Parser checks using a completed CLI run; no performance measurements."""
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("lane_parser", "scripts/lane_launchio_parse.py")
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)
SOURCE = Path(sys.argv.pop(1)).resolve()


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir="target/scratch")
        self.root = Path(self.temp.name)
        for name in ("results.json", "density.json", "density.bin"):
            shutil.copyfile(SOURCE / name, self.root / name)
        self.doc = json.loads((self.root / "results.json").read_text())

    def tearDown(self):
        self.temp.cleanup()

    def save(self):
        (self.root / "results.json").write_text(json.dumps(self.doc))

    def test_real_output_and_optional_external_resources(self):
        self.assertIsNone(parser.parse_run(self.root)["external"])
        (self.root / "launch-wall.json").write_text('{"wall_s": 100, "exit_code": 0}')
        (self.root / "resource-rank-0.json").write_text('{"exit_code": 0, "minor_faults": 12, "max_rss_kib": 128}')
        out = parser.parse_run(self.root)
        self.assertEqual(out["external"]["wall_s"], 100)
        self.assertEqual(out["per_rank"][0]["resources"]["minor_faults"], 12)

    def test_rejects_missing_timer_incomplete_and_nonfinite(self):
        original = copy.deepcopy(self.doc)
        mutations = [lambda d: d["launch"]["seconds"].pop("scf_driver"),
                     lambda d: d["run"].update(completed=False),
                     lambda d: d["launch"].update(main_wall_s=float("nan")),
                     lambda d: d["launchio"]["sizes"].update(density_bytes=1),
                     lambda d: d["launch"].update(covered_s=1000)]
        for mutate in mutations:
            self.doc = copy.deepcopy(original)
            mutate(self.doc)
            self.save()
            with self.assertRaises(ValueError):
                parser.parse_run(self.root)

    def test_density_corruption_rejects(self):
        with (self.root / "density.bin").open("ab") as f:
            f.write(b"x")
        with self.assertRaisesRegex(ValueError, "density hash"):
            parser.parse_run(self.root)

    def test_duplicate_keys_and_truncation_reject(self):
        (self.root / "results.json").write_text('{"run": 1, "run": 2}')
        with self.assertRaisesRegex(ValueError, "duplicate JSON"):
            parser.parse_run(self.root)
        (self.root / "results.json").write_text('{"run":')
        with self.assertRaises(ValueError):
            parser.parse_run(self.root)

    def test_synthetic_two_rank_transport_and_mixed_run_rejection(self):
        # Transport-only fixture, never MPI execution evidence.
        self.doc["launch"].update(ranks=2, mpi_backend="mpi")
        self.doc["parallel"] = {"model": "mpi_stick_rows"}
        self.save()
        with self.assertRaisesRegex(ValueError, "missing/duplicate ranks"):
            parser.parse_run(self.root)
        rank_dir = self.root / "rank1"
        rank_dir.mkdir()
        for name in ("density.json", "density.bin"):
            shutil.copyfile(self.root / name, rank_dir / name)
        other = copy.deepcopy(self.doc)
        other["launch"]["rank"] = other["launchio"]["rank"] = 1
        (rank_dir / "results.json").write_text(json.dumps(other))
        self.assertEqual(parser.parse_run(self.root)["ranks"], 2)
        other["run"]["run_id"] = "different-run"
        (rank_dir / "results.json").write_text(json.dumps(other))
        with self.assertRaisesRegex(ValueError, "mixed run IDs"):
            parser.parse_run(self.root)


if __name__ == "__main__":
    unittest.main()
