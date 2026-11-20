#!/usr/bin/env python3
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location("lane_init_parse", ROOT / "scripts/lane_init_parse.py")
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def parse(self, counters, rank=0):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "results.json"
            path.write_text(json.dumps({"run": {"rank": rank, "ranks": 8, "counters": counters}}))
            return parser.parse(path)

    def test_disabled_and_count_only_are_not_zero_time_measurements(self):
        counters = {"init_lane.problem.calls": 0, "init_lane.problem.ns": 0}
        self.assertEqual(self.parse(counters)["state"], "disabled_or_not_reached")
        counters["init_lane.problem.calls"] = 2
        row = self.parse(counters)
        self.assertEqual(row["state"], "counts_only")
        self.assertIsNone(row["spans"][0]["seconds"])

    def test_rank_local_timer_units_and_payload(self):
        row = self.parse({"init_lane.problem.calls": 2, "init_lane.problem.ns": 6_000_000,
                          "init_lane.problem.grid_points": 100}, rank=5)
        self.assertEqual(row["rank"], 5)
        self.assertEqual(row["spans"][0]["seconds_per_call"], 0.003)
        self.assertEqual(row["payload_floors"]["one_complex_grid_bytes_per_problem_sum"], 1600)

    def test_missing_malformed_and_incomplete_records_reject(self):
        for counters in ({}, {"init_lane.problem.calls": 1},
                         {"init_lane.problem.calls": -1, "init_lane.problem.ns": 0},
                         {"init_lane.problem.calls": 0, "init_lane.problem.ns": 1},
                         {"init_lane.problem.calls": True, "init_lane.problem.ns": 0}):
            with self.assertRaises(ValueError):
                self.parse(counters)


if __name__ == "__main__":
    unittest.main()
