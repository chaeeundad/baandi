"""Phase-0 parser contract checks; run with python3 -I this file."""
import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location("lane_fft_parse", ROOT / "scripts/lane_fft_parse.py")
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


class ParserTests(unittest.TestCase):
    def setUp(self):
        self.fixture = {"phase_profile": {"ranks": 1, "fft_lane": {
            "schema": "lane-fft-1", "enabled": True, "stage_timers": "on",
            "per_rank": [{"rank": 0, "buckets": [{"bucket": "fine", "maps": 2,
                "stages": [{"stage": "multiply", "calls": 4, "seconds": 0.008,
                            "points": 100, "logical_bytes": 4000, "fft_flops_estimate": 0}]}]}]}}}

    def parse(self, data):
        with tempfile.TemporaryDirectory() as directory:
            file = Path(directory) / "results.json"
            file.write_text(json.dumps(data))
            return list(parser.rows(file))

    def test_units_and_normalization(self):
        row, = self.parse(self.fixture)
        self.assertEqual(row["us_per_call"], 2000)
        self.assertEqual(row["ms_per_map"], 4)
        self.assertEqual(row["logical_bytes"], 4000)

    def test_count_only_and_zero_calls_are_not_timings(self):
        lane = self.fixture["phase_profile"]["fft_lane"]
        lane["stage_timers"] = "off"
        stage = lane["per_rank"][0]["buckets"][0]["stages"][0]
        stage["seconds"] = 0
        row, = self.parse(self.fixture)
        self.assertIsNone(row["seconds"])
        self.assertIsNone(row["us_per_call"])
        stage["calls"] = 0
        row, = self.parse(self.fixture)
        self.assertEqual(row["status"], "unusable-zero-calls")

    def test_invalid_or_missing_data_rejected(self):
        for mode in ["missing", "disabled", "rank", "negative", "nan", "threads", "timer"]:
            data = copy.deepcopy(self.fixture)
            lane = data["phase_profile"]["fft_lane"]
            if mode == "missing": lane.clear()
            if mode == "disabled": lane["enabled"] = False
            if mode == "threads": data["threading"] = {"threads": 2}
            if mode == "timer": lane["stage_timers"] = "unknown"
            if mode == "rank": lane["per_rank"][0]["rank"] = 2
            if mode in ("negative", "nan"):
                lane["per_rank"][0]["buckets"][0]["stages"][0]["seconds"] = -1 if mode == "negative" else float("nan")
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.parse(data)


if __name__ == "__main__":
    unittest.main()
