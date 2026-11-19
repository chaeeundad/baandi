#!/usr/bin/env python3
"""Fixture alias regression tests; no UPF content or DFT run needed."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("handoff", ROOT / "scripts/ladder_handoff.py")
handoff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(handoff)


class FixtureAliases(unittest.TestCase):
    def test_si64_alias_loads_shipped_gamma_fixture(self):
        with patch.object(handoff, "sha", return_value="test-only-no-upf-read"):
            alias, _ = handoff.fixture("si64", 601, 0., Path("/unused"))
            explicit, _ = handoff.fixture("si64-gamma", 601, 0., Path("/unused"))
        self.assertEqual(alias, explicit)
        self.assertEqual(sum(line.startswith("atom Si ") for line in alias.splitlines()), 64)
        self.assertIn("cutoff 15.0", alias)

    def test_all_size_aliases_match_explicit_fixture(self):
        with patch.object(handoff, "sha", return_value="test-only-no-upf-read"):
            for n in (8, 16, 32, 64):
                with self.subTest(n=n):
                    alias, _ = handoff.fixture(f"si{n}", 603, .001, Path("/unused"))
                    explicit, _ = handoff.fixture(f"si{n}-gamma", 603, .001, Path("/unused"))
                    self.assertEqual(alias, explicit)
                    self.assertEqual(sum(line.startswith("atom Si ") for line in alias.splitlines()), n)


class CrossedPairOrder(unittest.TestCase):
    def test_full_pairs_alternate_and_keep_aa_separate(self):
        spec = importlib.util.spec_from_file_location("pairs", ROOT / "scripts/ladder_2x2_pairs.py")
        pairs = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(pairs)
        arms, manifest = pairs.arms_for([1.2345678901234567e-8, 2.9876543210987654e-6], 3)
        self.assertEqual(len(arms.split(",")), len(manifest))
        for pair in range(3):
            rows = [r for r in manifest if r["pair"] == pair]
            self.assertEqual(sum(r["kind"] == "aa" for r in rows), 2)
            cells = [(r["state"], r["seed"]) for r in rows if r["kind"] == "crossed"]
            expected = [(0, 0), (0, 1), (1, 0), (1, 1)]
            self.assertEqual(cells, expected[::-1] if pair % 2 else expected)
            self.assertTrue(all(r["discard_warmup"] == (pair == 0) for r in rows))
        self.assertEqual(float(arms.split(",")[1].split(":")[3]), 1.2345678901234567e-8)


if __name__ == "__main__":
    unittest.main()
