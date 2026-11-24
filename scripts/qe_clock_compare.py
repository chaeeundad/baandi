#!/usr/bin/env python3
"""Rank-zero QE/BAANDI clock diagnostics. Never a headline speed comparison."""
import argparse
import json
import math
import re
import sys
from pathlib import Path

CLOCK = re.compile(
    r"^\s*(\S+)\s*:\s*(.*?)\s+CPU\s+(.*?)\s+WALL\s*\(\s*(\d+)\s+calls?\s*\)",
    re.MULTILINE,
)
DURATION = re.compile(r"\s*(\d+(?:\.\d*)?|\.\d+)\s*([dhms])")
ALIASES = {"calbec": ("h_psi:calbec", "calbec")}


def seconds(value):
    """QE prints seconds or composite durations, e.g. 1h 2m 3.40s."""
    pos, total = 0, 0.0
    for match in DURATION.finditer(value):
        if value[pos:match.start()].strip():
            raise ValueError(f"invalid QE duration: {value!r}")
        total += float(match[1]) * {"d": 86400, "h": 3600, "m": 60, "s": 1}[match[2]]
        pos = match.end()
    if not pos or value[pos:].strip() or not math.isfinite(total):
        raise ValueError(f"invalid QE duration: {value!r}")
    return total


def parse_qe(text):
    runs = len(re.findall(r"Program\s+PWSCF", text))
    if runs > 1 or (runs == 1 and "JOB DONE." not in text):
        raise ValueError("supply one complete QE run, not concatenated or truncated output")
    rows = {}
    for match in CLOCK.finditer(text):
        name, cpu, wall, calls = match.groups()
        row = {"cpu_s": seconds(cpu), "wall_s": seconds(wall), "calls": int(calls), "rank": 0}
        if name in rows and rows[name] != row:
            raise ValueError(f"conflicting QE clock {name}; supply one complete run")
        rows[name] = row
    if not rows:
        raise ValueError("no QE clock rows with wall time and call counts")
    return rows


def per_call(row):
    if row is None or row.get("wall_s") is None or not row.get("calls"):
        return None
    return row["wall_s"] / row["calls"]


def compare(result, qe):
    reference = result.get("reference_scf", {})
    rows = reference.get("qe_clock")
    if not rows:
        raise ValueError("missing rank-zero reference_scf.qe_clock; enable davidson_stages")
    if reference.get("diagnostics", {}).get("level") != "full":
        raise ValueError("BAANDI diagnostics must be full")
    names, out = set(), []
    for row in rows:
        name = row["name"]
        if name in names or row["rank"] != 0:
            raise ValueError(f"duplicate or non-root BAANDI clock: {name}")
        names.add(name)
        wall, calls = row["wall_s"], row["calls"]
        if (wall is None) != (calls is None):
            raise ValueError(f"incomplete BAANDI clock: {name}")
        if wall is not None and (not math.isfinite(wall) or wall < 0 or type(calls) is not int or calls < 0):
            raise ValueError(f"invalid BAANDI clock: {name}")
        candidates = ALIASES.get(name, (name,))
        qe_name = next((label for label in candidates if label in qe), None)
        other = qe.get(qe_name)
        boundary = row["boundary"]
        if boundary not in ("exact", "unavailable", "baandi_only"):
            raise ValueError(f"unknown boundary status: {boundary}")
        # A global calbec clock may also include operations outside H.
        if name == "calbec" and qe_name == "calbec":
            boundary = "inexact: QE calbec is not H-scoped"
        if boundary == "baandi_only":
            other, qe_name = None, None
        if boundary == "exact" and other is None:
            boundary = "QE clock missing"
        out.append({"name": name, "boundary": boundary, "qe_label": qe_name,
                    "baandi_wall_s": wall, "baandi_calls": calls,
                    "baandi_s_per_call": per_call(row),
                    "qe_wall_s": other["wall_s"] if other else None,
                    "qe_calls": other["calls"] if other else None,
                    "qe_s_per_call": per_call(other)})
    # Global h_psi/rdiaghg include startup (and possibly other parents).
    # Never subtract them from regterg to fabricate a Davidson remainder.
    return out


def markdown(rows):
    def fmt(value):
        return "unavailable" if value is None else f"{value:.9g}" if isinstance(value, float) else str(value)
    keys = ("name", "boundary", "baandi_wall_s", "baandi_calls", "baandi_s_per_call",
            "qe_label", "qe_wall_s", "qe_calls", "qe_s_per_call")
    lines = ["Diagnosis only; rank 0 wall seconds. Inclusive rows overlap. No speed ratios.",
             "QE printed precision limits per-call precision; zero printed time is not free work.", "",
             "| " + " | ".join(keys) + " |", "| " + " | ".join("---" for _ in keys) + " |"]
    lines.extend("| " + " | ".join(fmt(row[key]) for key in keys) + " |" for row in rows)
    return "\n".join(lines) + "\n"


def self_test():
    import unittest

    class ParserTests(unittest.TestCase):
        def test_committed_qe_output(self):
            path = Path(__file__).resolve().parents[1] / "tests/reference/qe/i25-si8/si8-np1/k2-mix/stdout.txt"
            text = path.read_text()
            rows = parse_qe(text)
            self.assertEqual(rows["h_psi"]["calls"], 7)
            self.assertEqual(rows["h_psi:calbec"]["calls"], 7)
            self.assertEqual(rows["fftw"]["calls"], 156)
            with self.assertRaises(ValueError):
                parse_qe(text + text)
            with self.assertRaises(ValueError):
                parse_qe(text.replace("JOB DONE.", ""))

        def test_durations_and_repeated_sections(self):
            line = " h_psi:calbec : 1m 2.30s CPU 1h 2m 3.40s WALL ( 12 calls)\n"
            rows = parse_qe(line * 2)
            self.assertEqual(rows["h_psi:calbec"]["wall_s"], 3723.4)
            self.assertEqual(rows["h_psi:calbec"]["calls"], 12)
            with self.assertRaises(ValueError):
                parse_qe(line + line.replace("12 calls", "13 calls"))
            with self.assertRaises(ValueError):
                seconds("garbage 1s")

        def test_unavailable_zero_and_scoped_alias(self):
            rows = [{"name": n, "wall_s": w, "calls": c, "rank": 0, "boundary": b}
                    for n, w, c, b in [("calbec", 2., 4, "exact"),
                                       ("fftw", None, None, "unavailable"),
                                       ("output", 1., 1, "baandi_only"),
                                       ("regterg", 0., 0, "exact")]]
            result = {"reference_scf": {"qe_clock": rows, "diagnostics": {"level": "full"}}}
            qe = parse_qe("calbec : 0.10s CPU 0.20s WALL ( 2 calls)\nfftw : 0.00s CPU 0.00s WALL ( 4 calls)")
            out = compare(result, qe)
            self.assertTrue(out[0]["boundary"].startswith("inexact"))
            self.assertEqual(out[0]["baandi_s_per_call"], .5)
            self.assertIsNone(out[1]["baandi_s_per_call"])
            self.assertEqual(out[1]["qe_s_per_call"], 0.)
            self.assertIsNone(out[2]["qe_wall_s"])
            self.assertIsNone(out[3]["baandi_s_per_call"])
            qe["h_psi:calbec"] = qe["calbec"]
            self.assertEqual(compare(result, qe)[0]["boundary"], "exact")
            rows[0]["rank"] = 1
            with self.assertRaises(ValueError):
                compare(result, qe)

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ParserTests)
    return unittest.TextTestRunner().run(suite).wasSuccessful()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baandi", nargs="?", type=Path, help="rank-zero results.json")
    parser.add_argument("qe", nargs="?", type=Path, help="one complete QE stdout")
    parser.add_argument("--json", action="store_true", help="emit JSON rather than Markdown")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return 0 if self_test() else 1
    if not args.baandi or not args.qe:
        parser.error("baandi and qe paths are required")
    try:
        rows = compare(json.loads(args.baandi.read_text()), parse_qe(args.qe.read_text()))
        print(json.dumps({"purpose": "diagnosis_only", "statistic": "rank_0_wall_seconds",
                          "rows": rows}, indent=2, allow_nan=False) if args.json else markdown(rows), end="\n")
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"qe_clock_compare: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
