"""Unit-factor guard: python3 -I scripts/tests/test_no_hardcoded_units.py.

All workspace Rust sources (including build scripts and inline tests) are scanned.
Only explicit reference-test files below are exempt; new references must be placed
under tests/ and reviewed here. Documentation comments in docs/ are outside the
Rust-source scan. Build output, vendored dependencies and git metadata are not
workspace source. This is a literal-prefix guard, not a Rust constant evaluator.
"""
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE_TESTS = {"crates/dft-units/tests/units.rs"}
EXCLUDED_DIRS = {"target", "vendor", ".git", ".agents", ".codex"}
# Match decimal Rust literals, including digit separators and scientific notation.
LITERAL = re.compile(r"(?<![\w.])(?:\d[\d_]*(?:\.[\d_]*)?|\.\d[\d_]*)(?:[eE][+-]?\d[\d_]*)?")
PREFIXES = ("5291772", "2721138", "1360569", "1602176", "4359744", "2941", "2942")


def violations(root):
    found = []
    for path in sorted(root.rglob("*.rs")):
        relative = path.relative_to(root)
        if EXCLUDED_DIRS.intersection(relative.parts) or relative.as_posix() in REFERENCE_TESTS:
            continue
        for number, line in enumerate(path.read_text().splitlines(), 1):
            for match in LITERAL.finditer(line):
                mantissa = re.split("[eE]", match.group().replace("_", ""))[0]
                digits = mantissa.replace(".", "").lstrip("0")
                if digits.startswith(PREFIXES):
                    found.append(f"{relative}:{number}: forbidden unit-factor literal {match.group()}")
    return found


class NoHardcodedUnits(unittest.TestCase):
    def test_workspace(self):
        self.assertEqual(violations(ROOT), [])

    def test_planted_literals_in_temporary_copy(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            # Copy workspace Rust files only, so no generated sources or metadata.
            for source in ROOT.rglob("*.rs"):
                rel = source.relative_to(ROOT)
                if EXCLUDED_DIRS.intersection(rel.parts):
                    continue
                target = root / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
            self.assertEqual(violations(root), [])
            for directory in ("crates/dft-units/src", "crates/dft-cli/src"):
                target = root / directory / "planted.rs"
                target.parent.mkdir(parents=True, exist_ok=True)
                for literal in ("0.529177210903", "27.211_386", "13.60569", "1.602176e-19", "4.359744e-18", "2.941e4", "29_421.015", "5.291772e-11"):
                    target.write_text(f"const PLANTED: f64 = {literal};\n")
                    result = subprocess.run([sys.executable, "-I", str(pathlib.Path(__file__).resolve()), "--root", str(root)], capture_output=True, text=True)
                    self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                    self.assertIn(f"{directory}/planted.rs:1:", result.stdout)
                target.unlink()
            self.assertEqual(violations(root), [])


if __name__ == "__main__":
    root = pathlib.Path(sys.argv[2]) if len(sys.argv) == 3 and sys.argv[1] == "--root" else ROOT
    errors = violations(root)
    print("\n".join(errors) if errors else "PASS no hardcoded unit factors")
    sys.exit(bool(errors))
