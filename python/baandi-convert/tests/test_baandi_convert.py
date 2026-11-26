"""Tests for baandi-convert (standard library unittest; Python 3.9+).

Run: python3 -I python/baandi-convert/tests/test_baandi_convert.py -v
ASE-dependent tests are skipped with a printed note when ASE is missing.
No test touches the network: urlopen is replaced by a function that fails.
"""
import contextlib
import io
import json
import math
import os
import re
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import baandi_convert as bc  # noqa: E402

FIXTURES = ROOT / "crates/dft-cli/tests/fixtures/structure_convert"

try:
    import ase  # noqa: F401

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False
    print("note: ASE is not installed; ASE reader tests are skipped (pip install ase)", file=sys.stderr)


def no_network(*args, **kwargs):
    raise AssertionError("test attempted network access")


def parse_structure(text):
    """Minimal independent reader of the emitted subset (no tomllib on 3.9)."""
    rows = re.search(r"lattice_vectors = \[\n(.*?)\n\]", text, re.S).group(1)
    lattice = [[float(x) for x in r.strip(" [],").split(",")] for r in rows.splitlines()]
    species = re.findall(r'^species = "([^"]+)"$', text, re.M)
    fractional = [
        [float(x) for x in m.split(",")] for m in re.findall(r"^fractional = \[(.*)\]$", text, re.M)
    ]
    unit = re.search(r'^lattice_unit = "(\w+)"$', text, re.M).group(1)
    return unit, lattice, list(zip(species, fractional))


def canonical_sites(sites):
    return sorted((s, tuple(round(x % 1.0, 9) % 1.0 for x in f)) for s, f in sites)


SI_PMG = {
    "lattice": {"matrix": [[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]},
    "sites": [
        {"abc": [0.0, 0.0, 0.0], "species": [{"element": "Si", "occu": 1}]},
        {"abc": [0.25, 0.25, 0.25], "species": [{"element": "Si", "occu": 1}]},
    ],
}


class Writer(unittest.TestCase):
    def test_float_format_round_trips_and_is_a_toml_float(self):
        for x in [0.0, -0.0, 1.0, 2.715, 5.43 / 3, 1e-20, -123456.789, 0.1 + 0.2]:
            text = bc.toml_float(x)
            self.assertEqual(float(text), x + 0.0)
            self.assertTrue(any(c in text for c in ".eE"), text)
            self.assertFalse(text.startswith("-0") and float(text) == 0.0)
        self.assertEqual(bc.toml_float(2.715), "2.7149999999999999")
        with self.assertRaises(bc.ConversionError):
            bc.toml_float(math.nan)

    def test_units_and_skeleton(self):
        s = bc.structure_from_pmg_dict(SI_PMG, "test")
        s.provenance = {"source_kind": "test"}
        unit, lattice, sites = parse_structure(bc.to_toml(s, units="bohr"))
        self.assertEqual(unit, "bohr")
        self.assertAlmostEqual(lattice[0][1], 2.715 / bc.BOHR_IN_ANGSTROM, places=14)
        text = bc.to_toml(s, skeleton=True)
        for needle in ["[pseudopotentials]", 'Si = "/absolute/path/to/pseudo/Si.UPF"', "[basis]", "[xc]", "[scf]"]:
            self.assertIn(needle, text)
        self.assertIn('name = "Si2"', text)


class Rejections(unittest.TestCase):
    def test_partial_occupancy_mixed_sites_and_overlap(self):
        partial = json.loads(json.dumps(SI_PMG))
        partial["sites"][1]["species"][0]["occu"] = 0.5
        with self.assertRaisesRegex(bc.ConversionError, "partial occupancy"):
            bc.structure_from_pmg_dict(partial, "x")
        mixed = json.loads(json.dumps(SI_PMG))
        mixed["sites"][1]["species"].append({"element": "Ge", "occu": 0.5})
        with self.assertRaisesRegex(bc.ConversionError, "disordered"):
            bc.structure_from_pmg_dict(mixed, "x")
        overlap = json.loads(json.dumps(SI_PMG))
        overlap["sites"][1]["abc"] = [0.999999, 0.0, 0.0]
        with self.assertRaisesRegex(bc.ConversionError, "overlapping sites"):
            bc.structure_from_pmg_dict(overlap, "x")

    def test_missing_cell(self):
        flat = json.loads(json.dumps(SI_PMG))
        flat["lattice"]["matrix"][2] = [0.0, 0.0, 0.0]
        with self.assertRaisesRegex(bc.ConversionError, "missing cell"):
            bc.structure_from_pmg_dict(flat, "x")
        with self.assertRaisesRegex(bc.ConversionError, "missing cell"):
            bc.structure_from_pmg_dict({"sites": []}, "x")


class MaterialsProjectOffline(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="baandi-convert-"))
        self.addCleanup(shutil.rmtree, self.tmp)
        patcher = mock.patch("urllib.request.urlopen", no_network)
        patcher.start()
        self.addCleanup(patcher.stop)
        home = mock.patch.dict(os.environ, {"HOME": str(self.tmp / "home")})
        home.start()
        self.addCleanup(home.stop)
        os.environ.pop("MP_API_KEY", None)
        os.environ.pop("BAANDI_CACHE_DIR", None)

    def run_main(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = bc.main(argv)
        return code, out.getvalue(), err.getvalue()

    def test_cached_record_converts_offline_and_reproduces_the_committed_fixture(self):
        cache = self.tmp / "cache"
        shutil.copytree(FIXTURES / "mp", cache / "mp")
        os.environ["MP_API_KEY"] = "secret-test-key-123"
        with mock.patch.dict(os.environ, {"SOURCE_DATE_EPOCH": "1795651200", "BAANDI_CACHE_DIR": str(cache)}):
            code, out, err = self.run_main(["--mp", "mp-0"])
        self.assertEqual(code, 0, err)
        self.assertNotIn("secret-test-key-123", out + err)
        self.assertEqual(out, (FIXTURES / "si_mp0_converted.toml").read_text())
        self.assertIn('source_sha256 = "' + bc.sha256_file(cache / "mp/mp-0.json") + '"', out)

    def test_missing_cache_without_key_is_a_clear_error(self):
        code, out, err = self.run_main(["--mp", "mp-149", "--cache-dir", str(self.tmp / "empty")])
        self.assertEqual(code, 1)
        self.assertIn("no cached record", err)
        self.assertIn("MP_API_KEY", err)
        self.assertIn("mp_api_key", err)

    def test_key_file_is_read_and_never_echoed_on_failure(self):
        keyfile = self.tmp / "home/.config/baandi/mp_api_key"
        keyfile.parent.mkdir(parents=True)
        keyfile.write_text("secret-file-key-456\n")
        self.assertEqual(bc.api_key(), "secret-file-key-456")

        def failing(request, timeout):
            raise OSError("refused for key secret-file-key-456")

        with self.assertRaises(bc.ConversionError) as ctx:
            bc.fetch_mp("mp-149", cache_dir=self.tmp / "c", opener=failing)
        self.assertNotIn("secret-file-key-456", str(ctx.exception))
        self.assertIn("***", str(ctx.exception))
        self.assertFalse((self.tmp / "c/mp/mp-149.json").exists())

    def test_rest_fetch_writes_a_cache_without_the_key(self):
        responses = {
            "summary": {"data": [{"material_id": "mp-0", "structure": SI_PMG}], "meta": {}},
            "heartbeat": {"db_version": "test-db-version"},
        }
        seen = []

        def fake(request, timeout):
            seen.append(request.get_header("X-api-key"))
            body = responses["heartbeat" if "heartbeat" in request.full_url else "summary"]
            return _Response(body)

        record, path = bc.fetch_mp("mp-0", cache_dir=self.tmp / "c", key="k-789", opener=fake)
        self.assertEqual(seen, ["k-789", "k-789"])
        self.assertEqual(record["database_version"], "test-db-version")
        self.assertNotIn("k-789", path.read_text())
        # A second call is served from the cache with no key and no network.
        again, _ = bc.fetch_mp("mp-0", cache_dir=self.tmp / "c")
        self.assertEqual(again["retrieved_utc"], record["retrieved_utc"])

    def test_invalid_id_and_corrupt_cache(self):
        with self.assertRaisesRegex(bc.ConversionError, "not a Materials Project id"):
            bc.fetch_mp("../etc/passwd", cache_dir=self.tmp)
        bad = self.tmp / "mp/mp-1.json"
        bad.parent.mkdir(parents=True)
        bad.write_text("{}")
        with self.assertRaisesRegex(bc.ConversionError, "--refresh"):
            bc.fetch_mp("mp-1", cache_dir=self.tmp)


class _Response:
    def __init__(self, body):
        self.body = json.dumps(body).encode()

    def read(self):
        return self.body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class MissingReader(unittest.TestCase):
    def test_missing_ase_gives_install_advice(self):
        tmp = Path(tempfile.mkdtemp(prefix="baandi-convert-"))
        self.addCleanup(shutil.rmtree, tmp)
        path = tmp / "Si.cif"
        path.write_text("data_x\n")
        with mock.patch.dict(sys.modules, {"ase": None, "ase.io": None}):
            with self.assertRaisesRegex(bc.ConversionError, "pip install ase"):
                bc.read_file(path)


@unittest.skipUnless(HAVE_ASE, "ASE is not installed")
class AseFormats(unittest.TestCase):
    """Si written by ASE's own writers in four formats converts to one structure."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="baandi-convert-ase-"))
        self.addCleanup(shutil.rmtree, self.tmp)

    def test_cif_poscar_extxyz_xsf_agree(self):
        from ase.build import bulk
        from ase.io import write

        si = bulk("Si", "diamond", a=5.43, cubic=True)
        files = {
            "Si.cif": "cif",
            "POSCAR": "vasp",
            "Si.extxyz": "extxyz",
            "Si.xsf": "xsf",
        }
        parsed = {}
        for name, fmt in files.items():
            write(str(self.tmp / name), si, format=fmt)
            s = bc.read_file(self.tmp / name, fmt=fmt)
            unit, lattice, sites = parse_structure(bc.to_toml(s))
            self.assertEqual(unit, "angstrom")
            parsed[name] = (lattice, canonical_sites(sites), s)
        ref_lattice, ref_sites, _ = parsed["POSCAR"]
        for name, (lattice, sites, s) in parsed.items():
            for r, q in zip(lattice, ref_lattice):
                for a, b in zip(r, q):
                    self.assertLess(abs(a - b), 1e-10, name)
            self.assertEqual([x[0] for x in sites], [x[0] for x in ref_sites], name)
            for (_, f), (_, g) in zip(sites, ref_sites):
                for a, b in zip(f, g):
                    d = abs(a - b)
                    self.assertLess(min(d, 1 - d), 1e-10, name)
            self.assertEqual(s.provenance["reader"], "ase")
            self.assertEqual(len(s.provenance["source_sha256"]), 64)
        self.assertEqual(len(ref_sites), 8)

    def test_cif_partial_occupancy_is_rejected(self):
        path = self.tmp / "disordered.cif"
        path.write_text(
            "data_SiGe\n_cell_length_a 5.5\n_cell_length_b 5.5\n_cell_length_c 5.5\n"
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
            "_symmetry_space_group_name_H-M 'P 1'\n"
            "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n"
            "_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
            "Si1 Si 0 0 0 0.5\nGe1 Ge 0 0 0 0.5\nSi2 Si 0.5 0.5 0.5 1\n"
        )
        with self.assertRaisesRegex(bc.ConversionError, "occupancy|disorder"):
            bc.read_file(path)

    def test_nonperiodic_needs_assume_periodic(self):
        from ase import Atoms
        from ase.io import write

        path = self.tmp / "mol.xyz"
        write(str(path), Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]]), format="xyz")
        with self.assertRaisesRegex(bc.ConversionError, "missing cell|periodic"):
            bc.read_file(path)
        boxed = self.tmp / "boxed.extxyz"
        write(
            str(boxed),
            Atoms("H2", positions=[[5, 5, 5], [5, 5, 5.74]], cell=[10, 10, 10], pbc=False),
            format="extxyz",
        )
        with self.assertRaisesRegex(bc.ConversionError, "--assume-periodic"):
            bc.read_file(boxed)
        s = bc.read_file(boxed, assume_periodic=True)
        self.assertEqual(s.provenance["pbc_in_file"], [False, False, False])


if __name__ == "__main__":
    unittest.main()
