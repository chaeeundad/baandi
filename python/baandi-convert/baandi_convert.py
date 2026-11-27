#!/usr/bin/env python3
"""baandi-convert: write a BAANDI ``[system]`` block (or, with ``--compact``, a
``[structure]`` block) from a structure file or a Materials Project id
(workstation-side tool, Python 3.9+).

The BAANDI CLI reads only POSCAR and extended XYZ files natively
(``[structure] file``, Rust readers, no Python) and never contacts a database:
compute nodes may have no Python, ASE or network. This tool runs where those
exist (CIF with symmetry expansion, Materials Project, any ASE format) and
writes the explicit TOML form (lattice rows plus fractional
positions, 17 significant digits) with a ``[system.provenance]`` table that the
CLI copies into results.json without using it for the physics.

Readers: ASE (``ase.io.read``) by default, pymatgen with ``--reader pymatgen``.
Materials Project: ``mp-api`` when installed, else plain HTTPS (standard
library). The API key comes from ``MP_API_KEY`` or
``~/.config/baandi/mp_api_key`` and is never printed or written. Fetched
records are cached (``BAANDI_CACHE_DIR`` or ``~/.cache/baandi``) under
``mp/<id>.json``; later conversions of the same id are offline.

The functions ``read_file``, ``fetch_mp``, ``structure_from_mp_record`` and
``to_toml`` are the programmatic API (for example for an interactive shell).
"""

import argparse
import datetime
import hashlib
import json
import math
import os
import re
import sys
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

CONVERTER = "baandi-convert"
CONVERTER_VERSION = "1"
# CODATA 2018 Bohr radius in Angstrom (only used for --units bohr).
BOHR_IN_ANGSTROM = 0.529177210903
MP_ENDPOINT = "https://api.materialsproject.org"
CACHE_SCHEMA = "baandi-mp-cache-1"
OCCUPANCY_TOLERANCE = 1e-6
# Two sites closer than this (Angstrom) are treated as a disordered/duplicated site.
MIN_SITE_DISTANCE_ANGSTROM = 0.5
MIN_VOLUME_ANGSTROM3 = 1e-6
_MP_ID = re.compile(r"^(mp|mvc)-[0-9]+$")


class ConversionError(Exception):
    """A structure that cannot be written as a BAANDI periodic input."""


class Structure:
    """Periodic structure: lattice rows (Angstrom), symbols, fractional positions."""

    def __init__(self, lattice, symbols, fractional, provenance, fixed=None, magmoms=None):
        self.lattice = [[float(x) for x in row] for row in lattice]
        self.symbols = [str(s) for s in symbols]
        self.fractional = [[float(x) for x in p] for p in fractional]
        self.provenance = dict(provenance)
        # Per-atom Cartesian freeze masks (None = free) and initial moments
        # (mu_B); written only in the compact form (--compact).
        self.fixed = fixed
        self.magmoms = magmoms

    def formula(self):
        counts = {}
        for s in self.symbols:
            counts[s] = counts.get(s, 0) + 1
        return "".join(f"{s}{n}" for s, n in counts.items())


# ---------------------------------------------------------------- utilities


def utc_now():
    """UTC time as ISO 8601; honours SOURCE_DATE_EPOCH for reproducible output."""
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    moment = (
        datetime.datetime.fromtimestamp(int(epoch), datetime.timezone.utc)
        if epoch
        else datetime.datetime.now(datetime.timezone.utc)
    )
    return moment.strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def toml_float(x):
    """17 significant digits (exact round trip) in a form TOML reads as a float."""
    x = float(x) + 0.0  # turns -0.0 into 0.0
    if not math.isfinite(x):
        raise ConversionError(f"non-finite number {x!r} in structure")
    text = "%.17g" % x
    if not any(c in text for c in ".eE"):
        text += ".0"
    return text


def toml_string(s):
    # JSON string escapes are valid TOML basic-string escapes.
    return json.dumps(str(s), ensure_ascii=False)


def toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return toml_float(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    return toml_string(v)


def _det(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def check_structure(s):
    """Reject inputs BAANDI cannot represent; messages name the cause."""
    if len(s.lattice) != 3 or any(len(r) != 3 for r in s.lattice):
        raise ConversionError("missing cell: three lattice vectors are required")
    if not all(math.isfinite(x) for r in s.lattice for x in r):
        raise ConversionError("missing cell: lattice vectors are not finite")
    volume = abs(_det(s.lattice))
    if volume < MIN_VOLUME_ANGSTROM3:
        raise ConversionError(
            "missing cell: the structure has no three-dimensional cell "
            f"(volume {volume:.3g} A^3); BAANDI is a periodic plane-wave code"
        )
    if not s.symbols:
        raise ConversionError("the structure contains no atoms")
    if len(s.symbols) != len(s.fractional):
        raise ConversionError("internal: symbol and position counts differ")
    for sym in s.symbols:
        if not re.fullmatch(r"[A-Z][a-z]?", sym):
            raise ConversionError(
                f"site label {sym!r} is not a chemical element symbol "
                "(vacancies, dummy atoms and mixed sites are not supported)"
            )
    for p in s.fractional:
        if not all(math.isfinite(x) for x in p):
            raise ConversionError("atomic positions are not finite")
    # Minimum-image distance over neighbouring images; flags duplicated sites
    # that disordered files often contain (two species on one position).
    images = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for a in range(len(s.fractional)):
        for b in range(a + 1, len(s.fractional)):
            d = [s.fractional[a][c] - s.fractional[b][c] for c in range(3)]
            d = [x - round(x) for x in d]
            best = math.inf
            for shift in images:
                f = [d[c] + shift[c] for c in range(3)]
                cart = [sum(f[i] * s.lattice[i][c] for i in range(3)) for c in range(3)]
                best = min(best, math.sqrt(sum(x * x for x in cart)))
            if best < MIN_SITE_DISTANCE_ANGSTROM:
                raise ConversionError(
                    f"sites {a + 1} ({s.symbols[a]}) and {b + 1} ({s.symbols[b]}) are "
                    f"{best:.3g} A apart: overlapping sites usually mean a disordered "
                    "structure; resolve the disorder (choose one ordered configuration) first"
                )


# ---------------------------------------------------------------- readers


def _require_ase():
    try:
        import ase  # noqa: F401
        import ase.io  # noqa: F401
    except ImportError as e:
        raise ConversionError(
            "ASE is required to read structure files: pip install ase "
            f"(import failed: {e})"
        ) from None
    return sys.modules["ase"]


def _cif_partial_occupancy(path):
    """Occupancies below one in any CIF data block, via ASE's CIF tokenizer."""
    from ase.io import cif

    found = []
    for block in cif.parse_cif(str(path)):
        if isinstance(block, tuple):  # older ASE: (name, tags)
            tags = block[1]
            values = tags.get("_atom_site_occupancy")
        else:
            values = block.get("_atom_site_occupancy")
        for v in values or []:
            text = str(v).split("(")[0].strip()
            if text in ("", "?", "."):
                continue
            try:
                occ = float(text)
            except ValueError:
                continue
            if occ < 1.0 - OCCUPANCY_TOLERANCE:
                found.append(occ)
    return found


def read_file(path, fmt=None, index=-1, reader="ase", assume_periodic=False):
    """Read one structure from a file with ASE (default) or pymatgen."""
    path = Path(path)
    if not path.is_file():
        raise ConversionError(f"{path}: no such file")
    if reader == "pymatgen":
        return _read_pymatgen(path, fmt)
    if reader != "ase":
        raise ConversionError(f"unknown reader {reader!r}; expected ase or pymatgen")
    ase = _require_ase()
    import ase.io

    detected = fmt
    if detected is None:
        try:
            from ase.io.formats import filetype

            detected = filetype(str(path))
        except Exception:  # detection failure is reported by ase.io.read below
            detected = None
    if detected == "cif":
        occupancies = _cif_partial_occupancy(path)
        if occupancies:
            raise ConversionError(
                f"{path}: partial site occupancy {min(occupancies):g} < 1 "
                "(disordered structure); BAANDI needs one ordered configuration"
            )
    try:
        atoms = ase.io.read(str(path), index=index, format=fmt)
    except Exception as e:
        raise ConversionError(f"{path}: ASE could not read the file: {e}") from None
    if isinstance(atoms, list):
        raise ConversionError(f"{path}: --index must select a single structure")
    info_occ = atoms.info.get("occupancy") if hasattr(atoms, "info") else None
    if isinstance(info_occ, dict):
        for site in info_occ.values():
            if isinstance(site, dict) and (
                len(site) > 1 or any(float(o) < 1.0 - OCCUPANCY_TOLERANCE for o in site.values())
            ):
                raise ConversionError(
                    f"{path}: partial or mixed site occupancy (disordered structure); "
                    "BAANDI needs one ordered configuration"
                )
    pbc = [bool(x) for x in atoms.pbc]
    if not all(pbc) and not assume_periodic:
        raise ConversionError(
            f"{path}: periodic boundary conditions are {pbc}; BAANDI treats every cell "
            "as periodic in three directions. Give the file a full cell with enough "
            "vacuum and pass --assume-periodic to accept it as a periodic box"
        )
    lattice = [[float(x) for x in row] for row in atoms.cell[:]]
    try:
        fixed = _ase_fixed(atoms)
    except ConversionError as e:  # raised only if the compact form is written
        fixed = e
    magmoms = None
    if atoms.has("initial_magmoms"):
        magmoms = [float(m) for m in atoms.get_initial_magnetic_moments()]
        if not any(magmoms):
            magmoms = None
    provenance = {
        "source_kind": "file",
        "source": str(path),
        "source_path": str(path.resolve()),
        "source_sha256": sha256_file(path),
        "source_format": detected if detected is not None else "unknown",
        "source_index": index,
        "reader": "ase",
        "reader_version": ase.__version__,
        "pbc_in_file": pbc,
    }
    s = Structure(
        lattice,
        atoms.get_chemical_symbols(),
        atoms.get_scaled_positions(wrap=False).tolist(),
        provenance,
        fixed,
        magmoms,
    )
    check_structure(s)
    return s


def _reciprocal_axis(cell, i):
    """Cartesian axis parallel to a_j x a_k (the direction frozen with fractional i), or None."""
    a, b = cell[(i + 1) % 3], cell[(i + 2) % 3]
    r = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    norm = math.sqrt(sum(x * x for x in r))
    for c in range(3):
        if all(o == c or abs(r[o]) <= 1e-12 * norm for o in range(3)):
            return c
    return None


def _ase_fixed(atoms):
    """Per-atom Cartesian freeze masks from ASE FixAtoms/FixCartesian/FixScaled; None if unconstrained."""
    constraints = getattr(atoms, "constraints", None) or []
    if not constraints:
        return None
    from ase.constraints import FixAtoms, FixCartesian, FixScaled

    cell = [[float(x) for x in row] for row in atoms.cell[:]]
    masks = [[False] * 3 for _ in range(len(atoms))]
    for c in constraints:
        if isinstance(c, FixAtoms):
            for i in c.index:
                masks[int(i)] = [True] * 3
        elif isinstance(c, FixCartesian):
            for i in c.index:
                masks[int(i)] = [m or bool(x) for m, x in zip(masks[int(i)], c.mask)]
        elif isinstance(c, FixScaled):
            for i in c.index:
                frozen = [bool(x) for x in c.mask]
                if all(frozen):
                    masks[int(i)] = [True] * 3
                    continue
                for d in range(3):
                    if not frozen[d]:
                        continue
                    axis = _reciprocal_axis(cell, d)
                    if axis is None:
                        raise ConversionError(
                            f"atom {int(i) + 1}: a constraint along lattice direction {d + 1} is not a "
                            "Cartesian axis in this cell; BAANDI constraints are Cartesian components"
                        )
                    masks[int(i)][axis] = True
        else:
            raise ConversionError(
                f"constraint {type(c).__name__} cannot be written as BAANDI fixed flags; remove it"
            )
    return [m if any(m) else None for m in masks]


def _read_pymatgen(path, fmt):
    try:
        import pymatgen.core
        from pymatgen.core import Structure as PmgStructure
    except ImportError as e:
        raise ConversionError(
            f"pymatgen is required for --reader pymatgen: pip install pymatgen ({e})"
        ) from None
    if fmt is not None:
        raise ConversionError("--format is only used by the ASE reader")
    try:
        structure = PmgStructure.from_file(str(path))
    except Exception as e:
        raise ConversionError(f"{path}: pymatgen could not read the file: {e}") from None
    s = structure_from_pmg_dict(structure.as_dict(), str(path))
    s.provenance = {
        "source_kind": "file",
        "source": str(path),
        "source_path": str(Path(path).resolve()),
        "source_sha256": sha256_file(path),
        "reader": "pymatgen",
        "reader_version": getattr(pymatgen.core, "__version__", "unknown"),
    }
    return s


def structure_from_pmg_dict(record, where):
    """Ordered structure from a pymatgen ``Structure.as_dict()`` record."""
    try:
        lattice = record["lattice"]["matrix"]
        sites = record["sites"]
    except (KeyError, TypeError):
        raise ConversionError(f"{where}: record has no lattice matrix or sites (missing cell)") from None
    symbols, fractional = [], []
    for i, site in enumerate(sites):
        species = site.get("species") or []
        if len(species) != 1:
            raise ConversionError(
                f"{where}: site {i + 1} holds {len(species)} species (disordered structure); "
                "BAANDI needs one ordered configuration"
            )
        occ = float(species[0].get("occu", 1.0))
        if occ < 1.0 - OCCUPANCY_TOLERANCE:
            raise ConversionError(
                f"{where}: site {i + 1} has partial occupancy {occ:g} < 1 (disordered "
                "structure); BAANDI needs one ordered configuration"
            )
        element = re.match(r"[A-Z][a-z]?", str(species[0].get("element", "")))
        symbols.append(element.group(0) if element else str(species[0].get("element")))
        fractional.append(site["abc"])
    s = Structure(lattice, symbols, fractional, {})
    check_structure(s)
    return s


# ---------------------------------------------------------------- Materials Project


def cache_root(cache_dir=None):
    if cache_dir is not None:
        return Path(cache_dir)
    if os.environ.get("BAANDI_CACHE_DIR"):
        return Path(os.environ["BAANDI_CACHE_DIR"])
    return Path.home() / ".cache" / "baandi"


def api_key():
    """MP_API_KEY, else ~/.config/baandi/mp_api_key (bare key). Never printed."""
    key = os.environ.get("MP_API_KEY", "").strip()
    if key:
        return key
    path = Path.home() / ".config" / "baandi" / "mp_api_key"
    try:
        key = path.read_text().strip()
    except OSError:
        return None
    return key or None


def _redact(text, key):
    return text.replace(key, "***") if key else text


def _rest_get(url, key, opener):
    request = urllib.request.Request(
        url,
        headers={
            "X-API-KEY": key,
            "Accept": "application/json",
            "User-Agent": f"{CONVERTER}/{CONVERTER_VERSION}",
        },
    )
    with opener(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def _fetch_rest(mp_id, key, opener):
    data = _rest_get(
        f"{MP_ENDPOINT}/materials/summary/?material_ids={mp_id}&_fields=material_id,structure&_limit=1",
        key,
        opener,
    )
    rows = data.get("data") or []
    if not rows or not rows[0].get("structure"):
        raise ConversionError(f"Materials Project returned no structure for {mp_id}")
    version = (data.get("meta") or {}).get("db_version")
    if not version:
        try:
            version = _rest_get(f"{MP_ENDPOINT}/heartbeat", key, opener).get("db_version")
        except Exception:
            version = None
    client = {"name": "urllib", "version": f"python {sys.version.split()[0]}"}
    return rows[0]["structure"], version or "unknown", client


def _fetch_mp_api(mp_id, key):
    import mp_api
    from mp_api.client import MPRester

    with MPRester(key) as rester:
        structure = rester.get_structure_by_material_id(mp_id)
        try:
            version = rester.get_database_version()
        except Exception:
            version = None
    client = {"name": "mp-api", "version": getattr(mp_api, "__version__", "unknown")}
    return structure.as_dict(), version or "unknown", client


def fetch_mp(mp_id, cache_dir=None, refresh=False, key=None, opener=None):
    """Cached Materials Project record for ``mp_id``; returns (record, cache path).

    A cached record is used without any network access unless ``refresh``.
    """
    if not _MP_ID.match(mp_id):
        raise ConversionError(f"{mp_id!r} is not a Materials Project id such as mp-149")
    path = cache_root(cache_dir) / "mp" / f"{mp_id}.json"
    if path.is_file() and not refresh:
        try:
            record = json.loads(path.read_text())
        except ValueError as e:
            raise ConversionError(f"{path}: unreadable cache record ({e}); use --refresh") from None
        if record.get("schema") != CACHE_SCHEMA or record.get("mp_id") != mp_id:
            raise ConversionError(f"{path}: not a {CACHE_SCHEMA} record for {mp_id}; use --refresh")
        return record, path
    key = key or api_key()
    if not key:
        raise ConversionError(
            f"no cached record at {path} and no Materials Project API key: set MP_API_KEY "
            "or write the bare key to ~/.config/baandi/mp_api_key"
        )
    try:
        try:
            import mp_api.client  # noqa: F401

            have_mp_api = opener is None
        except ImportError:
            have_mp_api = False
        if have_mp_api:
            structure, version, client = _fetch_mp_api(mp_id, key)
            endpoint = "mp-api MPRester.get_structure_by_material_id"
        else:
            structure, version, client = _fetch_rest(
                mp_id, key, opener or urllib.request.urlopen
            )
            endpoint = f"{MP_ENDPOINT}/materials/summary/"
    except ConversionError as e:
        raise ConversionError(_redact(str(e), key)) from None
    except Exception as e:  # network, HTTP or client errors
        raise ConversionError(
            _redact(f"Materials Project request for {mp_id} failed: {type(e).__name__}: {e}", key)
        ) from None
    record = {
        "schema": CACHE_SCHEMA,
        "mp_id": mp_id,
        "database_version": str(version),
        "retrieved_utc": utc_now(),
        "client": client,
        "endpoint": endpoint,
        "structure": structure,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=f".{mp_id}.", suffix=".tmp")
    with os.fdopen(fd, "w") as handle:
        json.dump(record, handle, indent=1, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, path)
    return record, path


def structure_from_mp_record(record, path):
    s = structure_from_pmg_dict(record.get("structure") or {}, record.get("mp_id", str(path)))
    client = record.get("client") or {}
    s.provenance = {
        "source_kind": "materials_project",
        "mp_id": record["mp_id"],
        "database_version": record.get("database_version", "unknown"),
        "retrieved_utc": record.get("retrieved_utc", "unknown"),
        "source_sha256": sha256_file(path),
        "reader": client.get("name", "unknown"),
        "reader_version": client.get("version", "unknown"),
    }
    return s


# ---------------------------------------------------------------- writer

SKELETON = """
[pseudopotentials]
# One norm-conserving UPF v2 file per species, matching [xc]; relative paths
# resolve against this file. Inspect each with `baandi inspect-upf FILE`.
{pseudos}

[basis]
# Placeholder cutoff: converge it for your pseudopotentials and property.
ecutwfc_ry = 30.0
# Gamma only; for a Monkhorst-Pack mesh replace kpoints by, for example, kgrid = [4, 4, 4].
kpoints = [[0.0, 0.0, 0.0]]

[xc]
functional = "lda-pz"

[scf]
# All keys are optional; see docs/cli.md. Examples:
# max_iterations = 100
# nbands = 16
"""


def _fixed_flag(mask):
    return "".join(a for a, m in zip("xyz", mask) if m)


def to_compact_toml(structure, units, skeleton, name, provenance):
    """BAANDI TOML text in the compact form: [structure] with one atom per line."""
    if isinstance(structure.fixed, ConversionError):
        raise structure.fixed
    scale = 1.0 if units == "angstrom" else 1.0 / BOHR_IN_ANGSTROM
    lines = [
        f"# Written by {CONVERTER} {CONVERTER_VERSION}; numbers carry 17 significant digits.",
        "# [structure.provenance] is recorded in results.json and never used for the physics.",
    ]
    if units == "bohr":
        lines.append(f"# Lattice converted with 1 Bohr = {BOHR_IN_ANGSTROM} Angstrom (CODATA 2018).")
    lines += [
        "[structure]",
        f"name = {toml_string(name or structure.formula())}",
        f'lattice_unit = "{units}"',
        "# Each row is one lattice vector.",
        "lattice = [",
    ]
    for row in structure.lattice:
        lines.append("  [" + ", ".join(toml_float(x * scale) for x in row) + "],")
    lines += [
        "]",
        'coordinates = "fractional"',
        "# species, three fractional coordinates, optional fixed=<Cartesian components> magmom=<mu_B>",
        "atoms = '''",
    ]
    for i, (sym, f) in enumerate(zip(structure.symbols, structure.fractional)):
        line = f"{sym} " + " ".join(toml_float(x) for x in f)
        if structure.fixed and structure.fixed[i]:
            line += f" fixed={_fixed_flag(structure.fixed[i])}"
        if structure.magmoms and structure.magmoms[i] != 0.0:
            line += f" magmom={toml_float(structure.magmoms[i])}"
        lines.append(line)
    lines += ["'''", "", "[structure.provenance]"]
    for k, v in provenance.items():
        lines.append(f"{k} = {toml_value(v)}")
    text = "\n".join(lines) + "\n"
    if structure.magmoms:
        text += (
            "\n# Nonzero magmom flags need [spin] mode = \"collinear\" with an electron sector;\n"
            "# their sum becomes spin.initial_moment_seed (docs/cli.md, Structure input).\n"
        )
    if skeleton:
        species = list(dict.fromkeys(structure.symbols))
        pseudos = "\n".join(f'{s} = "/absolute/path/to/pseudo/{s}.UPF"' for s in species)
        text += SKELETON.format(pseudos=pseudos)
    return text


def to_toml(structure, units="angstrom", skeleton=False, name=None, compact=False):
    """BAANDI TOML text: [system] (+ provenance and atoms), optionally a full skeleton.

    With ``compact`` the [structure] form is written instead, with per-atom
    fixed and magmom flags when the source carries constraints or moments
    (the [system] form never writes them).
    """
    if units not in ("angstrom", "bohr"):
        raise ConversionError(f"unknown units {units!r}; expected angstrom or bohr")
    scale = 1.0 if units == "angstrom" else 1.0 / BOHR_IN_ANGSTROM
    p = dict(structure.provenance)
    p.update(
        {
            "converter": CONVERTER,
            "converter_version": CONVERTER_VERSION,
            "converted_utc": utc_now(),
            "units_written": units,
        }
    )
    if compact:
        return to_compact_toml(structure, units, skeleton, name, p)
    lines = [
        f"# Written by {CONVERTER} {CONVERTER_VERSION}; numbers carry 17 significant digits.",
        "# [system.provenance] is recorded in results.json and never used for the physics.",
    ]
    if units == "bohr":
        lines.append(f"# Lattice converted with 1 Bohr = {BOHR_IN_ANGSTROM} Angstrom (CODATA 2018).")
    lines += [
        "[system]",
        f"name = {toml_string(name or structure.formula())}",
        f'lattice_unit = "{units}"',
        "# Each row is one lattice vector.",
        "lattice_vectors = [",
    ]
    for row in structure.lattice:
        lines.append("  [" + ", ".join(toml_float(x * scale) for x in row) + "],")
    lines += ["]", "", "[system.provenance]"]
    for k, v in p.items():
        lines.append(f"{k} = {toml_value(v)}")
    for sym, f in zip(structure.symbols, structure.fractional):
        lines += [
            "",
            "[[system.atoms]]",
            f"species = {toml_string(sym)}",
            "fractional = [" + ", ".join(toml_float(x) for x in f) + "]",
        ]
    text = "\n".join(lines) + "\n"
    if skeleton:
        species = list(dict.fromkeys(structure.symbols))
        pseudos = "\n".join(f'{s} = "/absolute/path/to/pseudo/{s}.UPF"' for s in species)
        text += SKELETON.format(pseudos=pseudos)
    return text


# ---------------------------------------------------------------- command line


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog=CONVERTER,
        description="Write a BAANDI [system] block (or a full input skeleton) from a "
        "structure file read by ASE/pymatgen or from a Materials Project id.",
    )
    parser.add_argument("source", nargs="?", help="structure file (any ASE-readable format)")
    parser.add_argument("--mp", metavar="ID", help="Materials Project id, for example mp-149")
    parser.add_argument("-o", "--output", help="output TOML file (default: standard output)")
    parser.add_argument("--format", help="ASE format name (default: detected by ASE)")
    parser.add_argument("--index", type=int, default=-1, help="image index (default -1, the last)")
    parser.add_argument("--reader", choices=("ase", "pymatgen"), default="ase")
    parser.add_argument("--skeleton", action="store_true", help="add [pseudopotentials]/[basis]/[xc]/[scf] placeholders")
    parser.add_argument("--units", choices=("angstrom", "bohr"), default="angstrom")
    parser.add_argument(
        "--compact",
        action="store_true",
        help="write the compact [structure] block (one atom per line, fixed/magmom flags) instead of [system]",
    )
    parser.add_argument("--name", help="system.name (default: formula)")
    parser.add_argument("--assume-periodic", action="store_true", help="accept a non-periodic file with a full cell as a periodic box")
    parser.add_argument("--refresh", action="store_true", help="refetch a cached Materials Project record")
    parser.add_argument("--cache-dir", help="cache root (default: BAANDI_CACHE_DIR or ~/.cache/baandi)")
    parser.add_argument("--force", action="store_true", help="overwrite an existing output file")
    args = parser.parse_args(argv)
    if (args.source is None) == (args.mp is None):
        parser.error("give exactly one of a structure file or --mp ID")
    try:
        if args.mp:
            record, path = fetch_mp(args.mp, args.cache_dir, args.refresh)
            structure = structure_from_mp_record(record, path)
        else:
            structure = read_file(
                args.source, args.format, args.index, args.reader, args.assume_periodic
            )
        text = to_toml(structure, args.units, args.skeleton, args.name, args.compact)
        if args.output:
            out = Path(args.output)
            if out.exists() and not args.force:
                raise ConversionError(f"{out} exists; pass --force to overwrite")
            out.write_text(text)
        else:
            sys.stdout.write(text)
    except ConversionError as e:
        print(f"{CONVERTER}: error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
