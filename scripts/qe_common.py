"""Shared, standard-library-only QE interchange contracts (no QE source used)."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

EV_PER_HA = 27.211386245988
UNVERIFIED_SHA = 'e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217'
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?'


class Error(ValueError):
    def __init__(self, message, code=1):
        super().__init__(message)
        self.code = code


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Error(message)


def number(value):
    value = float(str(value).replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(value):
        raise Error('expected finite number')
    return value


def positive(value):
    value = number(value)
    if value <= 0:
        raise Error('expected positive number')
    return value


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise Error(f'duplicate JSON key: {key}')
            result[key] = value
        return result
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(Error(f'nonfinite JSON: {x}')))


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def results(path):
    path = Path(path)
    if path.is_dir():
        path /= 'results.json'
    data = read_json(path)
    try:
        if data['schema_version'] != 1 or data['xc'] not in ('lda-pz', 'pbe'):
            raise Error('expected results schema_version=1 and xc=lda-pz or pbe')
        if len(data['species']) != 1:
            raise Error('single-species protocol required')
        species = data['species'][0]
        if data['xc'] == 'pbe':
            allow = data['allow_functional_override']
            applied = data['functional_override_applied']
            if type(allow) is not bool or type(applied) is not bool:
                raise Error('PBE functional override flags must be boolean')
            label = species['functional']
            if not isinstance(label, str):
                raise Error('expected UPF functional label')
            mismatch = ' '.join(label.upper().split()) not in ('PBE', 'SLA PW PBX PBC')
            if mismatch and not allow:
                raise Error('PBE UPF mismatch requires explicit allow_functional_override=true')
            if applied != mismatch:
                raise Error('functional_override_applied disagrees with UPF label')
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', species['name']):
            raise Error('invalid species label')
        if not re.fullmatch('[0-9a-f]{64}', species['sha256']):
            raise Error('expected species sha256')
        a = data['cell']['lattice_columns_bohr']
        if len(a) != 3 or any(len(row) != 3 for row in a):
            raise Error('expected 3x3 lattice_columns_bohr')
        for row in a:
            for x in row:
                number(x)
        if not data['atoms']:
            raise Error('expected atoms')
        for atom in data['atoms']:
            if atom['species'] != species['name'] or len(atom['fractional']) != 3:
                raise Error('expected single-species fractional positions')
            for x in atom['fractional']:
                number(x)
        for key in ('n_plane_waves',):
            if type(data['basis'][key]) is not int or data['basis'][key] <= 0:
                raise Error(f'expected positive basis.{key}')
        grid = data['basis']['fft_grid']
        if len(grid) != 3 or any(type(x) is not int or x <= 0 for x in grid):
            raise Error('expected three positive fft_grid integers')
        positive(data['cutoff']['ecutwfc_ry'])
        ne = positive(data['electrons'])
        bands = data['bands']
        if ne % 2 or type(bands) is not int or bands < ne / 2:
            raise Error('expected even electrons and sufficient integer bands')
        if data['occupations'] != [2.0] * int(ne / 2) + [0.0] * (bands - int(ne / 2)):
            raise Error('expected fixed spin-unpolarized occupations')
        if not isinstance(data['scf'], dict):
            raise Error('expected scf options')
    except (KeyError, TypeError, IndexError) as exc:
        raise Error(f'malformed results.json: expected field/shape {exc}') from exc
    return path, data


def fmt(x):
    return format(number(x), '.17g')


def quoted(text):
    if not text or any(c in text for c in "'\"\n\r!&"):
        raise Error('QE text must be nonempty and contain no quotes, newlines, ! or &')
    return "'" + text + "'"


def materials_results(path):
    """Full CLI records only; summaries cannot reconstruct a measured FFT grid."""
    path = Path(path)
    if path.is_dir():
        path /= 'results.json'
    r = read_json(path)
    if 'cases' in r:
        raise Error('materials summary lacks structure/grid/spectrum; use per-run results.json')
    try:
        if r['schema_version'] != 1 or r.get('spin', {'mode': 'unpolarized'}) != {'mode': 'unpolarized'}:
            raise Error('unsupported spin: only schema 1 unpolarized results supported')
        if r.get('occupation_policy', 'fixed_insulator') != 'fixed_insulator' or r.get('smearing') is not None:
            raise Error('unsupported smearing/occupation policy')
        if 'state_contracts' in r:
            raise Error('state_contracts spin/occupation extension not supported')
        if r['xc'] not in ('lda-pz', 'pbe'):
            raise Error('unsupported XC functional')
        species = r['species']
        names = [s['name'] for s in species]
        if not names or len(set(names)) != len(names):
            raise Error('missing or duplicate species')
        mismatch = False
        for s in species:
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', s['name']):
                raise Error('invalid species label')
            if not re.fullmatch('[0-9a-f]{64}', s['sha256']):
                raise Error('expected species sha256')
            if r['xc'] == 'pbe':
                mismatch |= ' '.join(s['functional'].upper().split()) not in ('PBE', 'SLA PW PBX PBC')
        if r['xc'] == 'pbe':
            if any(type(r[k]) is not bool for k in ('allow_functional_override', 'functional_override_applied')):
                raise Error('functional override flags must be boolean')
            if (mismatch and not r['allow_functional_override']) or r['functional_override_applied'] != mismatch:
                raise Error('inconsistent functional override provenance')
        a = r['cell']['lattice_columns_bohr']
        if len(a) != 3 or any(len(row) != 3 for row in a):
            raise Error('expected 3x3 lattice')
        for row in a:
            for x in row:
                number(x)
        if not r['atoms']:
            raise Error('expected atoms')
        for atom in r['atoms']:
            if atom['species'] not in names or len(atom['fractional']) != 3:
                raise Error('invalid atom species/position')
            for x in atom['fractional']:
                number(x)
        grid = r['basis']['fft_grid']
        if len(grid) != 3 or any(type(n) is not int or n <= 0 for n in grid):
            raise Error('expected three positive FFT dimensions')
        if 2 * positive(r['cutoff']['ecutwfc_ha']) != positive(r['cutoff']['ecutwfc_ry']):
            raise Error('inconsistent Ha/Ry cutoff')
        ne, nb = positive(r['electrons']), r['bands']
        if ne % 2 or type(nb) is not int or nb < ne / 2:
            raise Error('expected even electrons and sufficient bands')
        fixed = [2.0] * int(ne / 2) + [0.0] * (nb - int(ne / 2))
        if r['occupations'] != fixed or any(p.get('occupations', fixed) != fixed for p in r.get('kpoints', {}).get('points', [])):
            raise Error('unsupported smearing: expected fixed unpolarized occupations')
        if not re.fullmatch('[0-9a-f]{64}', r['input_sha256']):
            raise Error('expected input_sha256')
        if not isinstance(r['tool']['git_commit'], str) or not r['tool']['git_commit']:
            raise Error('expected BAANDI git_commit provenance')
    except (KeyError, TypeError, IndexError, AttributeError) as exc:
        raise Error(f'malformed full results.json: {exc}') from exc
    return path, r


def manifest_species(r, manifest):
    """Resolve identities from metadata only. Never open a pseudopotential."""
    entries = read_json(manifest)['datasets']
    selected = []
    for s in r['species']:
        hits = [d for d in entries if d.get('element') == s['name']]
        if len(hits) != 1:
            raise Error('missing or ambiguous manifest species: ' + s['name'])
        d = hits[0]
        if not re.fullmatch('[0-9a-f]{64}', d.get('sha256', '')) or d['sha256'] != s['sha256']:
            raise Error('missing/mismatched manifest hash: ' + s['name'])
        filename = d['filename']
        if Path(filename).name != filename or filename in ('.', '..') or any(c.isspace() for c in filename):
            raise Error('manifest filename must be a basename')
        quoted(filename)
        if 'path' in s and Path(s['path']).name != filename:
            raise Error('results/manifest UPF filename mismatch')
        if d['pseudo_type'] != 'NC' or d['relativistic'] not in ('scalar', 'no'):
            raise Error('only scalar/nonrelativistic NC datasets supported')
        selected.append(d)
    return selected
