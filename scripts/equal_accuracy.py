#!/usr/bin/env python3
"""Work versus delivered accuracy, never cross-code validation or certification.

Stdlib only. Each code uses its own tight reference. Missing/failed references
produce null accuracy, never zero. Exit 3 means incomplete scientific evidence.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from compare_materials_qe import parse_output
from qe_common import Error
from qe_rhog import FormatError, load_ours, cross, dot, read_qe
from verify_qe import require_strict_result

CASES = ('si8', 'k444s', 'si-20', 'mgo-40', 'syn20-skew',
         'syn4-elongated', 'syn4-slab', 'syn16-gamma', 'syn16-k211')
THRESHOLDS = (5e-8, 5e-9, 5e-10)  # experiment candidates, NOT frozen defaults


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def finite(x):
    require(type(x) in (int, float) and math.isfinite(x), 'nonfinite/non-numeric measurement')
    return x


def module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def scf_settings(text, settings):
    """Edit only the generated input's SCF table, preserving all physics."""
    lines = text.splitlines()
    start = lines.index('[scf]') + 1
    end = next((i for i in range(start, len(lines)) if lines[i].lstrip().startswith('[')), len(lines))
    body = [line for line in lines[start:end] if line.split('=', 1)[0].strip() not in settings]
    body.extend(f'{key} = {json.dumps(value)}' for key, value in settings.items())
    return '\n'.join(lines[:start] + body + lines[end:]) + '\n'


def case_input(name, temp, upf, threads):
    require(name in CASES, 'unknown case: ' + name)
    if name not in ('si-20', 'mgo-40'):
        replay = module(ROOT / 'tests/p15a_replay.py')
        text, pseudo = replay.inputs(name, temp, upf)
        hashes = {pseudo.name: sha(pseudo)}
    else:
        nc = Path(os.environ['DFT_NC_DIR']).resolve()
        source = ROOT / 'tests/reference/materials' / (name + '.toml')
        text = source.read_text()
        manifest = ROOT / 'docs/data/nc-datasets.json'
        require(sha(manifest) in text, 'material manifest hash mismatch')
        datasets = {d['element']: d for d in json.loads(manifest.read_text())['datasets']}
        hashes = {}
        for element in sorted(set(re.findall(r'species = "([A-Za-z]+)"', text))):
            pseudo = nc / (element + '.upf')
            digest = sha(pseudo)
            require(digest == datasets[element]['sha256'] and digest in text, 'material UPF hash mismatch')
            hashes[pseudo.name] = digest
            text = text.replace(json.dumps(pseudo.name), json.dumps(str(pseudo)))
    return scf_settings(text, {'threads': threads}), hashes


def profile(result):
    labels = [s['convergence_profile'] for s in
              (result, result.get('scf', {}), result.get('convergence', {}))
              if 'convergence_profile' in s]
    require(len(set(labels)) <= 1, 'conflicting convergence profile labels')
    value = labels[0] if labels else 'strict'
    require(value in ('strict', 'standard'), 'unknown convergence profile')
    return value


def baandi_measurement(result, expected_profile):
    require(profile(result) == expected_profile, 'CLI did not return the requested profile')
    require(result['convergence']['converged'] is True, 'SCF did not converge')
    history = result['convergence']['history']
    require(bool(history), 'missing convergence history')
    certificates = [h['fixed_point_certificate'] for h in history if h.get('fixed_point_certificate')]
    if expected_profile == 'strict':
        require_strict_result(result)
        require(history[-1].get('fixed_point_certificate') is not None, 'strict result lacks final certificate')
    else:
        require(not certificates, 'standard result performed a certificate solve')
        require(any(s.get('certified') is False for s in
                    (result, result.get('convergence', {}))), 'standard result must explicitly be uncertified')
        require(not any(s.get('certified') is True for s in
                        (result, result.get('scf', {}), result.get('convergence', {}))),
                'standard result claims certification')
    points = result.get('kpoints', {}).get('points', [dict(fractional=[0, 0, 0], weight=1, eigenvalues=result['eigenvalues'])])
    occupied = sum(finite(x) > 0 for x in result['occupations'])
    require(0 < occupied < result['bands'], 'occupied plus first empty roots required')
    spectra = [[finite(x) for x in p['eigenvalues']['ha'][:occupied+1]] for p in points]
    require(all(len(s) == occupied+1 for s in spectra), 'incomplete eigenvalues')
    work = history + certificates
    guard = result.get('spectral_completeness') or {}
    return dict(total_energy_ha=finite(result['total_energy']['ha']), eigenvalues_ha=spectra,
                occupied=occupied, kpoints=[dict(fractional=p['fractional'], weight=p['weight']) for p in points],
                scf_iterations=len(history), davidson_steps=sum(h['davidson_iterations'] for h in work),
                h_vectors=sum(h['operator_applications'] for h in work) + guard.get('extra_operator_applications', 0),
                certified=expected_profile == 'strict', convergence_profile=expected_profile)


def accuracy(row, reference):
    require(row['occupied'] == reference['occupied'], 'occupation count mismatch')
    require(row['kpoints'] == reference['kpoints'], 'k point order/weights mismatch')
    require(len(row['eigenvalues_ha']) == len(reference['eigenvalues_ha']), 'spectrum count mismatch')
    differences = []
    empty = []
    for a, b in zip(row['eigenvalues_ha'], reference['eigenvalues_ha']):
        require(len(a) == len(b) == row['occupied']+1, 'root count mismatch')
        differences.extend(abs(x-y) for x, y in zip(a[:-1], b[:-1]))
        empty.append(abs(a[-1]-b[-1]))
    return dict(energy_difference_ha=row['total_energy_ha']-reference['total_energy_ha'],
                max_occupied_eigenvalue_difference_ha=max(differences),
                max_first_empty_eigenvalue_difference_ha=max(empty),
                max_eigenvalue_difference_ha=max(differences+empty),
                density_l2_bohr_m3_over_2=None, density_hartree_energy_ha=None,
                density_hartree_norm_sqrt_ha=None)


def density_data(directory):
    coeff, meta, _ = load_ours(directory)
    raw = (Path(directory) / 'density.bin').read_bytes()
    return coeff, meta, [v[0] for v in struct.iter_unpack('<d', raw)]


def density_accuracy(a, b):
    ca, ma, ra = a
    cb, mb, rb = b
    for key in ('fft_grid', 'lattice_columns_bohr', 'volume_bohr3'):
        require(ma[key] == mb[key], 'density geometry mismatch: ' + key)
    require(ca.keys() == cb.keys(), 'density sphere mismatch')
    require(len(ra) == len(rb) == math.prod(ma['fft_grid']), 'density sample count mismatch')
    volume = ma['volume_bohr3']
    columns = list(zip(*ma['lattice_columns_bohr']))
    determinant = dot(columns[0], cross(columns[1], columns[2]))
    reciprocal = [[2*math.pi*x/determinant for x in cross(columns[j], columns[k])]
                  for j, k in ((1, 2), (2, 0), (0, 1))]
    terms = []
    for m in ca:
        if m == (0, 0, 0):
            continue
        g = [sum(m[j]*reciprocal[j][i] for j in range(3)) for i in range(3)]
        terms.append(abs(ca[m]-cb[m])**2 / dot(g, g))
    energy = 2*math.pi*volume*math.fsum(terms)
    return dict(density_l2_bohr_m3_over_2=math.sqrt(volume/len(ra)*math.fsum((x-y)**2 for x, y in zip(ra, rb))),
                density_hartree_energy_ha=energy, density_hartree_norm_sqrt_ha=math.sqrt(energy))


def run_one(binary, text, directory, label, expected_profile, timeout):
    directory.mkdir()
    inp = directory / 'input.toml'
    inp.write_text(text)
    output = directory / 'output'
    command = [str(binary), 'run', str(inp), '--output-dir', str(output)]
    row = dict(code='BAANDI', label=label, command=command, input=text, input_sha256=sha(inp),
               convergence_profile=expected_profile, status='failed', accuracy=None,
               binary_sha256=sha(binary) if binary.exists() else None)
    start = time.monotonic()
    try:
        proc = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        row.update(exit_code=proc.returncode, wall_s=time.monotonic()-start,
                   wall_scope='process, including startup',
                   stderr=(proc.stderr if len(proc.stderr) <= 4000 else
                           proc.stderr[:2000]+'\n[excerpt truncated]\n'+proc.stderr[-2000:]))
        (directory / 'stdout.txt').write_text(proc.stdout)
        (directory / 'stderr.txt').write_text(proc.stderr)
        if (output / 'results.json').exists():
            result = json.loads((output / 'results.json').read_text())
            history = result.get('convergence', {}).get('history', [])
            work = history + [h['fixed_point_certificate'] for h in history if h.get('fixed_point_certificate')]
            row.update(scf_iterations=len(history),
                       davidson_steps=sum(h['davidson_iterations'] for h in work),
                       h_vectors=sum(h['operator_applications'] for h in work),
                       converged=result.get('convergence', {}).get('converged', False))
        require(proc.returncode == 0, f'CLI exit {proc.returncode}')
        result = json.loads((output / 'results.json').read_text())
        row.update(baandi_measurement(result, expected_profile), status='ok',
                   result_sha256=sha(output / 'results.json'))
        density = density_data(output)
        return row, density
    except (ValueError, KeyError, OSError, Error, FormatError, subprocess.TimeoutExpired) as exc:
        row.update(status='failed', error=str(exc), wall_s=time.monotonic()-start)
        return row, None


def baandi_case(args, name, directory):
    with tempfile.TemporaryDirectory(prefix='p17a-pseudo-') as temp:
        text, hashes = case_input(name, Path(temp), args.upf.resolve(), args.threads)
        variants = [('tight', dict(energy_tolerance_ha=1e-12, density_residual_tolerance=1e-9, eigen_tolerance=1e-10)),
                    ('tighter', dict(energy_tolerance_ha=1e-13, density_residual_tolerance=1e-10, eigen_tolerance=1e-11)),
                    ('strict', {})]
        if not args.strict_only:
            variants += [(f'standard-{t:g}', dict(convergence_profile='standard', standard_threshold_ha=t))
                         for t in args.thresholds]
            variants += [('strict-hartree', dict(mixing_metric='hartree'))]
        rows, densities = [], []
        for label, settings in variants:
            row, density = run_one(args.binary.resolve(), scf_settings(text, settings),
                                   directory / label, label, settings.get('convergence_profile', 'strict'), args.timeout)
            row.update(case=name, upf_sha256=hashes)
            rows.append(row)
            densities.append(density)
        tight, tighter = rows[:2]
        verified = tight['status'] == tighter['status'] == 'ok' and abs(tight['total_energy_ha']-tighter['total_energy_ha']) <= 1e-12
        reference_check = dict(verified=verified, criterion='two strict runs; |E_tight-E_tighter| <= 1e-12 Ha',
                               energy_difference_ha=(tight['total_energy_ha']-tighter['total_energy_ha'])
                               if tight['status'] == tighter['status'] == 'ok' else None)
        for row, density in zip(rows, densities):
            row['reference_check'] = reference_check
            if verified and row['status'] == 'ok':
                row['accuracy'] = accuracy(row, tighter)
                row['accuracy'].update(density_accuracy(density, densities[1]))
                row['reference_label'] = 'tighter'
        return rows


def qe_row(path):
    text = path.read_text()
    q = parse_output(text)
    occupied = int(q['electrons']/2)
    require(occupied*2 == q['electrons'] and 0 < occupied < q['bands'], 'QE integer occupations and first empty required')
    # Exact timer counts identify initial solves plus Davidson expansions;
    # rounded avg-iteration records do not. h_psi calls remain block counts.
    timers = {name: int(count) for name, count in re.findall(
        r'^\s*([a-z_0-9:]+)\s*:\s*[^\n]*WALL\s*\(\s*(\d+) calls\)', text, re.M)}
    solvers = [name for name in ('cegterg', 'regterg') if name in timers and name+':over' in timers]
    steps = sum(timers[name]+timers[name+':over'] for name in solvers) if solvers else None
    wall = re.search(r'^\s*PWSCF\s*:\s*.+?CPU\s+(.+?)\s+WALL\s*$', text, re.M)
    seconds = None
    if wall:
        tokens = re.findall(r'([0-9.]+)\s*([dhms])', wall[1])
        if tokens:
            seconds = sum(float(v)*dict(d=86400, h=3600, m=60, s=1)[u] for v, u in tokens)
    # Explicit-output comparisons also bind printed geometry and PP identity.
    # A manifest additionally binds SHA-256 and the complete physical inputs.
    scalar = r'([-+0-9.eEdD]+)'
    alat = re.search(r'lattice parameter \(alat\)\s*=\s*'+scalar, text)
    axes = re.findall(r'a\(\s*[123]\s*\)\s*=\s*\(\s*'+scalar+r'\s+'+scalar+r'\s+'+scalar+r'\s*\)', text)
    atoms = re.findall(r'\b([A-Z][a-z]?)\s+tau\(\s*\d+\s*\)\s*=\s*\(\s*'+scalar+r'\s+'+scalar+r'\s+'+scalar+r'\s*\)', text)
    pseudo_md5 = re.findall(r'MD5 check sum:\s*([0-9a-f]+)', text)
    geometry = dict(alat=alat[1], axes=axes, atoms=atoms, pseudo_md5=pseudo_md5) if alat and len(axes)==3 and atoms and pseudo_md5 else None
    energy_token = re.search(r'!\s+total energy\s*=\s*([-+\d.EeDd]+)', text)[1]
    # Half-last-place precision in Ry, then converted to Ha.
    token = energy_token.lower().replace('d', 'e')
    mantissa, _, exponent = token.partition('e')
    decimal_places = len(mantissa.partition('.')[2])
    rounding_ha = 0.25 * 10**(int(exponent or '0')-decimal_places)
    requested = re.search(r'scf convergence threshold\s*=\s*([-+\d.eEdD]+)', text)
    estimates = [float(v.replace('D','e').replace('d','e'))/2 for v in re.findall(r'estimated scf accuracy\s*[<=>]+\s*([-+\d.eEdD]+)\s+Ry', text)]
    return dict(requested_scf_threshold_ry=float(requested[1].lower().replace('d','e')) if requested else None, estimated_scf_errors_ha=estimates, final_estimated_scf_error_ha=estimates[-1] if estimates else None, code='QE', convergence_profile='qe-conv-thr', certified=False, status='ok',
                output=str(path), output_sha256=sha(path), version=q['version'],
                total_energy_ha=q['energy_ry']/2, energy_rounding_bound_ha=rounding_ha,
                eigenvalues_ha=[v[:occupied+1] for v in q['eigenvalues_ha']], occupied=occupied,
                kpoints=[dict(fractional=k[:3], weight=k[3]/2) for k in q['crystal']],
                scf_iterations=q['iterations'], davidson_steps=steps, h_vectors=None,
                davidson_steps_definition='initial solve calls plus expansion timer calls; unavailable without timers',
                h_blocks=timers.get('h_psi'),
                wall_s=seconds, wall_scope='PWSCF internal WALL; rounded; excludes launcher',
                model=dict(printed_geometry_and_pp=geometry, grid=q['grid'], ecutwfc=q['ecutwfc'], ecutrho=q['ecutrho'], xc=q['xc'],
                           electrons=q['electrons'], bands=q['bands']),
                accuracy=None)


def qe_density_accuracy(path, reference, output_text):
    """Same-code reciprocal densities; G vectors from stdout alat/b(1..3)."""
    a, ma = read_qe(path)
    b, mb = read_qe(reference)
    require(a.keys() == b.keys(), 'QE density Miller sets differ')
    require(ma['b_vectors_as_stored'] == mb['b_vectors_as_stored'], 'QE density reciprocal headers differ')
    n = r'([-+0-9.eEdD]+)'
    number = lambda v: float(v.lower().replace('d', 'e'))
    alat_match = re.search(r'lattice parameter \(alat\)\s*=\s*'+n, output_text)
    require(alat_match is not None, 'QE density requires stdout alat')
    alat = number(alat_match[1])
    rows = re.findall(r'b\(\s*([123])\s*\)\s*=\s*\(\s*'+n+r'\s+'+n+r'\s+'+n+r'\s*\)', output_text)
    require([int(r[0]) for r in rows] == [1,2,3], 'QE density requires unique reciprocal basis in stdout')
    reciprocal = [[2*math.pi*number(x)/alat for x in row[1:]] for row in rows]
    volume = (2*math.pi)**3 / abs(dot(reciprocal[0], cross(reciprocal[1], reciprocal[2])))
    energy = 0.0
    for m in a:
        if m != (0,0,0):
            g = [sum(m[j]*reciprocal[j][i] for j in range(3)) for i in range(3)]
            energy += abs(a[m]-b[m])**2 / dot(g,g)
    energy *= 2*math.pi*volume
    return dict(density_l2_bohr_m3_over_2=math.sqrt(volume*math.fsum(abs(a[m]-b[m])**2 for m in a)),
                density_hartree_energy_ha=energy, density_hartree_norm_sqrt_ha=math.sqrt(energy))


def qe_pairs(manifest):
    """Manifest pairs bind raw inputs, outputs, commands and UPF hashes."""
    entries = json.loads(manifest.read_text())
    rows = []
    for entry in entries:
        require(bool(entry['upf_sha256']) and all(re.fullmatch('[0-9a-f]{64}', v)
                for v in entry['upf_sha256'].values()), 'QE UPF SHA-256 provenance required')
        pair = []
        models = []
        for key, threshold in (('standard', 1e-8), ('tight', 1e-14)):
            spec = entry[key]
            require(bool(spec['command']), 'QE run command required')
            inp = (manifest.parent / spec['input']).resolve()
            text = inp.read_text()
            match = re.findall(r'\bconv_thr\s*=\s*([-+\d.eEdD]+)', text, re.I)
            require(len(match) == 1 and float(match[0].lower().replace('d', 'e')) == threshold,
                    'QE input conv_thr mismatch')
            # Only convergence threshold and run-directory naming may differ.
            normalized = re.sub(r'\bconv_thr\s*=\s*[-+\d.eEdD]+', 'conv_thr=THRESHOLD', text, flags=re.I)
            normalized = re.sub(r"\b(prefix|outdir)\s*=\s*['\"][^'\"]*['\"]", r'\1=PATH', normalized, flags=re.I)
            models.append(normalized)
            row = qe_row((manifest.parent / spec['output']).resolve())
            row.update(case=entry['case'], label=key, conv_thr_ry=threshold, command=spec['command'],
                       input_sha256=sha(inp), upf_sha256=entry['upf_sha256'])
            pair.append(row)
        require(models[0] == models[1], 'QE inputs differ beyond conv_thr/prefix/outdir')
        require(pair[0]['model'] == pair[1]['model'] and pair[0]['version'] == pair[1]['version'], 'QE model/version mismatch')
        for row in pair:
            row['accuracy'] = accuracy(row, pair[1])
            row['reference_label'] = 'tight'
            row['energy_difference_rounding_bound_ha'] = row['energy_rounding_bound_ha'] + pair[1]['energy_rounding_bound_ha']
            row['reference_check'] = dict(verified=False, criterion='QE conv_thr=1e-14; stdout alone does not prove 1e-12 Ha convergence')
        rows.extend(pair)
    return rows


def markdown(rows):
    lines = ['# Work versus delivered accuracy', '',
             'Errors are against each code\'s own tight reference. Missing values are not zero.', '',
             '| Case | Code / profile | Status | SCF | Davidson | H vectors | Wall (s) | dE (Ha) | max dEigen (Ha) | density L2 | sqrt(E_H) |',
             '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    def show(x):
        return 'unavailable' if x is None else f'{x:.9g}' if isinstance(x, float) else str(x)
    for row in rows:
        a = row.get('accuracy') or {}
        values = [row['case'], row['code']+'/'+row['label'], row['status'], row.get('scf_iterations'),
                  row.get('davidson_steps'), row.get('h_vectors'), row.get('wall_s'),
                  a.get('energy_difference_ha'), a.get('max_eigenvalue_difference_ha'),
                  a.get('density_l2_bohr_m3_over_2'), a.get('density_hartree_norm_sqrt_ha')]
        lines.append('| ' + ' | '.join(show(v) for v in values) + ' |')
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT / 'target/release/dft-cli')
    parser.add_argument('--upf', type=Path, default=ROOT.parent / 'pseudo/Si.pz-vbc.UPF')
    parser.add_argument('--cases', default=','.join(CASES), help='comma-separated cases; empty for QE-only')
    parser.add_argument('--thresholds', type=float, nargs=3, default=THRESHOLDS)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--timeout', type=float, default=1800)
    parser.add_argument('--strict-only', action='store_true', help='explicit preparatory run; no standard/E4 acceptance')
    parser.add_argument('--qe-output', type=Path, action='append', default=[], help='explicit pw.x output; repeatable')
    parser.add_argument('--qe-reference', type=Path, help='same-model tight pw.x output for --qe-output')
    parser.add_argument('--qe-density', type=Path, help='optional charge-density.dat for a single --qe-output')
    parser.add_argument('--qe-reference-density', type=Path)
    parser.add_argument('--qe-pairs', type=Path, help='JSON manifest of raw QE output pairs')
    parser.add_argument('--output-dir', type=Path, required=True, help='fresh directory, never overwritten')
    args = parser.parse_args()
    require(all(name in CASES for name in filter(None, args.cases.split(','))), 'unknown case')
    require(args.threads > 0 and 0 < args.timeout < math.inf, 'positive threads and finite timeout required')
    require(all(0 < t < math.inf for t in args.thresholds) and len(set(args.thresholds)) == 3, 'three distinct positive finite thresholds required')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    rows = []
    for name in filter(None, args.cases.split(',')):
        directory = args.output_dir / name
        directory.mkdir()
        try:
            rows.extend(baandi_case(args, name, directory))
        except (ValueError, KeyError, OSError, Error, FormatError) as exc:
            rows.append(dict(case=name, code='BAANDI', label='setup', status='failed', error=str(exc), accuracy=None))
        print(name + ': ' + rows[-1]['status'], flush=True)
    if args.qe_pairs:
        try:
            rows.extend(qe_pairs(args.qe_pairs.resolve()))
        except (ValueError, KeyError, OSError, Error, FormatError) as exc:
            rows.append(dict(case='QE-pairs', code='QE', label='setup', status='failed', error=str(exc), accuracy=None))
    if args.qe_density or args.qe_reference_density:
        require(args.qe_density and args.qe_reference_density and args.qe_reference and len(args.qe_output) == 1,
                'QE densities require a single output, reference output and both density paths')
    for path in args.qe_output:
        try:
            row = qe_row(path)
            row.update(case=str(path.parent), label='stored-output', reference_available=False)
            if args.qe_reference:
                reference = qe_row(args.qe_reference)
                require(row['model']['printed_geometry_and_pp'] is not None, 'explicit QE comparison requires printed geometry and PP checksums; use a fully bound manifest otherwise')
                require(row['model'] == reference['model'] and row['version'] == reference['version'], 'QE model/version mismatch')
                row['accuracy'] = accuracy(row, reference)
                row['reference_available'] = True
                row['reference_output_sha256'] = reference['output_sha256']
                row['energy_difference_rounding_bound_ha'] = row['energy_rounding_bound_ha'] + reference['energy_rounding_bound_ha']
                row['reference_check'] = dict(verified=False, criterion='explicit user-supplied output; stdout is not a 1e-12 Ha proof')
                if args.qe_density:
                    row['accuracy'].update(qe_density_accuracy(args.qe_density, args.qe_reference_density, path.read_text()))
            rows.append(row)
        except (ValueError, KeyError, OSError, Error, FormatError) as exc:
            rows.append(dict(case=str(path), code='QE', label='stored-output', status='failed', error=str(exc), accuracy=None))
    complete = bool(rows) and all(r['status'] == 'ok' and r.get('accuracy') is not None for r in rows)
    report = dict(schema_version=1, command=[sys.executable, *sys.argv], rows=rows,
                  binary_sha256=sha(args.binary) if args.cases else None,
                  experiment_thresholds_ha=args.thresholds, frozen_standard_default_ha=5e-9,
                  complete=complete, strict_only=args.strict_only,
                  limitations=['A two-level energy plateau is empirical, not proof of the ground state or spectral completeness.',
                    'Hartree norm is sqrt(2*pi*Omega*sum_G!=0 |delta n_G|^2/G^2), density sphere only; no electron-count division.',
                    'L2 is sqrt(integral(delta n^2)), full real-space grid, in Bohr^-3/2.',
                    'QE Davidson steps use initial+expansion timer calls when available; H-vector counts remain unknown, h_psi calls are blocks. Optional density norms use rounded stdout geometry and the stored G set.',
                    'QE rounded equal energies do not establish zero delivered error; rounding bounds are reported.',
                    'Wall samples are local observations; no speedup acceptance or default calibration is implied.'])
    (args.output_dir / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    (args.output_dir / 'table.md').write_text(markdown(rows))
    return 0 if complete else 3


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
