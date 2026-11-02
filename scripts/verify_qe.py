#!/usr/bin/env python3
"""Fail-closed comparison of rust-dft results with a returned QE bundle."""
import datetime
import math
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_common import (Error, Parser, NUMBER, EV_PER_HA, UNVERIFIED_SHA, number,
                       positive, read_json, results, sha256, write_json)

DEFAULTS = dict(energy=1e-4, term=1e-4, eigenvalue=1e-3, gap=1e-3, electrons=1e-8)
REQUIRED_PROVENANCE = ('command', 'qe_version', 'host', 'mpi_processes', 'omp_threads',
                       'date_utc', 'input_sha256', 'upf_sha256', 'wall_time_s',
                       'returned_by', 'notes', 'is_format_fixture')


def match(text, pattern, label, multiple=False):
    found = list(re.finditer(pattern, text, re.MULTILINE))
    if not found or (not multiple and len(found) != 1):
        raise Error(f'expected {label}: /{pattern}/ (found {len(found)})', 3)
    return found[-1]


def parse_output(text):
    """Parse documented literal records, not arbitrary QE variants or multiple runs."""
    n = '(' + NUMBER + ')'
    q = {}
    q['version'] = match(text, r'^\s*(Program PWSCF\s+v\.\S+\s+starts[^\n]*)$', 'Program PWSCF version').group(1).strip()
    for key, label, pattern, cast in [
        ('electrons', 'number of electrons', r'number of electrons\s*=\s*' + n, number),
        ('bands', 'number of Kohn-Sham states', r'number of Kohn-Sham states\s*=\s*(\d+)', int),
        ('ecutwfc', 'kinetic-energy cutoff', r'kinetic-energy cutoff\s*=\s*' + n + r'\s+Ry', number),
        ('ecutrho', 'charge density cutoff', r'charge density cutoff\s*=\s*' + n + r'\s+Ry', number),
        ('nk', 'number of k points', r'number of k points\s*=\s*(\d+)', int),
        ('xc', 'Exchange-correlation', r'Exchange-correlation\s*=\s*([^\n]+)', str),
        ('iterations', 'convergence has been achieved', r'convergence has been achieved in\s+(\d+)\s+iterations', int),
    ]:
        q[key] = cast(match(text, r'^\s*' + pattern + r'\s*$', label).group(1))
    if q['nk'] != 1:
        raise Error('expected number of k points = 1', 3)
    grid = match(text, r'^\s*Dense\s+grid:[^\n]*FFT dimensions:\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*$', 'Dense grid FFT dimensions')
    q['grid'] = [int(x) for x in grid.groups()]
    end = match(text, r'^\s*End of self-consistent calculation\s*$', 'End of self-consistent calculation')
    final = text[end.end():]
    k = match(final, r'^\s*k\s*=\s*' + n + r'\s+' + n + r'\s+' + n +
              r'\s*\(\s*(\d+)\s+PWs\)\s+bands \(ev\):\s*$', 'single final k = ... ( N PWs) bands (ev):')
    if any(number(x) != 0 for x in k.groups()[:3]):
        raise Error('expected final Gamma k = 0 0 0', 3)
    q['npw'] = int(k.group(4))
    values = []
    lines = final[k.end():].splitlines()
    started = False
    for line in lines:
        if not line.strip():
            if started:
                break
            continue
        tokens = line.split()
        if any(not re.fullmatch(NUMBER, x) for x in tokens):
            if not started or len(values) != q['bands']:
                raise Error('expected eigenvalues: exactly nbnd numeric values in eV', 3)
            break
        started = True
        values.extend(number(x) for x in tokens)
    if len(values) != q['bands']:
        raise Error('expected eigenvalues: exactly nbnd numeric values in eV', 3)
    if any(x > y for x, y in zip(values, values[1:])):
        raise Error('expected sorted eigenvalues in eV', 3)
    q['eigenvalues_ha'] = [x / EV_PER_HA for x in values]
    energy = match(text, r'^!\s+total energy\s*=\s*' + n + r'\s+Ry\s*$', '! total energy')
    q['total'] = number(energy.group(1)) / 2
    decomposition = text[energy.end():]
    for key, label in [('one_electron', 'one-electron'), ('hartree', 'hartree'), ('xc', 'xc'), ('ewald', 'ewald')]:
        m = match(decomposition, r'^\s*' + label + r' contribution\s*=\s*' + n + r'\s+Ry\s*$', label + ' contribution')
        q[key + '_ha'] = number(m.group(1)) / 2
    q['levels_ev'] = None
    if 'highest occupied' in final:
        pair = 'highest occupied, lowest unoccupied' in final
        label = 'highest occupied, lowest unoccupied levels' if pair else 'highest occupied level'
        pattern = r'^\s*' + re.escape(label) + r' \(ev\):\s*' + n
        if pair:
            pattern += r'\s+' + n
        m = match(final, pattern + r'\s*$', label + ' (ev)')
        q['levels_ev'] = [number(x) for x in m.groups()]
    q['wall_time_text'] = None
    q['wall_time_s'] = None
    if re.search(r'^\s*PWSCF\s*:', text, re.M):
        m = match(text, r'^\s*PWSCF\s*:\s*(.+?)\s+CPU\s+(.+?)\s+WALL\s*$', 'PWSCF CPU/WALL timing')
        duration = r'(?:\d+d\s*)?(?:\d+h\s*)?(?:\d+m\s*)?(?:\d+(?:\.\d+)?s)?'
        if not m.group(2).strip() or not re.fullmatch(duration, m.group(2).strip()):
            raise Error('expected PWSCF wall duration Nd Nh Nm Ns', 3)
        q['wall_time_text'] = m.group(2)
        factors = {'d':86400, 'h':3600, 'm':60, 's':1}
        q['wall_time_s'] = sum(float(value)*factors[unit] for value, unit in
                               re.findall(r'(\d+(?:\.\d+)?)([dhms])', m.group(2)))
    match(text, r'^\s*JOB DONE\.\s*$', 'JOB DONE.')
    return q


def parse_input(text):
    """Strict generated-input subset. Unknown/duplicate keys and cards fail closed."""
    lines = []
    # Comments only outside quoted strings; quoted paths cannot contain !.
    for line in text.splitlines():
        line = line.split('!', 1)[0].strip()
        if line:
            lines.append(line)
    sections = {}
    i = 0
    allowed = {
        'CONTROL': {'calculation', 'prefix', 'pseudo_dir', 'outdir', 'disk_io', 'verbosity', 'tprnfor', 'tstress'},
        'SYSTEM': {'ibrav', 'nat', 'ntyp', 'ecutwfc', 'ecutrho', 'nbnd', 'occupations', 'nosym', 'noinv', 'nr1', 'nr2', 'nr3'},
        'ELECTRONS': {'diagonalization', 'diago_full_acc', 'conv_thr', 'mixing_beta', 'electron_maxstep'},
    }
    for section, keys in allowed.items():
        if i >= len(lines) or lines[i] != '&' + section:
            raise Error('expected &' + section, 3)
        i += 1
        data = {}
        while i < len(lines) and lines[i] != '/':
            m = re.fullmatch(r"([a-z][a-z0-9_]*)\s*=\s*('[^']*'|\.true\.|\.false\.|" + NUMBER + r'),?', lines[i])
            if not m or m.group(1) not in keys or m.group(1) in data:
                raise Error(f'expected unique supported {section} assignment: {lines[i]}', 3)
            data[m.group(1)] = m.group(2)
            i += 1
        if i >= len(lines):
            raise Error('expected namelist terminator /', 3)
        missing = keys - data.keys() - {'nr1', 'nr2', 'nr3'}
        if missing:
            raise Error(f'expected {section} keys: {sorted(missing)}', 3)
        sections[section] = data
        i += 1
    def card(name, count):
        nonlocal i
        if i >= len(lines) or lines[i] != name or i + count >= len(lines):
            raise Error(f'expected {name} and {count} rows', 3)
        rows = [line.split() for line in lines[i+1:i+1+count]]
        i += count + 1
        return rows
    s = sections['SYSTEM']
    for key in ('nat', 'ntyp', 'nbnd', 'ibrav'):
        if not re.fullmatch(r'\d+', s[key]):
            raise Error(f'expected integer {key}', 3)
    species = card('ATOMIC_SPECIES', int(s['ntyp']))
    lattice = card('CELL_PARAMETERS bohr', 3)
    positions = card('ATOMIC_POSITIONS crystal', int(s['nat']))
    kpoints = card('K_POINTS automatic', 1)
    if i != len(lines):
        raise Error('expected end of input after K_POINTS', 3)
    if any(len(row) != 3 for row in lattice) or any(len(row) != 4 for row in positions) or any(len(row) != 3 for row in species):
        raise Error('expected complete species/cell/position rows', 3)
    return sections, species, [[number(x) for x in row] for row in lattice], positions, kpoints


def provenance(bundle, ours):
    try:
        p = read_json(bundle / 'provenance.json')
        # Fixture marker takes precedence over missing fields and all numerical work.
        if p.get('is_format_fixture') is True:
            raise Error('FORMAT FIXTURE: not a reference', 4)
        missing = set(REQUIRED_PROVENANCE) - p.keys()
        if missing:
            raise Error(f'provenance incomplete: missing {sorted(missing)}', 2)
        if type(p['is_format_fixture']) is not bool:
            raise Error('provenance is_format_fixture must be Boolean', 2)
        for key in ('command', 'qe_version', 'host', 'date_utc', 'returned_by'):
            if not isinstance(p[key], str) or not p[key].strip():
                raise Error(f'provenance expected nonempty {key}', 2)
        if not isinstance(p['notes'], str):
            raise Error('provenance expected notes string', 2)
        for key in ('mpi_processes', 'omp_threads'):
            if type(p[key]) is not int or p[key] <= 0:
                raise Error(f'provenance expected positive integer {key}', 2)
        if p['omp_threads'] != 1:
            raise Error('provenance protocol requires omp_threads=1', 2)
        if type(p['wall_time_s']) not in (int, float) or number(p['wall_time_s']) < 0:
            raise Error('provenance expected nonnegative wall_time_s', 2)
        date = datetime.datetime.fromisoformat(p['date_utc'].replace('Z', '+00:00'))
        if date.utcoffset() != datetime.timedelta(0):
            raise Error('provenance expected UTC date_utc', 2)
        for key in ('input_sha256', 'upf_sha256'):
            if not isinstance(p[key], str) or not re.fullmatch('[0-9a-f]{64}', p[key]):
                raise Error(f'provenance expected SHA-256 {key}', 2)
        if sha256(bundle / 'qe.in') != p['input_sha256']:
            raise Error('provenance mismatch: input_sha256', 2)
        if p['upf_sha256'] != ours['species'][0]['sha256']:
            raise Error('provenance mismatch: upf_sha256', 2)
        return p
    except Error as exc:
        if exc.code == 1:
            raise Error(f'provenance invalid: {exc}', 2) from exc
        raise
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        raise Error(f'provenance incomplete or invalid: {exc}', 2) from exc


def conditions(r, q, input_text):
    sections, species, vectors, positions, kpoints = parse_input(input_text)
    failures = []
    def check(ok, description):
        if not ok:
            failures.append(description)
    c, s, e = (sections[x] for x in ('CONTROL', 'SYSTEM', 'ELECTRONS'))
    expected = [('calculation', "'scf'", c), ('verbosity', "'high'", c),
                ('ibrav', '0', s), ('ntyp', '1', s), ('occupations', "'fixed'", s),
                ('nosym', '.true.', s), ('noinv', '.true.', s),
                ('diagonalization', "'david'", e), ('diago_full_acc', '.true.', e)]
    for key, value, section in expected:
        check(section[key] == value, f'input {key} must equal {value}')
    check(kpoints == [['1', '1', '1', '0', '0', '0']], 'complex Gamma K_POINTS automatic 1 1 1 0 0 0')
    for key, target in [('ecutwfc', r['cutoff']['ecutwfc_ry']), ('ecutrho', 4*r['cutoff']['ecutwfc_ry'])]:
        check(number(s[key]) == target, f'input {key}: expected {target} Ry')
        check(q[key] == target, f'output {key}: expected {target} Ry, found {q[key]}')
    check(q['electrons'] == r['electrons'], 'number of electrons differs')
    check(q['bands'] == r['bands'] == int(s['nbnd']), 'number of bands differs')
    check(q['npw'] == r['basis']['n_plane_waves'], 'plane-wave count differs')
    grid_keys = [key for key in ('nr1', 'nr2', 'nr3') if key in s]
    check(len(grid_keys) in (0, 3), 'all or none of nr1/nr2/nr3 required')
    if grid_keys:
        requested = [number(s[key]) for key in grid_keys]
        check(requested == r['basis']['fft_grid'], 'input FFT dimensions differ from ours')
        check(q['grid'] == r['basis']['fft_grid'], 'output FFT dimensions differ from ours')
    check(bool(re.fullmatch(r'(?:SLA\s+PZ\s+NOGX\s+NOGC|SLA\s+PZ)\s*(?:\(\s*1\s+1\s+0\s+0\s*\))?', q['xc'].strip())),
          'functional must be LDA-PZ; QE printed: ' + q['xc'])
    check(len(species) == 1 and species[0][0] == r['species'][0]['name'], 'species label/count differs')
    check(len(positions) == len(r['atoms']), 'atom count differs')
    a = r['cell']['lattice_columns_bohr']
    check(all(abs(vectors[j][i]-a[i][j]) <= 1e-12 for i in range(3) for j in range(3)), 'lattice differs by >1e-12 Bohr')
    for pos, atom in zip(positions, r['atoms']):
        check(pos[0] == atom['species'], 'atomic species differs')
        delta = [sum(vectors[j][i]*number(pos[j+1]) - a[i][j]*atom['fractional'][j]
                     for j in range(3)) for i in range(3)]
        check(math.sqrt(sum(x*x for x in delta)) <= 1e-12,
              'atomic position differs by >1e-12 Bohr')
    check(r['convergence']['converged'] is True, 'our SCF is not certified converged')
    if failures:
        raise Error('conditions mismatch:\n- ' + '\n- '.join(failures), 3)
    return [] if grid_keys else ['QE automatic FFT grid: grid equality not required; inspect grid convergence']


def compare(ours, bundle, tolerances=None, decomposition_note=None):
    source, r = results(ours)
    bundle = Path(bundle)
    p = provenance(bundle, r)
    try:
        output_text = (bundle / 'qe.out').read_text()
        q = parse_output(output_text)
    except OSError as exc:
        raise Error(f'expected readable qe.out: {exc}', 3) from exc
    if q['version'] != p['qe_version'].strip():
        raise Error('provenance mismatch: qe_version versus Program PWSCF line', 2)
    warnings = conditions(r, q, (bundle / 'qe.in').read_text())
    tol = DEFAULTS.copy()
    if tolerances:
        tol.update(tolerances)
    for value in tol.values():
        positive(value)
    checks = []
    def add(name, ours_value, qe_value, unit, tolerance, informational=False):
        ours_value, qe_value = number(ours_value), number(qe_value)
        diff = abs(ours_value - qe_value)
        checks.append(dict(item=name, ours=ours_value, qe=qe_value, absolute_difference=diff,
                           unit=unit, tolerance=tolerance, passed=diff <= tolerance,
                           informational=informational))
    energy = r['energies']
    add('total energy per atom', energy['total']['ha']/len(r['atoms']), q['total']/len(r['atoms']), 'Ha/atom', tol['energy'])
    terms = dict(one_electron=sum(energy[k]['ha'] for k in ('kinetic', 'local', 'nonlocal')),
                 **{k: energy[k]['ha'] for k in ('hartree', 'xc', 'ewald')})
    for key, value in terms.items():
        add(key, value, q[key + '_ha'], 'Ha', tol['term'], bool(decomposition_note))
    eig = r['eigenvalues']['ha']
    if len(eig) != r['bands'] or any(number(x) > number(y) for x, y in zip(eig, eig[1:])):
        raise Error('expected sorted our eigenvalues of length bands', 3)
    for i, (x, y) in enumerate(zip(eig, q['eigenvalues_ha'])):
        add(f'eigenvalue band {i+1}', x, y, 'Ha', tol['eigenvalue'])
    occupied = int(r['electrons']/2)
    if occupied < r['bands']:
        add('HOMO-LUMO gap', eig[occupied]-eig[occupied-1],
            q['eigenvalues_ha'][occupied]-q['eigenvalues_ha'][occupied-1], 'Ha', tol['gap'])
    else:
        warnings.append('No empty band: HOMO-LUMO gap unavailable')
    # QE prints configured electrons, not an independent integrated density.
    add('electron count', r['electrons'], q['electrons'], 'electrons', tol['electrons']*abs(r['electrons']))
    last = r['convergence']['last_iteration']
    add('our integrated relative charge error', 0, last['relative_charge_error'], 'relative', tol['electrons'])
    add('our occupied eigen residual', 0, last['max_eigen_residual_ha'], 'Ha', 1e-7)
    passed = all(c['passed'] or c['informational'] for c in checks)
    maxima = {}
    for group, subset in [('eigenvalues', [c for c in checks if c['item'].startswith('eigenvalue')]),
                          ('energy terms', [c for c in checks if c['item'] in terms])]:
        maxima[group] = max(subset, key=lambda c: c['absolute_difference'])
    header = 'QE comparison'
    if 'synthetic self-test, not a reference' in p['notes']:
        header += ' — synthetic self-test, not a reference'
    if p['upf_sha256'] == UNVERIFIED_SHA:
        header += '\nPseudopotential provenance is unverified'
    worst = max(checks, key=lambda c:c['absolute_difference']/c['tolerance'])
    return dict(header=header, passed=passed, provenance=p, ours_results_sha256=sha256(source),
                qe_output_sha256=sha256(bundle/'qe.out'), qe=q, ours_grid=r['basis']['fft_grid'],
                qe_grid=q['grid'], warnings=warnings, tolerances=tol,
                maximum_tolerance_ratio=dict(item=worst['item'], ratio=worst['absolute_difference']/worst['tolerance']),
                nondefault_tolerances={k:v for k,v in tol.items() if v != DEFAULTS[k]},
                decomposition_note=decomposition_note, checks=checks, maximum_differences=maxima,
                limitations=['Orthogonality is not exported by results.json or this QE stdout protocol; not assessed.',
                             'QE electron count is configured count, not an independent density integral.',
                             '1 Ha = 2 Ry exactly; 1 Ha = 27.211386245988 eV (CODATA 2018). QE constant differences at 1e-9 relative are irrelevant to default eigenvalue tolerances.'])


def markdown(report):
    lines = ['# ' + report['header'], '', 'PASS' if report['passed'] else 'FAIL', '',
             '## Provenance', '']
    lines += [f'- {k}: {v}' for k, v in report['provenance'].items()]
    lines += ['', f"QE functional: {report['qe']['xc']}",
              f"FFT grids: ours {report['ours_grid']}; QE {report['qe_grid']}",
              f"Tolerances used: {report['tolerances']}",
              f"Nondefault tolerances: {report['nondefault_tolerances']}",
              f"Informational decomposition reason: {report['decomposition_note']}", '',
              '| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |',
              '| --- | ---: | ---: | ---: | --- | ---: | --- |']
    for c in report['checks']:
        status = ('PASS' if c['passed'] else 'FAIL') + (' (informational)' if c['informational'] else '')
        lines.append(f"| {c['item']} | {c['ours']:.17g} | {c['qe']:.17g} | {c['absolute_difference']:.9g} | {c['unit']} | {c['tolerance']:.9g} | {status} |")
    lines += ['', f"Our results SHA-256: {report['ours_results_sha256']}",
              f"QE output SHA-256: {report['qe_output_sha256']}", 'Maximum differences:']
    lines += [f"- {k}: {v['item']}, {v['absolute_difference']:.9g} {v['unit']}" for k,v in report['maximum_differences'].items()]
    lines += ['', *report['warnings'], *report['limitations'], '']
    return '\n'.join(lines)


def main(argv=None):
    parser = Parser(description=__doc__)
    parser.add_argument('--ours', required=True)
    parser.add_argument('--bundle', required=True)
    parser.add_argument('--json', default='qe-comparison.json')
    parser.add_argument('--markdown', default='qe-comparison.md')
    parser.add_argument('--decomposition-note', help='Explicit reason to treat all decomposition terms as informational')
    for key, value in DEFAULTS.items():
        parser.add_argument('--tol-' + key, type=positive, default=value)
    try:
        args = parser.parse_args(argv)
        report = compare(args.ours, args.bundle, {k: getattr(args, 'tol_' + k) for k in DEFAULTS}, args.decomposition_note)
        protected = {Path(args.ours).resolve(), *( (Path(args.bundle)/x).resolve() for x in ('qe.in', 'qe.out', 'provenance.json'))}
        if Path(args.ours).is_dir():
            protected.add((Path(args.ours)/'results.json').resolve())
        if Path(args.json).resolve() in protected or Path(args.markdown).resolve() in protected or Path(args.json).resolve() == Path(args.markdown).resolve():
            raise Error('report paths must be distinct and must not overwrite evidence')
        write_json(args.json, report)
        Path(args.markdown).write_text(markdown(report))
        print(report['header'] + ': ' + ('PASS' if report['passed'] else 'FAIL'))
        return 0 if report['passed'] else 5
    except Error as exc:
        print(exc, file=sys.stderr)
        return exc.code
    except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        print(f'invalid or missing input: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
