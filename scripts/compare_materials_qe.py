#!/usr/bin/env python3
"""Report BAANDI minus QE differences, without deciding MF05/MF06 acceptance caps.

Exit 0: report produced (NOT scientific acceptance); 3: inconsistent evidence;
1: usage/I/O/unsupported input. No pseudopotential files are opened.
"""
import json
import math
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_common import (Error, Parser, NUMBER, EV_PER_HA, materials_results,
                       number, sha256)
from make_qe_input import kpoint_input
from verify_qe import match, require_strict_result
from qe_rhog import load_ours, read_qe, validate_grid, FormatError
from compare_rhog import compare as compare_density


def require(ok, message):
    if not ok:
        raise Error(message, 3)


def parse_output(text):
    """Literal QE stdout records; reject partial, repeated or ambiguous SCFs."""
    n = '(' + NUMBER + ')'
    q = {'version': match(text, r'^\s*Program PWSCF\s+v\.(\S+)\s+starts[^\n]*$', 'version').group(1)}
    for key, pattern, cast in [
        ('electrons', r'number of electrons\s*=\s*' + n, number),
        ('bands', r'number of Kohn-Sham states\s*=\s*(\d+)', int),
        ('nk', r'number of k points\s*=\s*(\d+)', int),
        ('ecutwfc', r'kinetic-energy cutoff\s*=\s*' + n + r'\s+Ry', number),
        ('ecutrho', r'charge density cutoff\s*=\s*' + n + r'\s+Ry', number),
        ('iterations', r'convergence has been achieved in\s+(\d+)\s+iterations', int),
        ('energy_ry', r'!\s+total energy\s*=\s*' + n + r'\s+Ry', number),
        ('xc', r'Exchange-correlation\s*=\s*([^\n]+)', str),
    ]:
        # QE appends smearing information to its k-count record; fixed only here.
        if key == 'xc':
            # QE 7.6 prints this record twice when input_dft is set (once in the
            # enforced-functional banner, once in the summary). Accept repeats only
            # when every value is identical; differing values stay an error.
            found = match(text, r'^\s*' + pattern + r'\s*$', key, multiple=True)
            values = {m.group(1).strip() for m in re.finditer(
                r'^\s*' + pattern + r'\s*$', text, re.MULTILINE)}
            require(len(values) == 1, f'conflicting {key} records: {sorted(values)}')
            q[key] = found.group(1)
            continue
        q[key] = cast(match(text, r'^\s*' + pattern + r'\s*$', key).group(1))
    match(text, r'^\s*JOB DONE\.\s*$', 'JOB DONE')
    g = match(text, r'^\s*Dense\s+grid:[^\n]*FFT dimensions:\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*$', 'FFT grid')
    q['grid'] = list(map(int, g.groups()))
    end = match(text, r'^\s*End of self-consistent calculation\s*$', 'final SCF')
    row = re.compile(r'^\s*k\(\s*(\d+)\)\s*=\s*\(\s*' + n + r'\s+' + n + r'\s+' + n + r'\s*\),\s*wk\s*=\s*' + n + r'\s*$', re.M)
    def table(marker):
        start = match(text[:end.start()], marker, 'k coordinate table')
        rows = []
        for line in text[start.end():end.start()].splitlines():
            if not line.strip():
                if rows:
                    break
                continue
            m = row.fullmatch(line)
            if not m:
                break
            require(int(m[1]) == len(rows) + 1, 'k table indices not consecutive')
            rows.append([number(x) for x in m.groups()[1:]])
        require(len(rows) == q['nk'], 'incomplete k coordinate table; use verbosity high')
        return rows
    q['crystal'] = table(r'^\s*cryst\. coord\.\s*$')
    cart = table(r'^\s*cart\. coord\. in units 2pi/alat\s*$')
    final = text[end.end():]
    header = re.compile(r'^\s*k\s*=\s*' + n + r'\s*' + n + r'\s*' + n + r'\s*\(\s*\d+\s+PWs\)\s+bands \(ev\):\s*$', re.M)
    blocks = list(header.finditer(final))
    require(len(blocks) == q['nk'], 'final band block count mismatch')
    q['eigenvalues_ha'] = []
    for i, block in enumerate(blocks):
        require(all(abs(number(x)-y) <= 1e-4 for x, y in zip(block.groups(), cart[i][:3])),
                'final bands/cartesian k order mismatch')
        values = []
        for line in final[block.end():].splitlines():
            if not line.strip():
                if values:
                    break
                continue
            tokens = line.split()
            require(all(re.fullmatch(NUMBER, x) for x in tokens), 'malformed eigenvalue tokens')
            values.extend(number(x)/EV_PER_HA for x in tokens)
            if len(values) >= q['bands']:
                break
        require(len(values) == q['bands'], 'incomplete eigenvalues')
        require(all(x <= y for x, y in zip(values, values[1:])), 'unsorted eigenvalues')
        q['eigenvalues_ha'].append(values)
    require(all(k[3] > 0 for k in q['crystal']), 'nonpositive k weights')
    require(abs(math.fsum(k[3] for k in q['crystal'])-2) <= q['nk']*5e-7,
            'QE spin-degenerate weights must sum to two within printed precision')
    return q


def difference(ours, qe, unit):
    ours, qe = number(ours), number(qe)
    return dict(ours=ours, qe=qe, difference=ours-qe,
                absolute_difference=abs(ours-qe), unit=unit)


def compare(ours, qe_out, qe_density=None):
    source, r = materials_results(ours)
    require_strict_result(r)
    q = parse_output(Path(qe_out).read_text())
    require(r['convergence']['converged'] is True, 'BAANDI SCF is not converged')
    require(q['grid'] == r['basis']['fft_grid'], 'FFT grid mismatch')
    require(q['bands'] == r['bands'], 'band count mismatch')
    require(q['electrons'] == r['electrons'], 'configured electron count mismatch')
    require(q['ecutwfc'] == r['cutoff']['ecutwfc_ry'] and q['ecutrho'] == 4*q['ecutwfc'], 'cutoff mismatch')
    expected_xc = ('SLA PW PBX PBC', 'PBE') if r['xc'] == 'pbe' else ('SLA PZ NOGX NOGC', 'SLA PZ')
    label = ' '.join(q['xc'].split('(')[0].split())
    require(label in expected_xc, 'XC functional mismatch')
    _, noinv, expected = kpoint_input(r)
    points = r.get('kpoints', {}).get('points', [dict(fractional=[0, 0, 0], weight=1, eigenvalues=r['eigenvalues'])])
    require(len(expected) == q['nk'], 'k point count mismatch')
    # Fixed tolerances below describe stdout rounding, never scientific acceptance.
    remaining = set(range(q['nk']))
    pairs = []
    for p in expected:
        hits = [i for i in remaining if abs(p['weight']-q['crystal'][i][3]/2) <= 5e-7 and
                any(all(abs((a-sign*b)-round(a-sign*b)) <= 1e-6
                        for a, b in zip(p['fractional'], q['crystal'][i][:3]))
                    for sign in ((1,) if noinv else (1, -1)))]
        require(len(hits) == 1, 'ambiguous or mismatched k point/weight')
        i = hits[0]
        remaining.remove(i)
        pairs.append(i)
    occupied = int(r['electrons']/2)
    eigenvalues = []
    ours_spectra = []
    for index, (p, qi) in enumerate(zip(points, pairs)):
        values = [number(x) for x in p['eigenvalues']['ha']]
        require(len(values) == r['bands'] and all(x <= y for x, y in zip(values, values[1:])), 'invalid BAANDI eigenvalues')
        ours_spectra.append(values)
        for band in range(min(occupied+1, r['bands'])):
            eigenvalues.append(dict(k_index=index+1, qe_k_index=qi+1, fractional=p['fractional'],
                band=band+1, role='occupied' if band < occupied else 'first_empty',
                **difference(values[band], q['eigenvalues_ha'][qi][band], 'Ha')))
    def gaps(spectra):
        return dict(direct=min(e[occupied]-e[occupied-1] for e in spectra),
                    indirect=min(e[occupied] for e in spectra)-max(e[occupied-1] for e in spectra))
    gap = None
    if occupied < r['bands']:
        og, qg = gaps(ours_spectra), gaps(q['eigenvalues_ha'])
        gap = {key: difference(og[key], qg[key], 'Ha') for key in og}
    total = number(r['total_energy']['ha'])
    report = dict(status='REPORTED', acceptance='Not assessed; MF05/MF06 caps are frozen separately.',
        tolerances=None, difference_convention='BAANDI minus QE; eigenvalues use absolute energy zero without alignment',
        total_energy_ha=difference(total, q['energy_ry']/2, 'Ha'),
        total_energy_ry=difference(2*total, q['energy_ry'], 'Ry'),
        eigenvalues=eigenvalues, band_gap=gap,
        electrons=difference(r['electrons'], q['electrons'], 'electrons'),
        iterations=dict(ours=r['convergence']['iterations'], qe=q['iterations']),
        source=dict(ours_results_sha256=sha256(source), qe_output_sha256=sha256(qe_out),
                    baandi_commit=r['tool']['git_commit'], input_sha256=r['input_sha256'],
                    species=r['species'], qe_version=q['version']),
        density=None, limitations=[
            'QE stdout electron count is configured, not an independently integrated density.',
            'Stdout k matching uses rounding bounds 1e-6 coordinates and 5e-7 normalized weights.',
            'Stdout eigenvalues are rounded in eV; conversion uses 27.211386245988 eV/Ha.',
            'File hashes identify evidence; remote command, input and UPF provenance must be archived separately.',
            'No empty band means band gap is unavailable.'])
    if qe_density is not None:
        og, meta, validation = load_ours(source.parent)
        qg, metadata = read_qe(qe_density)
        validate_grid(meta['fft_grid'], qg)
        require(og.keys() == qg.keys(), 'density G sets mismatch')
        # Reuse established metrics with zero thresholds, then remove all gate fields.
        # No threshold or PASS result from this helper is an acceptance decision.
        density = compare_density(og, qg, dict(coefficient=0, rho0=0, parseval=0, relative_l2=0))
        for key in ('tolerances', 'non_default_tolerances', 'passed', 'exit_code', 'failures'):
            density.pop(key)
        density.update(validation=validation, qe_metadata=metadata, qe_sha256=sha256(qe_density),
                       qe_integral_electrons=qg[(0, 0, 0)].real*meta['volume_bohr3'])
        report['density'] = density
    return report


def main(argv=None):
    p = Parser(description=__doc__)
    p.add_argument('--ours', required=True)
    p.add_argument('--qe-out', required=True)
    p.add_argument('--qe-density', help='Optional charge-density.dat; requires sibling BAANDI density.json/bin')
    p.add_argument('--json', help='New report path; default stdout')
    try:
        a = p.parse_args(argv)
        report = compare(a.ours, a.qe_out, a.qe_density)
        text = json.dumps(report, indent=2, allow_nan=False) + '\n'
        if a.json:
            with Path(a.json).open('x') as stream:
                stream.write(text)
        else:
            print(text, end='')
        return 0
    except Error as exc:
        print(exc, file=sys.stderr)
        return exc.code
    except FormatError as exc:
        print(exc, file=sys.stderr)
        return 3
    except (OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        print(f'invalid or missing input: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
