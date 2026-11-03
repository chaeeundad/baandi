#!/usr/bin/env python3
"""Emit complex-wavefunction QE SCF inputs; never inspect pseudopotential contents."""
import sys
import json
import math
from pathlib import Path
# Explicit sibling import also works under python3 -I.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_common import Error, Parser, fmt, positive, quoted, results, sha256, write_json


def functional_input(r):
    """Explicitly select PBE even when the pseudopotential was made with LDA."""
    if r['xc'] == 'pbe':
        return [
            '! Requested XC: pbe; UPF functional: ' + json.dumps(r['species'][0]['functional'], ensure_ascii=True),
            '! Explicit functional override: allow_functional_override=' +
            str(r['allow_functional_override']).lower() + '; functional_override_applied=' +
            str(r['functional_override_applied']).lower(),
            "  input_dft = 'pbe',",
        ]
    if r['xc'] == 'lda-pz':
        # Preserve legacy input unless the user explicitly overrode a UPF label.
        return ["  input_dft = 'PZ',"] if r.get('allow_functional_override', False) else []
    raise Error('unsupported XC functional')


def kpoint_input(r):
    """Return QE card, noinv flag and the expected canonical weighted set."""
    k = r.get('kpoints')
    if k is None:
        return ['K_POINTS automatic', '1 1 1 0 0 0'], True, [dict(fractional=[0, 0, 0], weight=1)]
    try:
        points = k['points']
        if not points or k['count'] != len(points) or type(k['time_reversal']) is not bool:
            raise Error('invalid kpoint count/time_reversal')
        expected = []
        for p in points:
            xyz = p['fractional']
            if len(xyz) != 3:
                raise Error('expected three k coordinates')
            expected.append(dict(fractional=[float(fmt(x)) for x in xyz], weight=positive(p['weight'])))
        if abs(math.fsum(p['weight'] for p in expected) - 1) > 1e-14:
            raise Error('kpoint weights must sum to 1 without normalization')
        if k['grid'] is not None:
            grid, shift = k['grid'], k['shift']
            if len(grid) != 3 or any(type(n) is not int or n <= 0 for n in grid):
                raise Error('expected positive kgrid integers')
            if len(shift) != 3 or any(s not in (0, 0.5) for s in shift):
                raise Error('kshift must be 0 or 0.5')
            card = ['K_POINTS automatic', ' '.join(map(str, grid + [int(2*s) for s in shift]))]
            return card, not k['time_reversal'], expected
        # Explicit points already include any requested CLI time-reversal reduction.
        # Disable all further reduction so weights and the supplied set are retained.
        return ['K_POINTS crystal', str(len(expected))] + [
            ' '.join(map(fmt, p['fractional'] + [p['weight']])) for p in expected], True, expected
    except (KeyError, TypeError) as exc:
        raise Error(f'malformed kpoints block: {exc}') from exc


def compare_kpoints(expected, observed, time_reversal=False, tolerance=1e-10):
    """Compare normalized QE crystal-coordinate lists, independent of ordering.

    Callers must convert QE Cartesian coordinates to fractional reciprocal units
    and QE's spin weights (sum 2) to integration weights (sum 1) explicitly.
    Returns False for mismatched multiplicities, weights or coordinates.
    """
    for points in (expected, observed):
        if not points:
            return False
        for p in points:
            if len(p['fractional']) != 3 or any(not math.isfinite(x) for x in p['fractional']):
                return False
            if not math.isfinite(p['weight']) or p['weight'] <= 0:
                return False
        if abs(math.fsum(p['weight'] for p in points) - 1) > tolerance:
            return False
    remaining = list(observed)
    for p in expected:
        for i, q in enumerate(remaining):
            if abs(p['weight'] - q['weight']) > tolerance:
                continue
            signs = (1, -1) if time_reversal else (1,)
            if any(all(abs((a-sign*b)-round(a-sign*b)) <= tolerance
                       for a, b in zip(p['fractional'], q['fractional'])) for sign in signs):
                remaining.pop(i)
                break
        else:
            return False
    return not remaining


def generate(results_path, pseudo_dir, upf_name, out, fft_grid='ours', ecut_ry=None,
             conv_thr=1e-12, mixing_beta=0.7, prefix='rust_dft', outdir='./qe-work'):
    source, r = results(results_path)
    if fft_grid not in ('ours', 'auto'):
        raise Error('fft_grid must be ours or auto')
    if Path(upf_name).name != upf_name or upf_name in ('.', '..') or any(c.isspace() for c in upf_name):
        raise Error('upf-name must be a basename')
    for s in (str(pseudo_dir), upf_name, prefix, outdir):
        quoted(s)
    cutoff = positive(r['cutoff']['ecutwfc_ry'] if ecut_ry is None else ecut_ry)
    positive(conv_thr)
    if not 0 < positive(mixing_beta) <= 1:
        raise Error('QE mixing_beta must be <= 1')
    card, noinv, expected_kpoints = kpoint_input(r)
    species = r['species'][0]
    # The UPF is deliberately not opened: the source result carries its digest.
    lines = ["! Generated from rust-dft results; lattice columns transposed to QE rows.",
             "! Automatic 1x1x1 selects complex Gamma, not the real-wavefunction gamma trick.",
             '&CONTROL', "  calculation = 'scf',", f'  prefix = {quoted(prefix)},',
             f'  pseudo_dir = {quoted(str(pseudo_dir))},', f'  outdir = {quoted(outdir)},',
             "  disk_io = 'low',", "  verbosity = 'high',", '  tprnfor = .false.,',
             '  tstress = .false.,', '/', '&SYSTEM', '  ibrav = 0,',
             f"  nat = {len(r['atoms'])},", '  ntyp = 1,',
             f'  ecutwfc = {fmt(cutoff)},', f'  ecutrho = {fmt(4 * cutoff)},',
             f"  nbnd = {r['bands']},", "  occupations = 'fixed',",
             '  nosym = .true.,']
    lines += functional_input(r)
    if 'kpoints' in r:
        lines[1] = '! Gamma-centered (m+s)/n grid; only requested time reversal, no space-group reduction.'
    if noinv:
        lines += ['  noinv = .true.,']
    if fft_grid == 'ours':
        lines += [f'  nr{i+1} = {n},' for i, n in enumerate(r['basis']['fft_grid'])]
    lines += ['/', '&ELECTRONS', "  diagonalization = 'david',", '  diago_full_acc = .true.,',
              f'  conv_thr = {fmt(conv_thr)},', f'  mixing_beta = {fmt(mixing_beta)},',
              '  electron_maxstep = 200,', '/', 'ATOMIC_SPECIES',
              f"{species['name']} 1.0 {upf_name}", 'CELL_PARAMETERS bohr']
    a = r['cell']['lattice_columns_bohr']
    lines += [' '.join(fmt(a[i][j]) for i in range(3)) for j in range(3)]
    lines += ['ATOMIC_POSITIONS crystal']
    lines += [atom['species'] + ' ' + ' '.join(map(fmt, atom['fractional'])) for atom in r['atoms']]
    lines += card
    out = Path(out)
    protected = {source.resolve(), (Path(pseudo_dir)/upf_name).resolve()}
    if out.resolve() in protected or Path(str(out)+'.meta.json').resolve() in protected:
        raise Error('output must not overwrite source results or pseudopotential')
    out.write_text('\n'.join(lines) + '\n')
    meta = dict(input_sha256=sha256(out), source_results_sha256=sha256(source),
                upf_name=upf_name, upf_sha256=species['sha256'],
                generator_options=dict(results=str(source), pseudo_dir=str(pseudo_dir),
                    upf_name=upf_name, out=str(out), fft_grid=fft_grid, ecut_ry=ecut_ry,
                    conv_thr=conv_thr, mixing_beta=mixing_beta, prefix=prefix, outdir=outdir))
    meta['expected_kpoints'] = expected_kpoints
    meta['kpoint_comparison'] = 'Compare modulo reciprocal integers and time reversal when enabled; weights sum to one.'
    write_json(str(out) + '.meta.json', meta)
    return meta


def main(argv=None):
    p = Parser(description=__doc__)
    for arg in ('results', 'pseudo-dir', 'upf-name', 'out'):
        p.add_argument('--' + arg, required=True)
    p.add_argument('--fft-grid', choices=['ours', 'auto'], default='ours')
    p.add_argument('--ecut-ry', type=positive)
    p.add_argument('--conv-thr', type=positive, default=1e-12)
    p.add_argument('--mixing-beta', type=positive, default=0.7)
    p.add_argument('--prefix', default='rust_dft')
    p.add_argument('--outdir', default='./qe-work')
    try:
        args = vars(p.parse_args(argv))
        args['results_path'] = args.pop('results')
        generate(**args)
        return 0
    except (Error, OSError, ValueError, TypeError) as exc:
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
