#!/usr/bin/env python3
"""Emit complex-Gamma QE SCF inputs; never inspect pseudopotential contents."""
import sys
from pathlib import Path
# Explicit sibling import also works under python3 -I.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_common import Error, Parser, fmt, positive, quoted, results, sha256, write_json


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
             '  nosym = .true.,', '  noinv = .true.,']
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
    lines += ['K_POINTS automatic', '1 1 1 0 0 0']
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
