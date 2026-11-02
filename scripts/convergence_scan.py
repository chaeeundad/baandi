#!/usr/bin/env python3
"""Run independent convergence series and emit paired QE inputs; stdlib only."""
import csv
import itertools
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_common import Error, Parser, EV_PER_HA, positive, results, sha256, write_json
from make_qe_input import generate


def edit_table(text, table, updates):
    """Preserve source text, replacing only single-line values in a named table.

    Deliberately not a general TOML parser. The CLI validates the resulting TOML.
    Edited tables must use the documented bare keys and single-line values.
    """
    lines = text.splitlines()
    heads = [i for i, line in enumerate(lines) if re.fullmatch(r'\s*\[' + re.escape(table) + r'\]\s*(?:#.*)?', line)]
    if len(heads) > 1:
        raise Error(f'duplicate [{table}]')
    if not heads:
        if not any(value is not None for value in updates.values()):
            return text
        lines += ['', '[' + table + ']']
        heads = [len(lines)-1]
    start = heads[0] + 1
    end = next((i for i in range(start, len(lines)) if lines[i].lstrip().startswith('[')), len(lines))
    block = lines[start:end]
    found = set()
    replacement = []
    for line in block:
        if not line.strip() or line.lstrip().startswith('#'):
            replacement.append(line)
            continue
        m = re.fullmatch(r'\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*', line)
        if not m:
            raise Error(f'expected bare key and single-line value in [{table}]: {line}')
        key, value = m.groups()
        if key in found:
            raise Error(f'duplicate {table}.{key}')
        found.add(key)
        if key in updates:
            # Reject multiline arrays/strings rather than leaving continuation debris.
            simple = value.split('#', 1)[0].strip()
            if simple.count('[') != simple.count(']') or "'''" in simple or '"""' in simple:
                raise Error(f'expected single-line {table}.{key}')
            if updates[key] is not None:
                replacement.append(key + ' = ' + updates[key])
        else:
            replacement.append(line)
    replacement += [key + ' = ' + value for key, value in updates.items() if key not in found and value is not None]
    return '\n'.join(lines[:start] + replacement + lines[end:]) + '\n'


def absolute_pseudos(text, base):
    lines = text.splitlines()
    inside = False
    found = False
    for i, line in enumerate(lines):
        if line.lstrip().startswith('['):
            inside = bool(re.fullmatch(r'\s*\[pseudopotentials\]\s*(?:#.*)?', line))
            found |= inside
        elif inside and line.strip() and not line.lstrip().startswith('#'):
            m = re.fullmatch(r'''\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*("(?:[^"\\]|\\.)*"|'[^']*')\s*(?:#.*)?''', line)
            if not m:
                raise Error('expected single-line quoted pseudopotential path with bare species key')
            value = json.loads(m.group(2)) if m.group(2).startswith('"') else m.group(2)[1:-1]
            lines[i] = m.group(1) + ' = ' + json.dumps(str((base / value).resolve()))
    if not found:
        raise Error('expected [pseudopotentials]')
    return '\n'.join(lines) + '\n'


def smooth_at_least(n):
    def smooth(x):
        for factor in (2, 3, 5):
            while x % factor == 0:
                x //= factor
        return x == 1
    while not smooth(n):
        n += 1
    return n


def minimum_grid(r):
    """Enumerate density sphere inside the first run's automatic grid bounds."""
    b = r['cell']['reciprocal_columns_inverse_bohr']
    dims = r['basis']['fft_grid']
    required = [1, 1, 1]
    for m in itertools.product(*(range(-(n//2), n//2+1) for n in dims)):
        norm2 = sum(sum(b[i][j]*m[j] for j in range(3))**2 for i in range(3))
        if norm2 <= 4*r['cutoff']['ecutwfc_ry']*(1+1e-12):
            required = [max(n, 2*abs(v)+1) for n, v in zip(required, m)]
    return required


def run_scan(args):
    source = Path(args.input).resolve()
    binary = Path(args.dft_cli).resolve()
    root = Path(args.out_dir).resolve()
    # Exclusive destination prevents stale successes and accidental overwrites.
    root.mkdir(parents=True, exist_ok=False)
    original = source.read_text()
    base = absolute_pseudos(original, source.parent)
    if args.emit_qe_inputs and (not args.pseudo_dir or not args.upf_name):
        raise Error('--emit-qe-inputs requires --pseudo-dir and --upf-name')
    manifest = dict(source_input=str(source), source_input_sha256=sha256(source),
                    driver_command=[sys.executable, '-I', str(Path(__file__).resolve()), *sys.argv[1:]],
                    points=[], qe_inputs=[], status='running')
    write_json(root/'manifest.json', manifest)
    rows = []

    def point(series, label, updates, rank):
        point_id = f'{len(rows):03d}-{series}-{label}'
        folder = root/point_id
        folder.mkdir()
        text = base
        for table, fields in updates.items():
            text = edit_table(text, table, fields)
        inp = folder/'input.toml'
        inp.write_text(text)
        command = [str(binary), 'run', str(inp), '--output-dir', str(folder)]
        record = dict(point_id=point_id, command=command, input_sha256=sha256(inp),
                      input=str(inp), results=str(folder/'results.json'), status='running')
        manifest['points'].append(record)
        write_json(root/'manifest.json', manifest)
        start = time.monotonic()
        try:
            with (folder/'stdout.txt').open('w') as stdout, (folder/'stderr.txt').open('w') as stderr:
                process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=args.timeout, check=False)
            record['exit_code'] = process.returncode
            record['elapsed_wall_s'] = time.monotonic()-start
            if process.returncode != 0:
                raise Error(f'{point_id}: dft-cli exited {process.returncode}; see {folder}/stderr.txt')
            result_path, r = results(folder)
            if r['convergence']['converged'] is not True:
                raise Error(f'{point_id}: expected certified convergence')
            if r['input_sha256'] != record['input_sha256']:
                raise Error(f'{point_id}: input hash mismatch')
            record.update(status='done', results_sha256=sha256(result_path), commit_hash=r['tool']['git_commit'])
            occ = int(r['electrons']/2)
            eig = r['eigenvalues']['ha']
            gap = (eig[occ]-eig[occ-1])*EV_PER_HA if occ < r['bands'] else None
            rows.append(dict(point_id=point_id, scan=series, parameter=label, rank=rank,
                             total_energy_ha=r['energies']['total']['ha'],
                             energy_per_atom_ha=r['energies']['total']['ha']/len(r['atoms']),
                             difference_to_finest_ha=None, gap_ev=gap,
                             iterations=r['convergence']['iterations'], wall_time_s=r['timings']['total_s']))
            if args.emit_qe_inputs:
                for mode in ('ours', 'auto'):
                    out = folder/f'qe-{mode}.in'
                    meta = generate(result_path, args.pseudo_dir, args.upf_name, out,
                                    fft_grid=mode, prefix=point_id + '-' + mode)
                    manifest['qe_inputs'].append(dict(input=str(out), sha256=meta['input_sha256'],
                        fft_grid=mode, our_results=str(result_path), our_results_sha256=sha256(result_path),
                        upf_sha256=meta['upf_sha256']))
            write_json(root/'manifest.json', manifest)
            return r
        except (OSError, subprocess.TimeoutExpired, Error, ValueError, KeyError, TypeError, IndexError) as exc:
            record.update(status='failed', error=str(exc), elapsed_wall_s=time.monotonic()-start)
            manifest['status'] = 'failed'
            write_json(root/'manifest.json', manifest)
            raise Error(f'scan failed: {exc}') from exc

    if args.scan in ('ecut', 'all'):
        for cutoff in sorted(set(positive(x) for x in args.ecuts.split(','))):
            point('ecut', f'{cutoff:g}Ry', {'basis': {'ecutwfc_ry': repr(cutoff), 'fft_grid': None}}, cutoff)
    if args.scan in ('grid', 'all'):
        r = point('grid', 'automatic', {'basis': {'fft_grid': None}}, 0)
        dims = r['basis']['fft_grid']
        rows[-1]['rank'] = math.prod(dims)
        if args.grids:
            grids = [[int(x) for x in item.split('x')] for item in args.grids.split(',')]
            if any(len(g) != 3 or min(g) <= 0 for g in grids):
                raise Error('--grids requires comma-separated positive NxNxN triples')
        else:
            grids = [minimum_grid(r)] + [[smooth_at_least(math.ceil(f*n)) for n in dims] for f in (1.25, 1.5)]
        # Keep explicit minimal even if it equals automatic: records the override path.
        for i, grid in enumerate(grids):
            point('grid', f'{i}-' + 'x'.join(map(str, grid)), {'basis': {'fft_grid': json.dumps(grid)}}, math.prod(grid))
    if args.scan in ('tolerance', 'all'):
        for density in (1e-5, 1e-6, 1e-7, 1e-8, 1e-9):
            point('tolerance', f'{density:g}', {'scf': {
                'density_residual_tolerance': repr(density),
                'energy_tolerance_ha': repr(density*0.1),
                'eigen_tolerance': repr(density*0.1)}}, -density)
    for row in rows:
        finest = max((x for x in rows if x['scan'] == row['scan']), key=lambda x: x['rank'])
        row['difference_to_finest_ha'] = row['total_energy_ha']-finest['total_energy_ha']
        row['finest_point_id'] = finest['point_id']
    write_json(root/'scan.json', rows)
    with (root/'scan.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    keys = ['point_id', 'parameter', 'total_energy_ha', 'energy_per_atom_ha',
            'difference_to_finest_ha', 'gap_ev', 'iterations', 'wall_time_s']
    lines = ['# Convergence scan', '', 'No QE output or reference values are included.', '',
             '| ' + ' | '.join(keys) + ' |', '| ' + ' | '.join(['---']*len(keys)) + ' |']
    lines += ['| ' + ' | '.join(str(row[key]) for key in keys) + ' |' for row in rows]
    (root/'scan.md').write_text('\n'.join(lines)+'\n')
    manifest['status'] = 'done'
    write_json(root/'manifest.json', manifest)
    return rows


def main(argv=None):
    p = Parser(description=__doc__)
    for key in ('input', 'dft-cli', 'out-dir', 'scan'):
        p.add_argument('--'+key, required=True, **({'choices': ['ecut', 'grid', 'tolerance', 'all']} if key == 'scan' else {}))
    p.add_argument('--ecuts', default='20,25,30,40,50,60,80')
    p.add_argument('--grids', help='e.g. 24x24x24,30x30x30')
    p.add_argument('--timeout', type=positive, default=14400)
    p.add_argument('--emit-qe-inputs', action='store_true')
    p.add_argument('--pseudo-dir')
    p.add_argument('--upf-name')
    try:
        run_scan(p.parse_args(argv))
        return 0
    except (Error, OSError, ValueError, KeyError, TypeError) as exc:
        print(exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
