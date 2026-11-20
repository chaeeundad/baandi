"""Mac scientific identity and actual arrival JSON verification; no speed analysis."""
import argparse
import collections
import hashlib
import json
import runpy
from pathlib import Path

HERE = Path(__file__).resolve().parent
parse = runpy.run_path(str(HERE / 'wait_parse.py'))['parse']
exact = runpy.run_path(str(HERE / 'verify_cli.py'))['exact']


def records(path, prefix):
    return [json.loads(line[len(prefix):]) for line in path.read_text().splitlines() if line.startswith(prefix)]


def verify(root):
    arms = ('off', 'counts', 'timed')
    dav = [json.loads((root / f'arrival-dav-{a}.json').read_text()) for a in arms]
    keys = ('kind', 'eigen_fnv', 'eigenvalues', 'residuals', 'converged', 'h_columns',
            'davidson_iterations', 'counters')
    for d in dav:
        assert all(exact(d['kinds'][0][k]) == exact(dav[0]['kinds'][0][k]) for k in keys)
    assert dav[0]['phase_profile']['arrival_trace'] is None
    counts, timed = [d['phase_profile']['arrival_trace'] for d in dav[1:]]
    table = [parse(d) for d in dav[1:]]
    assert len(table[0]) == len(table[1]) > 0
    signature_keys = ('map_id', 'bucket', 'seq', 'op', 'site', 'payload_bytes', 'previous_stage',
                      'previous_stage_points', 'previous_stage_logical_bytes', 'previous_stage_shape', 'call_index')
    for a, b in zip(counts['ranks'], timed['ranks']):
        assert len(a['events']) == len(b['events']) > 0
        assert [[e[k] for k in signature_keys] for e in a['events']] == [[e[k] for k in signature_keys] for e in b['events']]
        assert all(e['enter_ns'] is None and e['exit_ns'] is None and e['previous_stage_ns'] is None for e in a['events'])
        assert all(e['enter_ns'] <= e['exit_ns'] for e in b['events'])
    scf = [records(root / f'arrival-scf-{a}.log', 'M2SCREEN_JSON ') for a in arms]
    keys = ('ranks', 'energy_bits', 'density_fnv', 'scf', 'davidson', 'h_columns')
    assert all(len(rows) == 3 for rows in scf)
    for rows in scf:
        assert [[r[k] for k in keys] for r in rows] == [[r[k] for k in keys] for r in scf[0]]
    cli = [json.loads((root / f'arrival-cli-{a}/results.json').read_text()) for a in arms]
    fields = ('total_energy', 'energies', 'eigenvalues', 'occupations', 'band_accuracy', 'certified',
              'verified', 'verification', 'basis', 'cell', 'atoms', 'bands', 'electrons', 'cutoff')
    density = [(root / f'arrival-cli-{a}/density.bin').read_bytes() for a in arms]
    assert all(d == density[0] for d in density)
    for d in cli:
        assert d['run']['completed'] and d['verified']
        assert all(exact(d[k]) == exact(cli[0][k]) for k in fields)
        assert exact(d['convergence']['history']) == exact(cli[0]['convergence']['history'])
    return dict(task_id='LANE-FFT', purpose='Mac correctness only; no performance conclusion',
        serial_energy_ha=cli[0]['total_energy']['ha'], serial_verified=True, serial_delta_ha=0,
        serial_density_sha256=hashlib.sha256(density[0]).hexdigest(),
        simulated_mpi=[{k:r[k] for k in keys} | {'energy_ha':r['energy_ha']} for r in scf[0]],
        simulated_mpi_off_count_timed_energy_delta_ha=0,
        mpi_density_evidence='FNV from existing SCF harness; raw density bytes checked on serial CLI',
        warm_eigen_fnv=dav[0]['kinds'][0]['eigen_fnv'],
        warm_eigenvalues_residuals_work_equal=True, all_trace_signatures_equal=True,
        trace_ranks=timed['rank_count'], events_per_rank=[len(r['events']) for r in timed['ranks']],
        fine_parser_rows=len(table[1]), clock_error_ns=timed['clock_error_ns'],
        collective_kinds=dict(collections.Counter(e['op'] for e in timed['ranks'][0]['events'])),
        sites=dict(collections.Counter(e['site'] for e in timed['ranks'][0]['events'])),
        raw_json_sha256={a:hashlib.sha256((root/f'arrival-dav-{a}.json').read_bytes()).hexdigest() for a in arms})


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('scratch', type=Path)
    args = ap.parse_args()
    print(json.dumps(verify(args.scratch), indent=2))
