"""Verify Phase-2 Mac correctness files; deliberately discard all timing data."""
import argparse
import hashlib
import json
import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
exact = runpy.run_path(str(Path(__file__).with_name('verify_cli.py')))['exact']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records(path, prefix):
    return [json.loads(line[len(prefix):]) for line in path.read_text().splitlines()
            if line.startswith(prefix)]


def verify(scratch, upf):
    arms = ('off', 'on')
    cli = [json.loads((scratch / f'fft-phase2-cli-{arm}/results.json').read_text()) for arm in arms]
    fields = ('total_energy', 'energies', 'eigenvalues', 'occupations', 'band_accuracy',
              'certified', 'verified', 'verification', 'basis', 'cell', 'atoms',
              'bands', 'electrons', 'cutoff', 'input_sha256')
    for r in cli:
        assert r['run']['completed'] and r['verified']
        for field in fields:
            assert exact(r[field]) == exact(cli[0][field]), field
        assert exact(r['convergence']['history']) == exact(cli[0]['convergence']['history'])
    densities = [(scratch / f'fft-phase2-cli-{arm}/density.bin').read_bytes() for arm in arms]
    assert densities[0] == densities[1]
    scf = [records(scratch / f'fft-phase2-scf-{arm}.log', 'M2SCREEN_JSON ') for arm in arms]
    assert len(scf[0]) == len(scf[1]) == 3
    keys = ('ranks', 'energy_ha', 'energy_bits', 'density_fnv', 'scf', 'davidson', 'h_columns')
    for a, b in zip(*scf):
        assert all(a[k] == b[k] for k in keys)
    assert [r['ranks'] for r in scf[0]] == [2, 3, 8]
    dav = [records(scratch / f'fft-phase2-dav-{arm}.log', 'DAVBENCH_JSON ') for arm in arms]
    assert all(len(r) == 1 for r in dav)
    dav = [r[0] for r in dav]
    a, b = [r['kinds'][0] for r in dav]
    for key in ('kind', 'eigen_fnv', 'eigenvalues', 'residuals', 'converged',
                'davidson_iterations', 'h_columns', 'counters'):
        assert exact(a[key]) == exact(b[key]), key
    assert a['kind'] == 'warm' and a['repeats_identical'] and b['repeats_identical']
    ranks = [r['phase_profile']['fft_lane']['per_rank'] for r in dav]
    assert len(ranks[0]) == len(ranks[1]) == 3
    counter_rows = []
    for ra, rb in zip(*ranks):
        assert ra['rank'] == rb['rank']
        for ba, bb in zip(ra['buckets'], rb['buckets']):
            assert ba['bucket'] == bb['bucket'] and ba['maps'] == bb['maps']
            removed = 0
            for sa, sb in zip(ba['stages'], bb['stages']):
                assert sa['stage'] == sb['stage']
                assert sa['seconds'] == sb['seconds'] == 0
                if sa['stage'] == 'kept_copy':
                    removed = sa['logical_bytes']
                    assert removed == 32 * sa['points']
                    assert sb['calls'] == sb['points'] == sb['logical_bytes'] == 0
                else:
                    assert sa == sb, sa['stage']
            assert sum(s['logical_bytes'] for s in ba['stages']) - sum(s['logical_bytes'] for s in bb['stages']) == removed
            if ba['bucket'] == 'fine':
                assert removed > 0
                counter_rows.append({'rank': ra['rank'], 'fine_removed_logical_bytes': removed})
    source_files = ['crates/dft-solvers/src/distributed_fft.rs', 'crates/dft-numerics/src/switches.rs']
    return {
        'task_id': 'LANE-FFT', 'phase': 2, 'purpose': 'Mac correctness only; no timing evidence',
        'cli_serial_control': {'energy_ha': cli[0]['total_energy']['ha'], 'delta_ha': 0,
            'verified': True, 'density_sha256': hashlib.sha256(densities[0]).hexdigest(),
            'scientific_fields_and_history_bit_identical': True,
            'scope': 'serial path unchanged; candidate active only in distributed M2'},
        'simulated_mpi_scf': [{k: r[k] for k in keys} for r in scf[0]],
        'simulated_mpi_off_on_delta_ha': 0,
        'simulated_mpi_rank_energy_spread_ha': max(r['energy_ha'] for r in scf[0]) - min(r['energy_ha'] for r in scf[0]),
        'scf_limit': 'M2SCREEN exports energy bits, density FNV and counts, not full density/history or verification fields',
        'fixed_h_real_upf': {'ranks': 3, 'kind': a['kind'], 'eigen_fnv': a['eigen_fnv'],
            'eigenvalues_residuals_and_counters_bit_identical': True, 'removed': counter_rows},
        'source_sha256': {p: sha(ROOT / p) for p in source_files},
        'binary_sha256': {p: sha(ROOT / p) for p in ['target/release/dft-cli', 'target/release/examples/davidson_bench']},
        'upf_sha256': sha(upf),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scratch', type=Path)
    parser.add_argument('--upf', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.scratch, args.upf), indent=2))
