"""Check actual OFF/timed/count-only CLI outputs. No performance analysis."""
import hashlib
import importlib.util
import json
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('lane_fft_parse', ROOT / 'scripts/lane_fft_parse.py')
parser = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parser)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(value):
    if isinstance(value, float):
        return struct.pack('>d', value).hex()
    if isinstance(value, list):
        return [exact(v) for v in value]
    if isinstance(value, dict):
        # Only reporting clocks are excluded, never a numerical tolerance.
        return {k: exact(v) for k, v in value.items() if k not in ('timings', 'wall_s')}
    return value


def verify(paths):
    reports = [json.loads((p / 'results.json').read_text()) for p in paths]
    fields = ['total_energy', 'energies', 'eigenvalues', 'occupations',
              'band_accuracy', 'certified', 'verified', 'verification',
              'basis', 'cell', 'atoms', 'bands', 'electrons', 'cutoff', 'input_sha256']
    for report in reports:
        assert report['run']['completed'] and report['verified']
        for field in fields:
            assert exact(report[field]) == exact(reports[0][field]), field
        assert exact(report['convergence']['history']) == exact(reports[0]['convergence']['history'])
    work_keys = ['certificate_davidson_iterations', 'certificate_h_vectors', 'certificate_solves',
                 'davidson_iterations_recorded', 'h_vectors_recorded', 'second_vector_columns']
    for report in reports:
        for key in work_keys:
            assert report['convergence']['work'][key] == reports[0]['convergence']['work'][key], key
    hashes = [sha(p / 'density.bin') for p in paths]
    assert len(set(hashes)) == 1, 'density bytes changed'
    timed, counts = [list(parser.rows(p / 'results.json')) for p in paths[1:]]
    assert len(timed) == len(counts) and timed
    keys = ['rank', 'bucket', 'maps', 'stage', 'calls', 'points', 'logical_bytes', 'fft_flops_estimate']
    assert [[r[k] for k in keys] for r in timed] == [[r[k] for k in keys] for r in counts]
    assert any(r['calls'] > 0 for r in counts), 'vacuous instrumentation'
    assert any(r['seconds'] > 0 for r in timed), 'timed mode did not read clocks'
    assert all(r['seconds'] is None for r in counts)
    summary = {'task_id': 'LANE-FFT', 'phase': 0, 'purpose': 'Mac correctness only; no performance conclusion',
        'max_energy_difference_ha': 0.0, 'energy_ha': reports[0]['total_energy']['ha'],
        'verified': True, 'scf_iterations': reports[0]['convergence']['iterations'],
        'density_sha256': hashes[0], 'scientific_fields_bit_identical': fields + ['convergence.history'],
        'parsed_rows_per_enabled_file': len(counts),
        'count_only_equals_timed_counters': True, 'timed_clock_populated': True,
        'work_counts': {k: reports[0]['convergence']['work'][k] for k in work_keys},
        'results_sha256': dict(zip(['off', 'timed', 'counts'], [sha(p / 'results.json') for p in paths])),
        'calls_by_bucket': {b: sum(r['calls'] for r in counts if r['bucket'] == b)
                            for b in sorted({r['bucket'] for r in counts})}}
    return summary


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('usage: verify_cli.py OFF_OUTPUT_DIR TIMED_OUTPUT_DIR COUNTS_OUTPUT_DIR')
    print(json.dumps(verify([Path(p) for p in sys.argv[1:]]), indent=2))
