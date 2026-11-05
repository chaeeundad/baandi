#!/usr/bin/env python3
"""Validate stored P14A accounting identities and regenerate the work tables."""
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent / 'reference/perf/p14a'


def read(name):
    return json.loads((ROOT / (name + '.json')).read_text())


def scientific(value):
    if isinstance(value, dict):
        return {k: scientific(v) for k, v in value.items()
                if k not in ('timings', 'phase_fractions')}
    if isinstance(value, list):
        return [scientific(v) for v in value]
    return value


def main():
    baseline = read('factor1')['cases']
    final = read('default-final')['cases']
    threaded = read('threads4')['cases']
    comparison = []
    for name, entry in final.items():
        assert scientific(entry['history']) == scientific(baseline[name]['history'])
        other = threaded[name]
        a, b = entry['measurements'], other['measurements']
        row = dict(case=name,
                   energy_identical=a['total_energy_ha'].hex() == b['total_energy_ha'].hex(),
                   density_identical=entry['density_sha256'] == other['density_sha256'],
                   iterations_identical=a['iterations'] == b['iterations'],
                   work_identical=a['operator_applications'] == b['operator_applications'])
        assert all(value for key, value in row.items() if key != 'case')
        comparison.append(row)
    (ROOT / 'thread-comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')
    lines = ['# P14A work accounting', '',
             'Vector applications count bands passed to H. Block calls count calls of H on an entire block; a block can shrink during Davidson. Steps count Rayleigh–Ritz iterations. None of these units is interchangeable.', '',
             'Timings are shared-host observations, not accepted speedups. The source, binary hashes, inputs and units are recorded in the JSON files.', '',
             '## Loosened schedule experiment', '',
             'Only the density-residual term is multiplied. Initial cap 1e-2 Ha and floor 1e-9 Ha remain unchanged. Final fixed-point certificate is unchanged. Steps exclude the separate fixed-point certificate; H totals include it.', '',
             '| Case | Factor | Converged | SCF count | Davidson steps | Total H vectors | Energy difference from factor 1 (Ha) | Same SCF count |',
             '| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |']
    for factor in (1, 3, 10, 30, 100):
        for name, entry in read('factor' + str(factor))['cases'].items():
            ref = baseline[name]['measurements']
            if 'measurements' not in entry:
                # These rejected runs retain complete Debug history in stderr.
                # Count only ordinary records; do not treat them as converged output.
                stderr = entry['stderr']
                iterations = re.findall(r'\biteration: (\d+)', stderr)
                steps = re.findall(r'\bdavidson_iterations: (\d+)', stderr)
                vectors = re.findall(r'\boperator_applications: (\d+)', stderr)
                assert len(iterations) == len(steps) == len(vectors)
                lines.append(f'| {name} | {factor} | diverged (exit {entry["exit_code"]}) | {len(iterations)} completed | {sum(map(int, steps))} | {sum(map(int, vectors))} | — | no |')
                continue
            m = entry['measurements']
            steps = sum(r['davidson_iterations'] for r in entry['history'])
            same = m['converged'] and m['iterations'] == ref['iterations']
            lines.append(f"| {name} | {factor} | {m['converged']} | {m['iterations']} | {steps} | {m['operator_applications']} | {m['total_energy_ha']-ref['total_energy_ha']:.6e} | {same} |")
    lines += ['', 'No tested factor above 1 preserves the SCF count on any of the four cases. This is a measured bound within this capped, fixed-factor family, not a bound on every possible adaptive schedule. Equal iteration count alone would not establish an identical trajectory.', '',
              '## H7: band-count experiments', '',
              '| Case | Extra bands | Converged | SCF count | Total H vectors | Energy difference (Ha) | Engine seconds |',
              '| --- | ---: | --- | ---: | ---: | ---: | ---: |']
    for file, extra in [('baseline', 4), ('extra1', 1), ('extra2', 2)]:
        for name, entry in read(file)['cases'].items():
            m = entry['measurements']
            difference = m['total_energy_ha'] - baseline[name]['measurements']['total_energy_ha']
            lines.append(f"| {name} | {extra} | {m['converged']} | {m['iterations']} | {m['operator_applications']} | {difference:.6e} | {m['engine_total_s']:.6f} |")
    for name, entry in read('accounting')['cases'].items():
        ref = baseline[name]
        assert entry['measurements']['total_energy_ha'] == ref['measurements']['total_energy_ha']
        assert entry['measurements']['operator_applications'] == ref['measurements']['operator_applications']
        lines += ['', '## Per-SCF baseline: ' + name, '',
                  'Early/late count band roots over all k points. Early means permanently below its threshold before the last Davidson iteration; late means first permanently below at the last iteration. Buffer fraction attributes corrections to their originating root before orthogonalization, not a physical decomposition of mixed vectors. Unknown lineage is reported separately.', '',
                  '| SCF | Tolerance (Ha residual) | Steps | H blocks | H vectors | H / band / k | Restarts | Early / late | Buffer-origin % | Unknown-origin H | Certificate H | Potential / diagonalization / density / energy / mixing (%) |',
                  '| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |']
        offset = 0
        for row in entry['history']:
            nk = len(row['operator_applications_per_k'])
            records = entry['davidson_work'][offset:offset+nk]
            offset += nk
            assert len(records) == nk and all(w['returned'] for w in records)
            assert [w['vector_applications'] for w in records] == row['operator_applications_per_k']
            assert [w['iterations'] for w in records] == row['davidson_iterations_per_k']
            assert all(w['tolerance_ha'] == row['requested_eigen_tolerance_ha'] for w in records)
            cert = row['fixed_point_certificate']
            if cert:
                checks = entry['davidson_work'][offset:offset+nk]
                offset += nk
                assert [w['vector_applications'] for w in checks] == cert['operator_applications_per_k']
            total = row['operator_applications']
            blocks = sum(w['block_calls'] for w in records)
            restarts = sum(w['restarts'] for w in records)
            early = sum(0 < first and last < w['iterations']-1 for w in records for first, last in zip(w['first_below'], w['last_above']))
            late = sum(0 < first and last == w['iterations']-1 for w in records for first, last in zip(w['first_below'], w['last_above']))
            buffer = sum(w['buffer_origin_applications'] for w in records)
            unknown = sum(w['unattributed_applications'] for w in records)
            phases = ' / '.join(f"{100*row['phase_fractions'][k]:.1f}" for k in ('potential_s','diagonalization_s','density_s','total_energy_s','mixing_s'))
            lines.append(f"| {row['iteration']} | {row['requested_eigen_tolerance_ha']:.4e} | {row['davidson_iterations']} | {blocks} | {total} | {total/entry['measurements']['bands']/nk:.4f} | {restarts} | {early} / {late} | {100*buffer/total:.2f} | {unknown} | {cert['operator_applications'] if cert else 0} | {phases} |")
        assert offset == len(entry['davidson_work'])
    lines += ['', '## si8 alignment with relayed QE data', '',
              'QE data are relayed, not rerun here. QE uses eigenvalue-change thresholds in Ry; our residual norm thresholds in Ha are different tests. Local checkout differs from relayed fdce016. The fifteenth local iteration must not be dropped to claim equal SCF convergence.', '',
              '| SCF | QE ethr (Ry) | QE steps | Relayed BAANDI steps | Local BAANDI steps | Local residual tolerance (Ha) |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    qe = [(1e-2,3,9),(1.36e-3,4,10),(6.03e-5,4,9),(2.55e-6,12,15),(2.55e-6,2,10),(4.35e-8,7,14),(4.35e-8,2,14),(2.45e-8,2,16),(5.82e-9,11,29),(1.37e-9,1,19),(3.36e-10,11,16),(3.36e-10,2,22),(9.94e-11,20,54),(5.94e-11,2,21)]
    for row, (threshold, steps, relayed) in zip(baseline['si8']['history'], qe):
        lines.append(f"| {row['iteration']} | {threshold:.3e} | {steps} | {relayed} | {row['davidson_iterations']} | {row['requested_eigen_tolerance_ha']:.4e} |")
    last = baseline['si8']['history'][14]
    lines.append(f"| 15 | — | — | 32 | {last['davidson_iterations']} | {last['requested_eigen_tolerance_ha']:.4e} |")
    lines += ['', 'QE totals: 83 steps and 98 h_psi block calls. QE block widths are not present in the relay, so vector-application totals cannot be reconstructed. The k444s relay is internally ambiguous about five versus six SCF iterations and omits one threshold/error pairing; it is not expanded into invented iteration rows.', '',
              '## Stored k222 / k222s QE outputs', '',
              '| Case | SCF | ethr (Ry) | Average iterations |', '| --- | ---: | ---: | ---: |']
    for row in read('baseline')['qe']:
        lines.append(f"| {Path(row['source']).parent.name} | {row['iteration']} | {row['ethr_ry']:.3e} | {row['average_iterations']:.1f} |")
    (ROOT / 'work-tables.md').write_text('\n'.join(lines) + '\n')
    print('P14A accounting identities and instrumented/default energy identity: PASS')


if __name__ == '__main__':
    main()
