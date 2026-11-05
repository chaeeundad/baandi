#!/usr/bin/env python3
"""Validate GEMM identities and generate P14B tables from recorded runs."""
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent / 'reference/perf/p14b'


def read(name):
    return json.loads((ROOT / (name + '.json')).read_text())['cases']


def records(entry):
    for line in entry['stderr'].splitlines():
        if not line.startswith('P14B_GEMM,'):
            continue
        parts = line.split(',')[1:]
        if len(parts) == 9:
            solve, site, m, n, k, calls, flops, minimum, ns = parts
            solve, minimum = int(solve), int(minimum)
        else:
            site, m, n, k, calls, flops, ns = parts
            solve, minimum = None, None
        m, n, k, calls, flops, ns = map(int, (m, n, k, calls, flops, ns))
        assert flops == 2*m*n*k*calls
        assert minimum is None or 0 <= minimum <= flops
        yield dict(solve=solve, site=site, m=m, n=n, k=k,
                   calls=calls, flops=flops, minimum=minimum, ns=ns)


def main():
    baseline = read('accounting-final')
    for name, entry in baseline.items():
        counts = lambda e: [{k:v for k,v in r.items() if k != 'ns'} for r in records(e)]
        assert counts(entry) == counts(read('accounting-detailed')[name])
    trial = read('incremental-trial')
    threads = read('incremental-threads4')
    comparisons = []
    for name, e in read('default-final').items():
        b, t = read('default-before')[name], read('default-threads4')[name]
        assert e['density_sha256'] == b['density_sha256'] == t['density_sha256']
        for key in ('total_energy_ha', 'iterations', 'operator_applications'):
            assert e['measurements'][key] == b['measurements'][key] == t['measurements'][key]
        comparisons.append(dict(case=name, energy_density_work_thread_identity=True))
    for name, e in read('suite-final').items():
        b = read('suite-before')[name]
        for key in ('status', 'converged', 'total_energy_ha', 'iterations', 'density_sha256', 'operator_applications'):
            assert e.get(key) == b.get(key), (name, key)
        # Timing diagnostics were emitted on stderr before the unchanged error.
        for a_run, b_run in zip(e['runs'], b['runs']):
            strip = lambda text: '\n'.join(line for line in (text or '').splitlines()
                                           if not line.startswith('P14B_GEMM,'))
            assert strip(a_run.get('error')) == strip(b_run.get('error'))
        comparisons.append(dict(case=name, baseline_final_scientific_identity=True))
    (ROOT / 'identity.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    lines = ['# P14B logical GEMM accounting', '',
             'Flops use 2mnk; complex real-operation equivalents are four times these counts. Times are summed logical-call wall times in seconds (includes panel scheduling and assembly), not backend-only CPU time. Shared Mac observations are not acceptance timings.', '',
             'The minimum column reports a conservative bound for the **same unlocked residual-based step**. Reduced matrices need new rows/columns only. Other categories preserve current certification and two-pass projection; their bounds do not assume unvalidated locking or adaptive QR.', '']
    for name, entry in baseline.items():
        original = read('default-before')[name]
        assert entry['density_sha256'] == original['density_sha256']
        assert entry['measurements']['total_energy_ha'] == original['measurements']['total_energy_ha']
        phases = {}
        offset = 0
        for row in entry['history']:
            nk = len(row['operator_applications_per_k'])
            for j in range(nk):
                w = entry['davidson_work'][offset+j]
                assert w['iterations'] == row['davidson_iterations_per_k'][j]
                phases[offset+j] = 'ordinary'
            offset += nk
            if row['fixed_point_certificate']:
                for j in range(nk):
                    phases[offset+j] = 'certificate'
                offset += nk
        assert offset == len(entry['davidson_work'])
        groups = defaultdict(lambda: [0, 0, 0, 0])
        rs = list(records(entry))
        assert {r['solve'] for r in rs} == set(phases)
        for r in rs:
            g = groups[(phases[r['solve']], r['site'])]
            for i, key in enumerate(('calls', 'flops', 'minimum', 'ns')):
                g[i] += r[key]
        total = sum(r['ns'] for r in rs)
        lines += ['## '+name, '', '| Phase | Site | Calls | 2mnk | Minimum 2mnk | Seconds | GEMM time % |',
                  '| --- | --- | ---: | ---: | ---: | ---: | ---: |']
        for (phase, site), (calls, flops, minimum, ns) in sorted(groups.items()):
            lines.append(f'| {phase} | {site} | {calls} | {flops} | {minimum} | {ns/1e9:.6f} | {100*ns/total:.2f} |')
        assert not any(r['site'] == 'other' for r in rs)
        lines += ['', 'Restart contraction: 0 extra GEMMs (saved Ritz H products). Internal block MGS: scalar dot/AXPY, 0 GEMMs. There is no separate generalized-overlap reduced matrix for this orthonormal-basis solver. Projection overlaps are listed explicitly.', '']
    lines += ['## Initial incremental-Gram experiment (subsequently adopted)', '',
              '| Case | Before 2mnk | Trial 2mnk | SCF before / trial | Energy change (Ha) | T1/T4 density identical |',
              '| --- | ---: | ---: | --- | ---: | --- |']
    for name, entry in baseline.items():
        a, b = entry['measurements'], trial[name]['measurements']
        same = trial[name]['density_sha256'] == threads[name]['density_sha256']
        assert same and b['total_energy_ha'] == threads[name]['measurements']['total_energy_ha']
        before = sum(r['flops'] for r in records(entry))
        after = sum(r['flops'] for r in records(trial[name]))
        lines.append(f"| {name} | {before} | {after} | {a['iterations']} / {b['iterations']} | {b['total_energy_ha']-a['total_energy_ha']:.6e} | {same} |")
    before = read('suite-before')
    after = read('suite-trial')
    lines += ['', 'The four requested cases pass the narrow iteration/energy test, with these toy iteration changes (subsequently allowed by the architect):', '',
              '| Case | Before SCF | Trial SCF | Energy difference (Ha) |', '| --- | ---: | ---: | ---: |']
    for name, a in before.items():
        b = after[name]
        if a.get('converged') and a.get('iterations') != b.get('iterations'):
            lines.append(f"| {name} | {a['iterations']} | {b['iterations']} | {b['total_energy_ha']-a['total_energy_ha']:.6e} |")
    lines += ['', 'This table records the initial pre-adoption snapshots. The architect subsequently adopted the trial as default; see adopted/acceptance.md for the current rule and results. Exact shape histograms, including Ritz column widths, remain in each raw JSON stderr. Backend microkernel packing is not separately timed. No new QE run or cross-code performance claim is made.', '']
    (ROOT / 'gemm-tables.md').write_text('\n'.join(lines))
    print('GEMM shape/flop identities, solve attribution, instrumentation identity and trial thread identity: PASS')


if __name__ == '__main__':
    main()
