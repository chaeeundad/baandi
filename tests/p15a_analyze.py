#!/usr/bin/env python3
"""Offline stopping-rule replay. No fitted residual-to-eigenvalue conversion.
A candidate stops at a prefix of an actually executed deterministic trajectory.
Unreached criteria are censored, never assigned the final recorded step.
"""
import json
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'tests/reference/perf/p15a'
QE_SI8_ETHR_RY = [1e-2,1.36e-3,6.03e-5,2.55e-6,2.55e-6,4.35e-8,4.35e-8,2.45e-8,5.82e-9,1.37e-9,3.36e-10,3.36e-10,9.94e-11,5.94e-11]
LADDER = [10.0**(-i) for i in range(2, 14)]
RATES = [10.0**(-i) for i in range(2, 7)]
ALLOW = [0, 1, 2, 5]


def schedule(metrics):
    previous = 0.005  # Ha: 1e-2 Ry
    rows = []
    for m in metrics:
        iteration = m['iteration']
        redo = iteration == 1 and m['dr2_ry'] < 2 * previous * max(1, m['electrons'])
        after = max(5e-14, min(previous, 0.05 * m['dr2_ry'] / max(1, m['electrons'])))
        rows.append(dict(iteration=iteration, ethr_ha=previous, next_ethr_ha=after,
                         first_redo_trigger=redo, redo_ethr_ha=after if redo else None,
                         dr2_ry=m['dr2_ry']))
        previous = after
    return rows


def changes(trace, index):
    if index == 0:
        return None  # no invented zero previous Ritz value
    return [abs(a-b) for a, b in zip(trace['steps'][index]['eigenvalues_ha'],
                                     trace['steps'][index-1]['eigenvalues_ha'])]


def stop(trace, kind, tolerance, allowed=0, cap_only=False, occupied=None, limits_override=None):
    for i, row in enumerate(trace['steps']):
        values = changes(trace, i) if kind == 'ritz' else row['residuals_ha']
        if values is None:
            continue
        limits = limits_override or [tolerance if occupied is None or b < occupied else max(5*tolerance, 5e-6)
                                     for b in range(len(values))]
        # QE uses a strict < for Ritz changes, BAANDI residual uses <=.
        nc = sum(v >= t if kind == 'ritz' else v > t for v, t in zip(values, limits))
        if nc == 0 or (nc <= allowed and (not cap_only or row['step'] == 21)):
            return dict(step=row['step'], expansion_steps=row['step']-1,
                        h_vectors=row['h_vectors'], h_blocks=row['h_blocks'],
                        unconverged=nc, residuals_ha=row['residuals_ha'],
                        eigenvalues_ha=row['eigenvalues_ha'])
    return None


def error(hit, reference, occupied):
    if hit is None:
        return None
    return max(abs(a-b) for a, b in zip(hit['eigenvalues_ha'][:occupied], reference[:occupied]))


def rates(trace):
    result = {}
    for t in RATES:
        roots = []
        for b in range(len(trace['steps'][0]['residuals_ha'])):
            below = [r['step'] for r in trace['steps'] if r['residuals_ha'][b] <= t]
            above = [r['step'] for r in trace['steps'] if r['residuals_ha'][b] > t]
            last_above = max(above, default=0)
            roots.append(dict(first=below[0] if below else None,
                              sustained=last_above+1 if last_above < len(trace['steps']) else None,
                              recrossed=bool(below and last_above > below[0])))
        result[str(t)] = roots
    return result


def groups(data):
    offset = 0
    for h in data['report']['convergence']['history']:
        nk = len(h['davidson_iterations_per_k'])
        group = data['solves'][offset:offset+nk]
        assert len(group) == nk
        assert [len(s['0']['steps']) for s in group] == h['davidson_iterations_per_k']
        yield h, offset, group
        offset += nk
        cert = h.get('fixed_point_certificate')
        if cert:
            assert [len(s['0']['steps']) for s in data['solves'][offset:offset+nk]] == cert['davidson_iterations_per_k']
            offset += nk
    assert offset == len(data['solves']), 'unmapped frozen solves'


def main():
    result = dict(units='Ha for residuals/eigenvalues; Ry for dr2; work is vectors or explicitly blocks',
                  ladder_ha=LADDER, allowed_roots=ALLOW, cases={})
    table = ['# P15A frozen-H replay measurements', '',
             'Current native build; ordinary SCF solves only. Candidate counts are prefixes of the baseline trajectory.',
             'QE-like columns use the ratcheted schedule on the **BAANDI** history and occupied/empty thresholds.',
             'Censored solves are excluded from candidate sums and explicitly counted. No new SCF policy was run.', '',
             '| Case | SCF status / iterations | Native steps / H vectors | QE-like steps / H vectors (0 allowed) | Censored | E_res < 5e-9 Ha first iteration |',
             '| --- | --- | ---: | ---: | ---: | ---: |']
    for path in sorted(OUT.glob('*.json')):
        data = json.loads(path.read_text())
        if 'trace_file' in data:
            data['solves'] = json.loads(gzip.decompress((OUT / data['trace_file']).read_bytes()))
        if 'solves' not in data:
            continue
        history = data['report']['convergence']['history']
        occupied = int(data['report']['electrons']/2)
        sched = schedule(data['metrics'])
        assert len(sched) == len(history)
        snapshots = sorted(set([1, (len(history)+1)//2, len(history)]))
        analysis = dict(status='converged' if data['exit_code'] == 0 else 'exhausted',
                        iterations=len(history), snapshots=snapshots, schedule=sched,
                        per_iteration=[], snapshots_detail=[], rate_curves=[],
                        equivalence=data['equivalence'])
        e1 = next((r['iteration'] for r in sched if r['dr2_ry'] < 1e-8), None)
        analysis['first_energy_metric_iteration'] = e1
        if e1:
            analysis['energy_error_at_metric_ha'] = abs(history[e1-1]['energies']['total']['ha'] - data['report']['total_energy']['ha'])
        totals = dict(native_steps=0, native_h_vectors=0, ritz_steps=0, ritz_h_vectors=0, censored=0)
        for h, offset, group in groups(data):
            iteration = h['iteration']
            epsilon = sched[iteration-1]['ethr_ha']
            per = dict(iteration=iteration, ethr_ha=epsilon, native_steps=h['davidson_iterations'],
                       native_h_vectors=h['operator_applications'], kpoints=[])
            totals['native_steps'] += h['davidson_iterations']
            totals['native_h_vectors'] += h['operator_applications']
            for k, solve in enumerate(group):
                native, reference = solve['0'], solve['6']
                assert reference['success'], 'reference did not converge'
                ref = solve.get('dense_eigenvalues_ha', reference['steps'][-1]['eigenvalues_ha'])
                ref_delta = max(abs(a-b) for a,b in zip(ref,reference['steps'][-1]['eigenvalues_ha']))
                hit = stop(native, 'ritz', epsilon, occupied=occupied)
                if hit:
                    totals['ritz_steps'] += hit['step']
                    totals['ritz_h_vectors'] += hit['h_vectors']
                else:
                    totals['censored'] += 1
                kr = dict(dense_available='dense_eigenvalues_ha' in solve, tight_vs_dense_max_error_ha=ref_delta if 'dense_eigenvalues_ha' in solve else None, k=k, solve_id=offset+k, native_final_error_ha=error(native['steps'][-1], ref, occupied),
                          scheduled_ritz=hit, scheduled_ritz_error_ha=error(hit, ref, occupied),
                          same_numeric_tol_ritz=stop(native, 'ritz', native['tolerance_ha'], limits_override=[
                              native['tolerance_ha'] if b < native['strict_bands'] else native['buffer_tolerance_ha']
                              for b in range(data['report']['bands'])]),
                          qe_empty_native_tol_ritz=stop(native, 'ritz', native['tolerance_ha'], occupied=occupied),
                          relayed_qe_threshold_ritz=(stop(native, 'ritz', QE_SI8_ETHR_RY[iteration-1]/2, occupied=occupied)
                              if data['case'] == 'si8' and iteration <= 14 else None),
                          reference_scheduled_ritz=stop(reference, 'ritz', epsilon, occupied=occupied),
                          allowances=[])
                for allowed in ALLOW:
                    early = stop(native, 'ritz', epsilon, allowed, occupied=occupied)
                    cap = stop(native, 'ritz', epsilon, allowed, cap_only=True, occupied=occupied)
                    kr['allowances'].append(dict(allowed=allowed, early=early, qe_cap=cap,
                                                early_error_ha=error(early, ref, occupied)))
                per['kpoints'].append(kr)
                analysis['rate_curves'].append(dict(iteration=iteration, k=k, solve_id=offset+k,
                                                    all_root_reference=rates(reference)))
                if iteration in snapshots:
                    entry = dict(iteration=iteration, k=k, solve_id=offset+k, variants=[])
                    for key, trace in solve.items():
                        if not key.isdigit():
                            continue
                        candidates = []
                        for kind in ('residual', 'ritz'):
                            for tol in LADDER:
                                for allowed in ALLOW:
                                    candidate = stop(trace, kind, tol, allowed)
                                    candidates.append(dict(kind=kind, tolerance_ha=tol, allowed=allowed,
                                                           stop=candidate, occupied_error_ha=error(candidate, ref, occupied)))
                        # Equal delivered accuracy, versus the same strict reference.
                        # Select only stops actually meeting the target, not just a threshold label.
                        matched = []
                        for target in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
                            by_kind = {}
                            for kind in ('residual', 'ritz'):
                                valid = [c for c in candidates if c['kind'] == kind and c['stop']
                                         and c['occupied_error_ha'] <= target]
                                best = min(valid, key=lambda c: (c['stop']['h_vectors'], c['stop']['step'])) if valid else None
                                by_kind[kind] = best
                            matched.append(dict(target_ha=target, **by_kind))
                        entry['variants'].append(dict(variant=int(key), tolerance_ha=trace['tolerance_ha'],
                                                      rates=rates(trace), candidates=candidates, matched_accuracy=matched))
                    analysis['snapshots_detail'].append(entry)
            analysis['per_iteration'].append(per)
        analysis['totals'] = totals
        result['cases'][data['case']] = analysis
        table.append(f"| {data['case']} | {analysis['status']} / {len(history)} | {totals['native_steps']} / {totals['native_h_vectors']} | {totals['ritz_steps']} / {totals['ritz_h_vectors']} | {totals['censored']} | {e1 or 'not reached'} |")
    table += ['', '## Dense spectral-completeness check', '',
              'Counts below cover ordinary solves. A mismatch is max sorted eigenvalue difference > 1e-8 Ha.',
              'Dense spectra pass an independent true-residual check <= 1e-10 Ha.', '',
              '| Case | Dense checks | Tight replay mismatches | Maximum tight vs dense difference (Ha) |',
              '| --- | ---: | ---: | ---: |']
    for name, case in result['cases'].items():
        rows = [k for iteration in case['per_iteration'] for k in iteration['kpoints'] if k['dense_available']]
        differences = [r['tight_vs_dense_max_error_ha'] for r in rows]
        table.append(f"| {name} | {len(rows)} | {sum(d > 1e-8 for d in differences)} | {max(differences, default=0):.9g} |")
    table += ['', '## Equal occupied-eigenvalue accuracy on selected k=0 snapshots', '',
              'Best scanned prefix on the native trajectory, selected in hindsight at error <= 1e-10 Ha',
              'against a dense oracle when available, otherwise the tight all-root replay. Counts include the initial reduced solve.',
              'A dash means no scanned stop delivered that accuracy; this includes missing low roots, not just a short trace.', '',
              '| Case | SCF iteration | Residual steps / H vectors / allowed | Ritz steps / H vectors / allowed |',
              '| --- | ---: | ---: | ---: |']
    for name, case in result['cases'].items():
        for snapshot in case['snapshots_detail']:
            if snapshot['k'] != 0:
                continue
            match = next(m for m in snapshot['variants'][0]['matched_accuracy'] if m['target_ha'] == 1e-10)
            def display(c):
                return f"{c['stop']['step']} / {c['stop']['h_vectors']} / {c['allowed']}" if c else '-'
            table.append(f"| {name} | {snapshot['iteration']} | {display(match['residual'])} | {display(match['ritz'])} |")
    table += ['', '## Per-root rate summary (k=0, tight all-root trajectory)', '',
              'Each entry is the range across roots of the first / sustained step below the residual threshold.',
              'Sustained means no subsequent crossing in the recorded trace, not a locking guarantee.', '',
              '| Case | SCF iteration | 1e-2 Ha | 1e-3 Ha | 1e-4 Ha | 1e-5 Ha | 1e-6 Ha |',
              '| --- | ---: | --- | --- | --- | --- | --- |']
    for name, case in result['cases'].items():
        for curve in case['rate_curves']:
            if curve['k'] != 0 or curve['iteration'] not in case['snapshots']:
                continue
            cells = []
            for roots in curve['all_root_reference'].values():
                first = [r['first'] for r in roots if r['first'] is not None]
                sustained = [r['sustained'] for r in roots if r['sustained'] is not None]
                cells.append(f"{min(first)}-{max(first)} / {min(sustained)}-{max(sustained)}")
            table.append(f"| {name} | {curve['iteration']} | " + ' | '.join(cells) + ' |')
    # Root arrays are in raw traces; candidate records reference their exact step.
    def compact(node):
        if isinstance(node, dict):
            if 'eigenvalues_ha' in node and 'step' in node:
                node.pop('eigenvalues_ha')
                node['max_residual_ha'] = max(node.pop('residuals_ha'))
            for value in node.values():
                compact(value)
        elif isinstance(node, list):
            for value in node:
                compact(value)
    compact(result)
    (OUT / 'analysis.json').write_text(json.dumps(result, separators=(',', ':'), allow_nan=False)+'\n')
    table += ['', 'Full per-root rate curves, scan hits, residuals at each stop, matched occupied-eigenvalue errors,',
              'and SCF schedules are in `analysis.json`. Raw per-step values, Ritz changes, residuals,',
              'unconverged-root counts and H vector/block/restart counts are in each `*-traces.json.gz`.',
              'Reference: full-space dense spectrum where available; otherwise a tight replay (1e-10 Ha true residual).',
              'Some tight replays miss low roots: residual convergence alone is not spectral completeness. Eigenvalue',
              'differences at 1e-12 Ha are floating-point comparisons, not a certified absolute error bound.',
              'A 5-root early allowance is a separate diagnostic; the `qe_cap` column applies it only',
              'after 20 expansion steps (BAANDI trace row 21). No repeated-call QE implementation is inferred.',
              'First-iteration redo flags use the stored native residual. A redo changes the residual',
              'and starting vectors; its actual work and subsequent SCF trajectory are not inferred.', '']
    (OUT / 'tables.md').write_text('\n'.join(table))
    print('\n'.join(table[:16]))


if __name__ == '__main__':
    main()
