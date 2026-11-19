#!/usr/bin/env python3
"""Summarize F-LADDER-FRACTION raw measurements, retaining failures and units."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
BASE = '0.5:1e-4:default'
LOW = '0.4:1e-4:default'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not a.out.resolve().is_relative_to((ROOT / 'target/scratch').resolve()):
        p.error('--out must be under target/scratch')
    a.out.mkdir(parents=True, exist_ok=True)
    rows = json.loads(a.rows.read_text())
    failures = json.loads((a.rows.parent / 'failures.json').read_text())
    assert len({(r['case'], r['seed'], r['arm']) for r in rows}) == len(rows)
    for r in rows:
        # Independent Parseval decomposition check and diagnostic range checks.
        assert math.isclose(r['handoff_density_error']**2,
                            r['density_low_pw_error']**2 + r['density_high_pw_error']**2,
                            rel_tol=1e-10, abs_tol=1e-25)
        assert 0 <= r['density_outside_coarse_density_error'] <= r['handoff_density_error']
        assert len(r['principal_angles_rad']) == r['occupied']
        assert all(0 <= angle <= math.pi / 2 for angle in r['principal_angles_rad'])
    baseline = {(r['case'], r['seed']): r for r in rows if r['arm'] == BASE}
    low = {(r['case'], r['seed']): r for r in rows if r['arm'] == LOW}
    cases = sorted({r['case'] for r in rows})
    arms = list(dict.fromkeys(r['arm'] for r in rows))
    summary = []
    for case in cases:
        for arm in arms:
            group = [r for r in rows if (r['case'], r['arm']) == (case, arm)]
            if not group:
                continue
            assert {r['seed'] for r in group} == {r['seed'] for r in rows if r['case'] == case and r['arm'] == BASE}
            ratio = lambda r, key: r[key] / baseline[case, r['seed']][key]
            fine = lambda r: r['fine_apply_h'] + r['verification_apply_h']
            med = lambda key: statistics.median(r[key] for r in group)
            item = dict(case=case, arm=arm, seeds=len(group), fraction=med('effective_fraction'),
                        coarse_maps=med('coarse_maps'), fine_maps=med('fine_maps'),
                        verification_maps=med('verification_maps'),
                        total_h=med('total_h'), total_h_eq=med('total_h_eq'),
                        total_h_ratio=statistics.median(ratio(r, 'total_h') for r in group),
                        total_h_eq_ratio=statistics.median(ratio(r, 'total_h_eq') for r in group),
                        fft_ratio=statistics.median(ratio(r, 'fft_column_eq') for r in group),
                        fine_ratio=statistics.median(fine(r)/fine(baseline[case,r['seed']]) for r in group),
                        seeds_more_fine_h=sum(fine(r)>fine(baseline[case,r['seed']]) for r in group),
                        seeds_more_fine_maps=sum(r['fine_maps']>baseline[case,r['seed']]['fine_maps'] for r in group),
                        worst_error_ha=max(r['error_ha'] for r in group))
            summary.append(item)
    watch_checks = 0
    for r in rows:
        if r['arm'].endswith(':watch'):
            b = low[r['case'],r['seed']]
            for key in ['energy_ha','coarse_maps','fine_maps','verification_maps','coarse_apply_h','fine_apply_h','verification_apply_h','coarse_residual']:
                assert r[key] == b[key], (r['case'],r['seed'],key)
            # Byte identity of the actual handoff arrays; both full ladders also
            # compare delivered density against their production-driver controls.
            job = a.rows.parent / f"{r['case']}-{r['seed']}"
            by_arm = {json.loads(f.read_text())['arm']: f.parent for f in job.glob('arm-*/measurement.json')}
            for name in ['handoff-density.bin','handoff-orbitals.bin','fine-history.json','coarse-history.json']:
                assert (by_arm[r['arm']]/name).read_bytes() == (by_arm[LOW]/name).read_bytes(), name
            watch_checks += 1
    verdicts = []
    for arm in arms:
        group = [s for s in summary if s['arm'] == arm]
        wins = sum(s['total_h_eq_ratio'] < 1 for s in group)
        raw_wins = sum(s['total_h_ratio'] < 1 for s in group)
        regressions = sum(s['seeds_more_fine_h'] > 0 or s['seeds_more_fine_maps'] > 0 for s in group)
        worst = max(s['worst_error_ha'] for s in group)
        verdicts.append(dict(arm=arm, winning_systems=wins, raw_h_winning_systems=raw_wins,
                             fine_regression_systems=regressions,
                             worst_error_ha=worst, passes=len(group)>=7 and wins>=5 and regressions==0 and worst<5e-9 and not failures))
    result = dict(rows=len(rows), failures=failures, rows_sha256=hashlib.sha256(a.rows.read_bytes()).hexdigest(),
                  watch_identity_checks=watch_checks, summary=summary, verdicts=verdicts,
                  worst=max(rows,key=lambda r:r['error_ha']))
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines = ['## Candidate matrix (Mac one core, three geometries per system)', '',
             'H includes coarse atomic start, fine redo and every verification. H-eq weights coarse H by the coarse/fine grid volume. Ratios are medians of paired ratios. Fine H includes verification; `increases` counts seeds with more fine H / more fine maps. W also charges delta-potential FFT columns. No wall claim is made.', '',
             '| System | Arm | Effective f | C/F/V maps | Total H | H-eq | H/base | H-eq/base | W/base | Fine H/base | Increases H/maps | Worst abs dE (Ha) |',
             '| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: |']
    for s in summary:
        arm=s['arm'].split(':')[-1] if s['arm'] not in (BASE,LOW) else s['arm'].split(':')[0]
        lines.append(f"| {s['case']} | {arm} | {s['fraction']:.7f} | {s['coarse_maps']:g}/{s['fine_maps']:g}/{s['verification_maps']:g} | {s['total_h']:g} | {s['total_h_eq']:.2f} | {s['total_h_ratio']:.4f} | {s['total_h_eq_ratio']:.4f} | {s['fft_ratio']:.4f} | {s['fine_ratio']:.4f} | {s['seeds_more_fine_h']}/{s['seeds_more_fine_maps']} | {s['worst_error_ha']:.4e} |")
    lines += ['', '### Per-geometry accounting', '', '| System | Seed | Arm | C/F/V maps | Hc/Hf/Hv | Total H | H-eq | abs dE (Ha) |', '| --- | ---: | --- | --- | --- | ---: | ---: | ---: |']
    for r in rows:
        lines.append(f"| {r['case']} | {r['seed']} | {r['arm']} | {r['coarse_maps']}/{r['fine_maps']}/{r['verification_maps']} | {r['coarse_apply_h']}/{r['fine_apply_h']}/{r['verification_apply_h']} | {r['total_h']} | {r['total_h_eq']:.3f} | {r['error_ha']:.6e} |")
    (a.out/'tables.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['rows','failures','watch_identity_checks','verdicts']},indent=2))


if __name__ == '__main__':
    main()
