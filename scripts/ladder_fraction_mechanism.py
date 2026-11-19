#!/usr/bin/env python3
"""Generate fraction mechanism tables from the ignored Rust experiment."""
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not a.out.resolve().is_relative_to((ROOT/'target/scratch').resolve()):
        p.error('--out must be under target/scratch')
    a.out.mkdir(parents=True, exist_ok=True)
    rows = json.loads((a.root/'rows.json').read_text())
    assert not json.loads((a.root/'failures.json').read_text())
    paths = {}
    for path in a.root.glob('*/arm-*/measurement.json'):
        r = json.loads(path.read_text())
        paths[path.parent.parent.name, r['arm']] = path.parent
    lines = ['## Fraction mechanism screen (seed 601)', '',
             'All density norms are Bohr^-3/2, energies and tolerances Ha, angles radians. Low G is the coarse orbital sphere; outside-density means beyond the 4Ec coarse density sphere. Raw histories contain every Pulay coefficient, conditioning and predicted norm.', '',
             '| System | f | C/F/V maps | Coarse residual | Coarse last standard error (Ha) | First fine tolerance (Ha) | H-eq | abs dE (Ha) |',
             '| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: |']
    histories = {}
    for r in rows:
        key=(f"{r['case']}-{r['seed']}",r['arm'])
        history=json.loads((paths[key]/'fine-history.json').read_text())
        coarse=json.loads((paths[key]/'coarse-history.json').read_text())
        assert sum(h['h'] for h in history)==r['fine_apply_h']
        assert sum(h['h'] for h in coarse)==r['coarse_h']
        assert len(history)==r['fine_maps'] and len(coarse)==r['coarse_maps']
        histories[r['case'],r['arm']]=history
        lines.append(f"| {r['case']} | {r['fraction']:g} | {r['coarse_maps']}/{r['fine_maps']}/{r['verification_maps']} | {r['coarse_residual']:.5e} | {coarse[-1]['standard_error_ha']:.5e} | {history[0]['eigen_tolerance_ha']:.5e} | {r['total_h_eq']:.2f} | {r['error_ha']:.3e} |")
    lines += ['', '### Handoff state compared with the strict solution', '',
              '| System | f | Density error | Low-G error | Outside-density error | Orbital energy error (Ha) | First-map energy error (Ha) | Max angle (rad) |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        h=histories[r['case'],r['arm']]
        assert math.isclose(r['handoff_density_error']**2,r['density_low_pw_error']**2+r['density_high_pw_error']**2,rel_tol=1e-10)
        lines.append(f"| {r['case']} | {r['fraction']:g} | {r['handoff_density_error']:.5e} | {r['density_low_pw_error']:.5e} | {r['density_outside_coarse_density_error']:.5e} | {r['handoff_orbital_energy_error_ha']:.5e} | {h[0]['energy_error_ha']:.5e} | {r['principal_angle_max_rad']:.5e} |")
    lines += ['', '### Handoff potential response', '',
              'The high-PW columns restore the high-PW part of the density error to its reference value and measure the actual nonlinear potential change. They are not bounds and are unavailable if that diagnostic density is negative.', '',
              '| System | f | Total Hartree RMS (Ha) | High-PW Hartree RMS (Ha) | Total XC RMS (Ha) | High-PW XC RMS (Ha) |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        fmt=lambda k:f"{r[k]:.5e}" if k in r else 'unavailable'
        lines.append(f"| {r['case']} | {r['fraction']:g} | {fmt('handoff_hartree_rms_ha')} | {fmt('high_pw_hartree_response_rms_ha')} | {fmt('handoff_xc_rms_ha')} | {fmt('high_pw_xc_response_rms_ha')} |")
    lines += ['', '### Fine-map trajectories: si32 and si64', '',
              'Reset is the cumulative feasibility-reset count; all arms start with empty fine history. Cpred is the predicted Coulomb norm ratio of the current L2 Pulay proposal, measured with the inert trace switch. `none` means no mixed proposal or no difference yet. All other systems have the same fields in their raw `fine-history.json` files.', '',
              '| System | f | Map | Standard error (Ha) | Eigen tolerance (Ha) | H | Davidson | Density residual | Resets | Pulay differences | Cpred |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in rows:
        if r['case'] not in ['si32','si64']:
            continue
        for h in histories[r['case'],r['arm']]:
            pulay=h['pulay'] or {}
            c=pulay.get('coulomb_predicted_ratio')
            c='none' if c is None else f'{c:.5g}'
            lines.append(f"| {r['case']} | {r['fraction']:g} | {h['map']} | {h['standard_error_ha']:.5e} | {h['eigen_tolerance_ha']:.5e} | {h['h']} | {h['davidson']} | {h['density_residual']:.5e} | {h['resets']} | {pulay.get('differences','none')} | {c} |")
    (a.out/'mechanism-tables.md').write_text('\n'.join(lines)+'\n')
    print(f'Validated {len(rows)} arms, {sum(len(h) for h in histories.values())} fine-map records')


if __name__ == '__main__':
    main()
