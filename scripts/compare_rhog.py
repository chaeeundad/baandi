#!/usr/bin/env python3
"""Compare independent density DFT coefficients with QE; run with python3 -I."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

# Isolated Python deliberately excludes the script directory from sys.path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qe_rhog import (FormatError, finite, load_ours, miller, read_json, read_qe,
                     require, validate_grid)

DEFAULTS = dict(coefficient=1e-7, rho0=1e-12, parseval=1e-6, relative_l2=1e-6)
ZERO = (0, 0, 0)


def pair(z):
    return dict(re=z.real, im=z.imag)


def norm2(values):
    return math.fsum(abs(z)**2 for z in values)


def relative(difference, denominator):
    # JSON-safe representation of a nonzero difference against zero.
    return difference / denominator if denominator else (0.0 if difference == 0 else None)


def noise(ours, qe):
    group = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    available = [m for m in group if m in ours and m in qe]
    return dict(informational_only=True, definition='maximum pairwise complex distance',
                millers=available, complete_triplet=len(available) == 3,
                ours_spread=max((abs(ours[a]-ours[b]) for a in available for b in available), default=None),
                qe_spread=max((abs(qe[a]-qe[b]) for a in available for b in available), default=None))


def compare(ours, qe, tolerances, summary=None):
    """Compare maps; summary mode uses aggregates and only listed coefficients."""
    require(ZERO in ours and ZERO in qe, 'missing rho(0)')
    failures = []
    report = dict(mode='summary' if summary else 'file', units='e/Bohr^3',
                  sum_abs2_units='e^2/Bohr^6', tolerances=tolerances,
                  non_default_tolerances={k: v for k, v in tolerances.items() if v != DEFAULTS[k]})
    only_ours, only_qe = sorted(ours.keys()-qe.keys()), sorted(qe.keys()-ours.keys())
    if not summary:
        report['set_comparison'] = dict(only_ours_count=len(only_ours), only_qe_count=len(only_qe),
                                        only_ours=only_ours, only_qe=only_qe)
        if only_ours or only_qe:
            failures.append('G set mismatch')
    else:
        require(not only_qe, f'summary samples outside our density sphere: {only_qe}')
        require(summary['ngm_g'] == len(ours), 'summary/ours density G count mismatch')
        report['set_comparison'] = dict(verified=False, note=
            'Summary omits the full QE Miller list. Count and listed indices checked; full set equality cannot be certified.')
    shared = sorted(ours.keys() & qe.keys())
    samples = [dict(miller=m, ours=pair(ours[m]), qe=pair(qe[m]),
                    abs_difference=abs(ours[m]-qe[m])) for m in shared]
    worst = max(samples, key=lambda x: x['abs_difference'])
    report['max_coefficient_difference'] = worst
    report['samples' if summary else 'coefficients'] = samples
    for row in samples:
        if row['abs_difference'] > tolerances['coefficient']:
            failures.append(f"coefficient {tuple(row['miller'])}: {row['abs_difference']:.17g}")
    delta0 = abs(ours[ZERO] - qe[ZERO])
    report['rho0'] = dict(ours=pair(ours[ZERO]), qe=pair(qe[ZERO]), abs_difference=delta0)
    if delta0 > tolerances['rho0']:
        failures.append('rho(0) tolerance exceeded')
    ours_sum = norm2(ours.values())
    qe_sum = finite(summary['sum_abs2']) if summary else norm2(qe.values())
    require(qe_sum >= 0, 'negative sum_abs2')
    delta_sum = abs(ours_sum - qe_sum)
    rel_sum = relative(delta_sum, qe_sum)
    report['parseval'] = dict(ours_sum_abs2=ours_sum, qe_sum_abs2=qe_sum,
                              abs_difference=delta_sum, relative_difference=rel_sum,
                              denominator='QE sum_abs2')
    if rel_sum is None or rel_sum > tolerances['parseval']:
        failures.append('Parseval sum tolerance exceeded')
    if summary:
        ours_m = max((m for m in ours if m != ZERO), key=lambda m: abs(ours[m]))
        qe_m = miller(summary['max_abs_nonzero']['miller'])
        require(qe_m != ZERO and qe_m in ours, 'invalid maximum nonzero Miller')
        qe_max = finite(summary['max_abs_nonzero']['value'])
        require(qe_max >= 0, 'negative max_abs_nonzero')
        delta_max = abs(abs(ours[ours_m]) - qe_max)
        report['max_abs_nonzero'] = dict(ours_value=abs(ours[ours_m]), ours_miller=ours_m,
            qe_value=qe_max, qe_miller=qe_m, abs_difference=delta_max,
            note='Only the magnitude is gated; differing symmetry-equivalent maximum triplets do not fail.')
        if delta_max > tolerances['coefficient']:
            failures.append('maximum nonzero magnitude tolerance exceeded')
    else:
        squared = math.fsum(row['abs_difference']**2 for row in samples)
        rel_l2_sq = relative(squared, ours_sum)
        rel_l2 = math.sqrt(rel_l2_sq) if rel_l2_sq is not None else None
        report['rms_difference'] = math.sqrt(squared / len(samples))
        report['relative_l2_difference'] = rel_l2
        report['metrics_scope'] = 'Common G only; set mismatch is an inconsistent-input failure.'
        report['max_nonzero_difference'] = max((r for r in samples if r['miller'] != ZERO),
                                               key=lambda r: r['abs_difference'], default=None)
        if rel_l2 is None or rel_l2 > tolerances['relative_l2']:
            failures.append('relative L2 tolerance exceeded')
    report['symmetry_noise'] = noise(ours, qe)
    report['failures'] = failures
    report['exit_code'] = 3 if not summary and (only_ours or only_qe) else (5 if failures else 0)
    report['passed'] = report['exit_code'] == 0
    return report


def load_summary(path, grid):
    summary = read_json(path)
    require(summary['gamma_only'] == 0, 'unsupported gamma_only: expected 0')
    require(summary['nspin'] == 1, 'unsupported nspin: expected 1')
    require(type(summary['ngm_g']) is int and summary['ngm_g'] > 0, 'invalid ngm_g')
    require(isinstance(summary['samples'], list) and summary['samples'], 'summary sample list must be non-empty')
    samples = {}
    for row in summary['samples']:
        m = miller(row['miller'])
        require(m not in samples, 'duplicate Miller triplet in summary')
        samples[m] = complex(finite(row['re']), finite(row['im']))
    zero = complex(finite(summary['rho_000']['re']), finite(summary['rho_000']['im']))
    require(ZERO not in samples or samples[ZERO] == zero, 'inconsistent summary rho(0)')
    samples[ZERO] = zero
    validate_grid(grid, samples)
    validate_grid(grid, [miller(summary['max_abs_nonzero']['miller'])])
    return samples, summary


def markdown(report):
    lines = ['# Fourier density comparison', '',
             ('PASS' if report['passed'] else 'FAIL') + f" ({report['mode']} mode)", '',
             'Coefficient differences and symmetry spreads are in e/Bohr^3.', '',
             '| Metric | Measured difference | Tolerance |',
             '| --- | ---: | ---: |']
    for name, value, gate in [
            ('Largest coefficient', report['max_coefficient_difference']['abs_difference'], 'coefficient'),
            ('rho(0)', report['rho0']['abs_difference'], 'rho0'),
            ('Relative sum of squares', report['parseval']['relative_difference'], 'parseval')]:
        lines.append(f"| {name} | {value} | {report['tolerances'][gate]} |")
    if report['mode'] == 'file':
        lines.append(f"| Relative L2 | {report['relative_l2_difference']} | {report['tolerances']['relative_l2']} |")
    lines += ['', f"Non-default tolerances: {report['non_default_tolerances']}", '',
              f"Largest difference Miller triplet: {report['max_coefficient_difference']['miller']}", '',
              f"Set comparison: {report['set_comparison']}", '']
    if report['mode'] == 'summary':
        maximum = report['max_abs_nonzero']
        lines += [maximum['note'], '',
                  f"Maximum |rho|: ours {maximum['ours_value']} at {maximum['ours_miller']}; "
                  f"QE {maximum['qe_value']} at {maximum['qe_miller']}.", '']
    spread = report['symmetry_noise']
    lines += [f"Informational symmetry spread: ours {spread['ours_spread']}; QE {spread['qe_spread']}.", '',
              'Full measurements, coefficient values, validation and source metadata:', '',
              '```json', json.dumps(report, indent=2, allow_nan=False), '```', '']
    return '\n'.join(lines)


class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(1, f'usage error: {message}\n')


def main(argv=None):
    parser = Parser(description=__doc__)
    parser.add_argument('--ours-dir', required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--qe-summary', type=Path)
    source.add_argument('--qe-file', type=Path)
    parser.add_argument('--json-report', required=True, type=Path)
    parser.add_argument('--markdown-report', required=True, type=Path)
    for key, value in DEFAULTS.items():
        parser.add_argument('--' + key.replace('_', '-') + '-tol', type=float, default=value)
    args = parser.parse_args(argv)
    try:
        outputs = [args.json_report, args.markdown_report]
        require(outputs[0].resolve() != outputs[1].resolve(), 'report paths must differ')
        for path in outputs:
            if path.exists() or path.is_symlink():
                raise OSError(f'refusing to overwrite existing report/evidence: {path}')
        tolerances = {k: finite(getattr(args, k + '_tol')) for k in DEFAULTS}
        require(all(v >= 0 for v in tolerances.values()), 'tolerances must be nonnegative')
        ours, meta, validation = load_ours(args.ours_dir)
        summary = None
        source_path = args.qe_summary or args.qe_file
        if args.qe_summary:
            qe, summary = load_summary(source_path, meta['fft_grid'])
            metadata = summary
        else:
            qe, metadata = read_qe(source_path)
            validate_grid(meta['fft_grid'], qe)
        report = compare(ours, qe, tolerances, summary)
        report['validation'] = validation
        report['source'] = dict(path=str(source_path), sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                                metadata=metadata)
        report['ours_dir'] = str(args.ours_dir)
        report['override_flags_supplied'] = [arg for arg in (sys.argv[1:] if argv is None else argv)
                                              if arg.startswith('--') and '-tol' in arg]
        # Reserve both paths exclusively before writing either report.
        with outputs[0].open('x') as js:
            try:
                md = outputs[1].open('x')
            except OSError:
                outputs[0].unlink()
                raise
            with md:
                js.write(json.dumps(report, indent=2, allow_nan=False) + '\n')
                md.write(markdown(report))
        print(json.dumps({k: report[k] for k in ('passed', 'exit_code', 'tolerances',
              'non_default_tolerances', 'max_coefficient_difference', 'rho0', 'parseval', 'failures')}, indent=2))
        return report['exit_code']
    except OSError as exc:
        print(f'I/O error: {exc}', file=sys.stderr)
        return 1
    except (FormatError, KeyError, TypeError, ValueError, OverflowError) as exc:
        print(f'inconsistent inputs: {exc}', file=sys.stderr)
        return 3


if __name__ == '__main__':
    sys.exit(main())
