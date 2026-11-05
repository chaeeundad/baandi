"""Publish measured rebuild/polish results and verify the recorded acceptance."""
import argparse
import json
from pathlib import Path
import shutil


def publish(run_dir, output):
    data = json.loads((run_dir / 'summary.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    lines = ['# Rebuild and fresh-basis polish validation', '',
             'Lattice constants and differences are in Bohr; pressures in Ha/Bohr^3; wall times in seconds.', '',
             '| Case | Main a | Main rebuilt P | Final a | Final rebuilt P | Polish cycles | Total time |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for case in data['cases']:
        stages = case['stages']
        assert 3 <= len(stages) <= 4, 'require two or three completed polish cycles'
        first, last = stages[0], stages[-1]
        for i, stage in enumerate(stages):
            assert stage['exit_code'] == 0 and stage['status'] == 'converged'
            result_path = run_dir / stage['stage'] / 'results.json'
            result = json.loads(result_path.read_text())
            relaxation = result['relaxation']
            settings = relaxation['settings']
            assert settings['stress_tol_ha_bohr3'] == 1e-7
            assert settings['energy_tol_ha'] == 1e-9
            assert settings['rebuild_every'] == (case['rebuild_every'] if i == 0 else 0)
            check = relaxation['rebuilt_pressure_check']
            assert check['warning'] == (abs(stage['rebuilt_pressure_ha_bohr3']) >= 1e-7)
            assert relaxation['frames'][0]['provenance'] == 'fresh'
            for extension in ('toml', 'stdout', 'stderr'):
                source = run_dir / (stage['stage'] + '.' + extension)
                shutil.copyfile(source, output / source.name)
            shutil.copyfile(result_path, output / (stage['stage'] + '.results.json'))
        assert abs(last['rebuilt_pressure_ha_bohr3']) < 1e-7
        total = sum(stage['elapsed_s'] for stage in stages)
        case['total_elapsed_s'] = total
        case['polish_cycles'] = len(stages) - 1
        lines.append('| %s | %.12f | %.9e | %.12f | %.9e | %d | %.3f |' %
                     (case['case'], first['lattice_bohr'], first['rebuilt_pressure_ha_bohr3'],
                      last['lattice_bohr'], last['rebuilt_pressure_ha_bohr3'], len(stages)-1, total))
    lines += ['', '## Start comparisons', '',
              'Signed differences are expanded minus compressed. The tolerance estimates are unchanged from the frozen-only study.', '',
              '| Material | N | Main difference | After polish difference | Unchanged tolerance | Pass |',
              '| --- | ---: | ---: | ---: | ---: | --- |']
    for pair in data['comparisons']:
        assert pair['passed'] and pair['pressure_agreement'] and pair['start_agreement']
        lines.append('| %s %d Ha | %d | %.9e | %.9e | %.9e | yes |' %
                     (pair['material'], pair['cutoff_ha'], pair['rebuild_every'],
                      pair['before_polish_difference_bohr'], pair['final_difference_bohr'],
                      pair['unchanged_estimated_tolerance_bohr']))
    assert len(data['comparisons']) == len(data['cases']) // 2
    lines += ['', '## Each actual stage', '',
              '| Stage | Lattice | Frozen P | Rebuilt P | Time | Moves | Intermediate rebuilds |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for case in data['cases']:
        for stage in case['stages']:
            lines.append('| %s | %.12f | %.9e | %.9e | %.3f | %d | %d |' %
                         (stage['stage'], stage['lattice_bohr'], stage['frozen_pressure_ha_bohr3'],
                          stage['rebuilt_pressure_ha_bohr3'], stage['elapsed_s'], stage['steps'],
                          stage['basis_rebuild_count']))
    lines += ['', 'Each polish is a fresh optimizer calculation, not an exact electronic restart. '
              'Original frozen-only failures and cutoff dependence remain in real-followup/. '
              'This is finite-cutoff validation on the specified meshes, without comparison to experiment or other codes. '
              'Times are shared-Mac observations; the frozen-only rerun ran concurrently for part of this sequence.', '']
    (output / 'validation.json').write_text(json.dumps(data, indent=2) + '\n')
    (output / 'validation.md').write_text('\n'.join(lines))
    print('PASS: %d cases, %d stage executions, %d start pairs' %
          (len(data['cases']), sum(len(c['stages']) for c in data['cases']), len(data['comparisons'])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    publish(args.run_dir, args.output)
