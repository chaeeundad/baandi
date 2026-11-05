"""Stress availability/schema, sign, term sum, and deliberate rejection probes."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
spec = importlib.util.spec_from_file_location('spin_schema', Path(__file__).with_name('validate_spin_results.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
validate, SCHEMAS = module.validate, module.SCHEMAS
schema = json.loads((SCHEMAS / 'results-state.proposal.json').read_text())
result = json.loads(Path(sys.argv[1]).read_text())
s = result['stress']
validate(s, schema['properties']['stress'])
if 'state_contracts' in result:
    validate(result['state_contracts'], schema)
for a in range(3):
    for b in range(3):
        assert abs(s['values'][a][b] - sum(t[a][b] for t in s['per_term'].values())) < 1e-12
        assert s['values'][a][b] == s['values'][b][a]
assert abs(s['pressure_ha_bohr3'] + sum(s['values'][a][a] for a in range(3))/3) < 1e-12
assert s['antisymmetric_norm_ha_bohr3'] < 1e-10
for label, mutate in [
    ('units', lambda x: x.update(units='GPa')),
    ('unavailable numeric', lambda x: x.update(status='unavailable', reason='missing')),
    ('incomplete terms', lambda x: x['per_term'].pop('core')),
    ('standard', lambda x: x['preconditions'].update(convergence_profile='standard')),
    ('certificate', lambda x: x['preconditions'].update(fixed_point_certificate_passed=False)),
    ('tensor shape', lambda x: x.update(values=[[0, 0]])),
    ('cutoff claim', lambda x: x.update(cutoff_converged=True)),
]:
    invalid = copy.deepcopy(s)
    mutate(invalid)
    try:
        validate(invalid, schema['properties']['stress'])
    except AssertionError:
        continue
    raise AssertionError(f'accepted {label}')
base = json.loads((Path(__file__).parent / 'fixtures/force_schema/accept-v2-unavailable.json').read_text())
base['stress'] = s
base['stress']['quantity'] = 'F' if base['occupations']['policy']=='fermi_dirac' else 'E'
validate(base, schema)
base['extension_version'] = 1
try:
    validate(base, schema)
except AssertionError:
    pass
else:
    raise AssertionError('version 1 accepted stress')
print('PASS stress schema and physical reporting semantics')
