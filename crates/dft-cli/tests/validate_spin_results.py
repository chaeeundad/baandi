"""Independent checks of the frozen proposal schema subset and standard cubes.

No third-party dependency. Unsupported schema keywords fail closed; this is
only the subset used by the three MF02 result proposal documents.
"""
import importlib.util
import json
import math
from pathlib import Path
import re
import sys
ROOT = Path(__file__).resolve().parents[3]
SCHEMAS = ROOT / 'docs/schemas'

def validate(v, s):
    allowed = {'$schema', '$comment', 'title', '$ref', 'type', 'const', 'enum',
               'required', 'properties', 'additionalProperties', 'items', 'prefixItems',
               'minItems', 'maxItems', 'minimum', 'maximum', 'exclusiveMinimum',
               'pattern', 'oneOf', 'allOf', 'if', 'then'}
    assert not set(s) - allowed, set(s) - allowed
    if '$ref' in s:
        validate(v, json.loads((SCHEMAS / s['$ref']).read_text()))
    if 'oneOf' in s:
        successes = 0
        for sub in s['oneOf']:
            try:
                validate(v, sub)
                successes += 1
            except AssertionError:
                pass
        assert successes == 1, (successes, v)
    for sub in s.get('allOf', []):
        validate(v, sub)
    if 'if' in s:
        try:
            validate(v, s['if'])
        except AssertionError:
            pass
        else:
            validate(v, s.get('then', {}))
    types = s.get('type', [])
    if isinstance(types, str):
        types = [types]
    checks = {'null': v is None, 'number': type(v) in (int, float),
              'integer': type(v) is int, 'array': isinstance(v, list),
              'object': isinstance(v, dict), 'string': isinstance(v, str), 'boolean': type(v) is bool}
    assert not types or any(checks[t] for t in types)
    if 'const' in s:
        assert v == s['const'] and (isinstance(s['const'], bool) == isinstance(v, bool)), (v, s['const'])
    if 'enum' in s:
        assert v in s['enum']
    if type(v) in (int, float):
        assert math.isfinite(v)
        assert v >= s.get('minimum', -math.inf) and v <= s.get('maximum', math.inf)
        assert v > s.get('exclusiveMinimum', -math.inf)
    if isinstance(v, str) and 'pattern' in s:
        assert re.search(s['pattern'], v)
    if isinstance(v, list):
        assert s.get('minItems', 0) <= len(v) <= s.get('maxItems', math.inf)
        prefix = s.get('prefixItems', [])
        for i, x in enumerate(v):
            validate(x, prefix[i] if i < len(prefix) else s.get('items', {}))
    if isinstance(v, dict):
        assert not set(s.get('required', [])) - v.keys()
        for k, x in v.items():
            if k in s.get('properties', {}):
                validate(x, s['properties'][k])
            else:
                assert s.get('additionalProperties') is not False, k

def main(path):
    result = json.loads(path.read_text())
    extension = result['state_contracts']
    validate(extension, json.loads((SCHEMAS / 'results-state.proposal.json').read_text()))
    # Negative control: validator must reject a hidden capacity-of-two error.
    bad = json.loads(json.dumps(extension))
    bad['spin']['channel_capacity_electrons'] = [2, 2]
    try:
        validate(bad, json.loads((SCHEMAS / 'results-state.proposal.json').read_text()))
    except AssertionError:
        pass
    else:
        raise AssertionError('schema negative control accepted')
    occ = extension['occupations']
    ns, nk, nb = occ['occupations_electrons']['shape']
    values = occ['occupations_electrons']['values']
    assert len(values) == ns * nk * nb
    counts = [sum(occ['k_weights'][k] * sum(values[(s*nk+k)*nb:(s*nk+k+1)*nb])
                  for k in range(nk)) for s in range(ns)]
    assert max(abs(a-b) for a, b in zip(counts, occ['channel_electrons'])) < 1e-10
    spec = importlib.util.spec_from_file_location('cube', ROOT / 'scripts/read_cube.py')
    reader = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reader)
    cubes = [reader.read_cube(path.parent / f'{name}.cube') for name in
             ('density-up', 'density-down', 'density-total', 'magnetization')]
    targets = counts + [sum(counts), counts[0]-counts[1]]
    for cube, target in zip(cubes, targets):
        assert abs(cube['integral']-target) < 1e-10
    for u, d, n, m in zip(*(c['values'] for c in cubes)):
        assert abs(u+d-n) < 1e-14 and abs(u-d-m) < 1e-14
    assert cubes[0]['values'] == result['density']['channel_values'][0]
    print('MF03 result proposal schemas, N/M sums and independent cube reader: PASS')

if __name__ == '__main__':
    main(Path(sys.argv[1]))
