"""Dependency-free checks of production force schemas and version compatibility."""
import json
from pathlib import Path
import sys
import importlib.util
spec = importlib.util.spec_from_file_location('spin_schema', Path(__file__).with_name('validate_spin_results.py'))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
validate, SCHEMAS = module.validate, module.SCHEMAS

def force_semantics(f):
    validate(f, json.loads((SCHEMAS / 'forces-available.proposal.json').read_text()))
    n = len(f['values'])
    assert all(len(v) == n for v in f['per_term'].values())
    for i in range(n):
        for a in range(3):
            assert abs(f['values'][i][a] - sum(v[i][a] for v in f['per_term'].values())) < 1e-12
    assert all(abs(f['net_force'][a] - sum(v[a] for v in f['values'])) < 1e-12 for a in range(3))
    assert abs(f['maximum_absolute_component'] - max(abs(x) for v in f['values'] for x in v)) < 1e-12

def main():
    schema = json.loads((SCHEMAS / 'results-state.proposal.json').read_text())
    directory = Path(__file__).parent / 'fixtures/force_schema'
    count = 0
    for path in sorted(directory.glob('*.json')):
        try:
            validate(json.loads(path.read_text()), schema)
            accepted = True
        except AssertionError:
            accepted = False
        assert accepted == path.name.startswith('accept-'), path
        count += 1
    if len(sys.argv) > 1:
        result = json.loads(Path(sys.argv[1]).read_text())
        if 'state_contracts' in result:
            validate(result['state_contracts'], schema)
        force_semantics(result['forces'])
    print(f'PASS force schema: {count} accepted/rejected fixtures and requested output')
if __name__ == '__main__':
    main()
