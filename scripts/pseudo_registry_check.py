#!/usr/bin/env python3
"""Validate the repository schema's explicit subset with only the standard library.

Unsupported schema keywords fail closed; this is not a general JSON Schema engine.
"""
import argparse
import datetime
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / 'docs/schemas/pseudo-registry.schema.json'


def validate_schema(value, schema, path='$'):
    allowed = {'$schema', '$id', 'title', 'type', 'enum', 'minLength', 'pattern',
               'minimum', 'exclusiveMinimum', 'minItems', 'items', 'properties',
               'required', 'additionalProperties', 'patternProperties'}
    if set(schema) - allowed:
        raise ValueError(f'{path}: unsupported schema keywords {set(schema) - allowed}')
    types = schema.get('type', [])
    types = [types] if isinstance(types, str) else types
    checks = {'null': value is None, 'boolean': type(value) is bool,
              'integer': type(value) is int, 'number': type(value) in (int, float),
              'string': isinstance(value, str), 'array': isinstance(value, list),
              'object': isinstance(value, dict)}
    if types and not any(checks[t] for t in types):
        raise ValueError(f'{path}: expected {types}')
    if 'enum' in schema and value not in schema['enum']:
        raise ValueError(f'{path}: invalid enum')
    if isinstance(value, str):
        if len(value) < schema.get('minLength', 0) or ('pattern' in schema and not re.search(schema['pattern'], value)):
            raise ValueError(f'{path}: invalid string')
    if type(value) in (int, float):
        import math
        if not math.isfinite(value) or value < schema.get('minimum', -float('inf')) or value <= schema.get('exclusiveMinimum', -float('inf')):
            raise ValueError(f'{path}: invalid number')
    if isinstance(value, list):
        if len(value) < schema.get('minItems', 0):
            raise ValueError(f'{path}: too few items')
        for i, item in enumerate(value):
            validate_schema(item, schema['items'], f'{path}[{i}]')
    if isinstance(value, dict):
        if set(schema.get('required', [])) - value.keys():
            raise ValueError(f'{path}: missing required fields')
        for key, item in value.items():
            subs = ([schema['properties'][key]] if key in schema.get('properties', {}) else [])
            subs += [s for p, s in schema.get('patternProperties', {}).items() if re.search(p, key)]
            if not subs and schema.get('additionalProperties') is False:
                raise ValueError(f'{path}: unexpected field {key}')
            for sub in subs:
                validate_schema(item, sub, f'{path}.{key}')


def compatibility_reasons(e):
    reasons = []
    if e['status'] != 'verified' or e['sha256'] is None:
        reasons.append('File identity is unverified or deprecated')
    if 'UPF2' not in e['formats']:
        reasons.append('Reader requires UPF v2')
    if e['pseudo_type'] != 'NC':
        reasons.append('Only norm-conserving NC is supported; PAW/ultrasoft need augmentation')
    for key, label in [('nlcc', 'nonlinear core correction'), ('spin_orbit', 'spin-orbit'), ('paw', 'PAW'), ('ultrasoft', 'ultrasoft')]:
        if e[key] is not False:
            reasons.append(f'{label} unsupported or unverified')
    if e['relativistic'] not in ('scalar', 'nonrelativistic', 'no'):
        reasons.append('Relativistic treatment unsupported or unverified')
    counts = e['projector_counts']
    if counts is None or sum(counts.values()) > 3 or any(int(l) > 3 for l in counts):
        reasons.append('Reader supports at most three projectors TOTAL and l <= 3')
    xc = e['xc'] or ''
    # Match the current engine classifiers, including ASCII-only LDA whitespace.
    lda = re.split(r'[ \t\n\r\v\f]+', xc.strip(' \t\n\r\v\f'))
    if lda not in (['SLA', 'PZ'], ['SLA', 'PZ', 'NOGX', 'NOGC']) and ' '.join(xc.upper().split()) not in ('PBE', 'SLA PW PBX PBC'):
        reasons.append('XC is not an engine-recognized LDA-PZ or PBE label')
    if e['valence_electrons'] is None:
        reasons.append('Valence is unverified')
    return reasons


def validate(registry, schema=None):
    validate_schema(registry, schema or json.loads(SCHEMA.read_text()))
    ids = set()
    for e in registry['entries']:
        if e['id'] in ids:
            raise ValueError('duplicate id: ' + e['id'])
        ids.add(e['id'])
        datetime.date.fromisoformat(e['date_checked'])
        if e['sha256'] is None and e['status'] != 'unverified':
            raise ValueError('null hash requires unverified status')
        if e['status'] == 'verified' and (e['size_bytes'] is None or e['source']['file_url'] is None):
            raise ValueError('verified identity requires size and exact file URL')
        if e['license']['verified'] and e['license']['identifier'] is None:
            raise ValueError('verified license requires identifier')
        if e['supported_by_engine']:
            if compatibility_reasons(e) or e['support_reasons']:
                raise ValueError(e['id'] + ': contradictory engine support')
        elif not e['support_reasons']:
            raise ValueError(e['id'] + ': unsupported entries require reasons')
    return registry


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('registry', nargs='?', type=Path, default=ROOT / 'registry/pseudopotentials.json')
    args = p.parse_args()
    try:
        registry = validate(json.loads(args.registry.read_text()))
    except (ValueError, OSError) as exc:
        p.exit(1, f'Invalid registry: {exc}\n')
    print(f'PASS: {len(registry["entries"])} entries')


if __name__ == '__main__':
    main()
