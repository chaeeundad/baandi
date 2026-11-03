"""Synthetic schema, consistency and builder tests; no external files or network."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]

def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

check = load('pseudo_registry_check')
build = load('pseudo_registry_build')


def supported():
    e = json.loads((ROOT / 'registry/pseudopotentials.json').read_text())['entries'][0]
    e.update(id='synthetic', element='Xx', formats=['UPF2'], xc='SLA PZ',
             valence_electrons=4, relativistic='scalar', pseudo_type='NC',
             nlcc=False, spin_orbit=False, paw=False, ultrasoft=False,
             projector_counts={'0': 2, '1': 1}, sha256='a' * 64,
             size_bytes=12, status='verified', supported_by_engine=True, support_reasons=[])
    e['source']['file_url'] = 'https://example.invalid/synthetic.upf'
    return dict(schema_version=1, entries=[e])


class RegistryTests(unittest.TestCase):
    def test_repository(self):
        check.validate(json.loads((ROOT / 'registry/pseudopotentials.json').read_text()))
        check.validate(supported())

    def test_schema_rejects_broken_fields(self):
        for key, value in [('sha256', 'made-up'), ('size_bytes', True), ('size_bytes', 0),
                           ('status', 'trusted'), ('valence_electrons', float('nan')),
                           ('projector_counts', {'s': 1}), ('date_checked', '2026-02-30'),
                           ('id', '../evil'), ('extra', 1), ('nlcc', 'false'),
                           ('citations', []), ('element', 'silicon')]:
            with self.subTest(key=key):
                r = supported(); r['entries'][0][key] = value
                with self.assertRaises(ValueError): check.validate(r)
        r = supported(); del r['entries'][0]['license']
        with self.assertRaises(ValueError): check.validate(r)
        r = supported(); r['entries'] *= 2
        with self.assertRaises(ValueError): check.validate(r)

    def test_all_consistency_guards(self):
        for key, value in [('sha256', None), ('nlcc', True), ('nlcc', None),
                           ('spin_orbit', True), ('paw', True), ('ultrasoft', True),
                           ('pseudo_type', 'USPP'), ('relativistic', 'full'),
                           ('projector_counts', {'0': 2, '1': 2}),
                           ('projector_counts', {'4': 1}), ('projector_counts', None),
                           ('xc', 'PBEsol'), ('formats', ['UPF']),
                           ('status', 'deprecated'), ('support_reasons', ['contradiction'])]:
            with self.subTest(key=key):
                r = supported(); r['entries'][0][key] = value
                with self.assertRaises(ValueError): check.validate(r)
        for label in ['SLA PZ', ' SLA\tPZ NOGX NOGC ', 'PBE', 'SLA PW PBX PBC']:
            r = supported(); r['entries'][0]['xc'] = label; check.validate(r)
        r = supported(); r['entries'][0]['sha256'] = None; r['entries'][0]['status'] = 'deprecated'
        with self.assertRaises(ValueError): check.validate(r)

    def test_schema_subset_fails_closed(self):
        with self.assertRaises(ValueError): check.validate_schema({}, {'oneOf': []})

    def test_header_facts_not_support(self):
        data = b'<UPF version="2.0.1"><PP_HEADER element="Xx" z_valence="4" functional="SLA PZ" core_correction="T" has_so="F"/><PP_NONLOCAL><PP_BETA.1 angular_momentum="0"/><PP_BETA.2 angular_momentum="0"/></PP_NONLOCAL></UPF>'
        facts = build.header_facts(data)
        self.assertTrue(facts['nlcc']); self.assertFalse(facts['spin_orbit'])
        self.assertIsNone(facts['paw']); self.assertEqual(facts['projector_counts'], {'0': 2})
        with self.assertRaises(ValueError): build.header_facts(b'<!DOCTYPE UPF []>' + data)
        with self.assertRaises(ValueError): build.header_facts(data.replace(b'</UPF>', b'<PP_HEADER/></UPF>'))
        # Even a parseable header cannot bypass strict numerical rejection.
        template = supported()['entries'][0]
        result = mock.Mock(stdout=json.dumps({'supported_by_engine': False, 'support_reasons': ['missing radial arrays']}))
        with mock.patch.object(build.subprocess, 'run', return_value=result):
            entry = build.build_entry(template, data, Path('/synthetic/cli'))
        self.assertFalse(entry['supported_by_engine']); self.assertIn('missing radial arrays', entry['support_reasons'])
        self.assertEqual(len(entry['sha256']), 64)

    def test_network_requires_https_and_rejects_redirects(self):
        with self.assertRaises(ValueError): build.download('http://example.invalid/test')
        with self.assertRaises(ValueError): build.NoRedirect().redirect_request(None, None, 302, '', {}, 'https://example.invalid/redirect')


if __name__ == '__main__':
    unittest.main()
