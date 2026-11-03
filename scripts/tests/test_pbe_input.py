"""PBE interchange, explicit overrides, and byte-identical legacy LDA input."""
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from make_qe_input import functional_input, generate
from qe_common import Error, read_json, results, write_json

FIXTURES = Path(__file__).resolve().parents[2] / 'tests/fixtures_format'


class FunctionalInputTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.source = self.root / 'results.json'
        self.output = self.root / 'qe.in'
        self.data = read_json(FIXTURES / 'results.json')

    def generate(self):
        write_json(self.source, self.data)
        generate(self.source, '/synthetic/pseudo', 'synthetic.UPF', self.output)
        return self.output.read_bytes()

    def pbe(self, label='SLA PZ', allow=True, applied=True):
        self.data.update(xc='pbe', allow_functional_override=allow,
                         functional_override_applied=applied)
        self.data['species'][0]['functional'] = label

    def test_lda_byte_identical(self):
        # Stored pre-PBE golden output: compare bytes, including comments/newlines.
        golden = (FIXTURES / 'golden.in').read_bytes()
        self.assertEqual(self.generate(), golden)
        self.data.update(allow_functional_override=False, functional_override_applied=False)
        self.assertEqual(self.generate(), golden)
        # Preserve the already implemented explicit LDA override path too.
        self.data.update(allow_functional_override=True, functional_override_applied=True)
        self.assertEqual(self.generate(), golden.replace(
            b'  nosym = .true.,\n', b"  nosym = .true.,\n  input_dft = 'PZ',\n"))

    def test_pbe_same_pseudopotential_override(self):
        original_species = self.data['species'].copy()
        self.pbe()
        text = self.generate().decode()
        self.assertIn("  input_dft = 'pbe',", text)
        self.assertIn('! Requested XC: pbe; UPF functional: "SLA PZ"', text)
        self.assertIn('! Explicit functional override: allow_functional_override=true; '
                      'functional_override_applied=true', text)
        self.assertEqual(results(self.source)[1]['species'], original_species)
        meta = read_json(str(self.output) + '.meta.json')
        self.assertEqual(meta['upf_sha256'], original_species[0]['sha256'])

    def test_matching_pbe_labels(self):
        for label in ('PBE', ' sla pw  pbx pbc '):
            for allow in (False, True):
                self.pbe(label, allow, False)
                text = self.generate().decode()
                self.assertIn("input_dft = 'pbe'", text)
                self.assertIn('functional_override_applied=false', text)

    def test_mismatch_without_explicit_override_rejected(self):
        self.pbe(allow=False)
        with self.assertRaisesRegex(Error, 'requires explicit'):
            self.generate()
        self.assertFalse(self.output.exists())

    def test_missing_invalid_or_inconsistent_provenance_rejected(self):
        for field in ('allow_functional_override', 'functional_override_applied'):
            for value in (None, 'true', 1):
                self.pbe()
                self.data[field] = value
                with self.assertRaisesRegex(Error, 'must be boolean'):
                    self.generate()
            self.pbe()
            del self.data[field]
            with self.assertRaisesRegex(Error, 'malformed results'):
                self.generate()
        self.pbe(applied=False)
        with self.assertRaisesRegex(Error, 'disagrees'):
            self.generate()
        self.pbe(label='PBE', applied=True)
        with self.assertRaisesRegex(Error, 'disagrees'):
            self.generate()

    def test_comment_label_cannot_inject_input(self):
        self.pbe(label="SLA PZ\n/\n&CONTROL")
        text = self.generate().decode()
        self.assertIn('UPF functional: "SLA PZ\\n/\\n&CONTROL"', text)
        self.assertEqual(text.count('&CONTROL'), 2)  # One escaped comment occurrence.
        self.assertEqual(text.splitlines().count('&CONTROL'), 1)

    def test_unknown_rejected(self):
        self.data['xc'] = 'pbe0'
        with self.assertRaises(Error):
            self.generate()
        with self.assertRaises(Error):
            functional_input(self.data)
