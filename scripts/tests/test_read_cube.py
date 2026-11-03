"""Reader checks independent of the Rust writer."""
import importlib.util
import pathlib
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('read_cube', pathlib.Path(__file__).resolve().parents[1] / 'read_cube.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CubeTests(unittest.TestCase):
    fixture = 'comment\nunits\n1 0 0 0\n2 1 0 0\n1 0.5 2 0\n2 0 0 -3\n14 4 1 2 3\n1 2 3 4\n'

    def read(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / 'test.cube'
            path.write_text(text)
            return MODULE.read_cube(path)

    def test_skew_left_handed(self):
        cube = self.read(self.fixture)
        self.assertEqual(cube['integral'], 60)
        self.assertEqual(cube['lattice_vectors_bohr'], [[2, 0, 0], [0.5, 2, 0], [0, 0, -6]])
        self.assertEqual(cube['atoms'][0]['atomic_number'], 14)
        self.assertEqual(cube['atoms'][0]['charge'], 4)
        self.assertEqual(cube['values'], [1, 2, 3, 4])

    def test_reject_malformed(self):
        for text in [self.fixture + '5', self.fixture.replace('1 2 3 4', '1 2 3'),
                     self.fixture.replace('1 2 3 4', '1 2 3 nan'),
                     self.fixture.replace('2 1 0 0', '-2 1 0 0'),
                     self.fixture.replace('1 0 0 0', '-1 0 0 0', 1),
                     self.fixture.replace('14 4', '119 4'),
                     self.fixture.replace('2 0 0 -3', '2 0 0 0')]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.read(text)
