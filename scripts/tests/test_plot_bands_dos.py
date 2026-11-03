import importlib.util
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

SPEC = importlib.util.spec_from_file_location("plot", Path(__file__).parents[1] / "plot_bands_dos.py")
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


class PlotTests(unittest.TestCase):
    def test_xml_bands_labels_and_dos(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bands = root / "bands.csv"
            bands.write_text("# energies eV, distance 1/Bohr\nk_index,distance_inv_bohr,kx,ky,kz,band_1,band_2\n0,0,0,0,0,-1,1\n1,1,0.5,0,0,0,2\n")
            (root / "bands_labels.csv").write_text("label,k_index,distance\nGamma,0,0\nX,1,1\n")
            dos = root / "dos.csv"
            dos.write_text("energy_eV,dos_states_per_eV_per_cell,integrated_dos\n-1,0,0\n0,2,1\n2,0,2\n")
            tree = ET.fromstring(PLOT.render(bands, dos))
            ns = {"s": "http://www.w3.org/2000/svg"}
            self.assertEqual(len(tree.findall("s:polyline", ns)), 3)
            labels = [t.text for t in tree.findall("s:text", ns)]
            self.assertIn("Gamma", labels)
            self.assertIn("X", labels)

    def test_breaks_escaping_and_reference_alignment(self):
        import json
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bands = root / 'bands.csv'
            bands.write_text('k_index,distance_inv_bohr,kx,ky,kz,band_1\n0,0,0,0,0,-1\n1,1,0,0,0,1\n2,1,0,0,0,2\n3,2,0,0,0,3\n')
            (root/'bands.json').write_text(json.dumps({'title':'A & <B>', 'path_break_indices':[2], 'input_sha256':'same', 'energy_reference_hartree':1, 'hartree_to_ev':27.211386245988034}))
            (root/'bands_labels.csv').write_text('label,k_index,distance\nA,0,0\nB,1,1\nC,2,1\nD,3,2\n')
            dos = root/'dos.csv'
            dos.write_text('energy_eV,dos_states_per_eV_per_cell,integrated_dos\n0,1,1\n1,2,2\n')
            (root/'dos.json').write_text(json.dumps({'input_sha256':'same','energy_reference_hartree':2}))
            tree = ET.fromstring(PLOT.render(bands,dos))
            ns = {'s':'http://www.w3.org/2000/svg'}
            self.assertEqual(len(tree.findall('s:polyline',ns)),3)
            labels=[t.text for t in tree.findall('s:text',ns)]
            self.assertIn('A & <B>',labels)
            self.assertIn('B|C',labels)
            # DOS coordinates move to the upper end of the shared energy axis.
            curve=tree.findall('s:polyline',ns)[-1].attrib['points']
            self.assertLess(float(curve.split()[0].split(',')[1]),100)
            (root/'dos.json').write_text(json.dumps({'input_sha256':'different','energy_reference_hartree':2}))
            with self.assertRaisesRegex(ValueError,'same input'):
                PLOT.render(bands,dos)

    def test_invalid_numeric_data(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bands.csv'
            for body in ['', 'distance_inv_bohr,x\n0,1\n', 'distance_inv_bohr,band_1\n0,nan\n1,0\n', 'distance_inv_bohr,band_1\n1,0\n0,1\n']:
                path.write_text(body)
                with self.assertRaises(ValueError):
                    PLOT.render(path)
