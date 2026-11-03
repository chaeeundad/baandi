"""Independent public-format checks: compare VASP/VTI exports with a cube.

Standard library only. This deliberately reads the scalar visualization subset,
not PAW restart files. Layouts: vasp.at/wiki/{POSCAR,CHGCAR,LOCPOT,ELFCAR};
docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html and vtkImageData
DirectionMatrix/TransformIndexToPhysicalPoint documentation. No code from those
programs is used. Execute with python3 -I SCRIPT CUBE VASP VTI QUANTITY.
"""
import importlib.util
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

SPEC = importlib.util.spec_from_file_location('cube_reader', Path(__file__).resolve().parents[3] / 'scripts/read_cube.py')
READER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(READER)
BOHR = 0.529177210903
# CODATA 2018 conversion used by dft-units/uom, independently stated here.
HARTREE_EV = 27.211386245988
SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
           'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm '
           'Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No '
           'Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()


def close(a, b):
    assert math.isfinite(a) and math.isfinite(b)
    assert math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-13), (a, b)


def read_vasp(path):
    lines = Path(path).read_text().splitlines()
    assert float(lines[1]) == 1.0
    lattice = [list(map(float, line.split())) for line in lines[2:5]]
    assert all(len(row) == 3 for row in lattice)
    symbols = lines[5].split()
    counts = list(map(int, lines[6].split()))
    assert len(symbols) == len(counts) and all(n > 0 for n in counts)
    assert lines[7] == 'Direct'
    n = sum(counts)
    positions = [list(map(float, line.split())) for line in lines[8:8+n]]
    assert all(len(row) == 3 for row in positions)
    assert not lines[8+n].strip()
    dims = list(map(int, lines[9+n].split()))
    assert len(dims) == 3 and min(dims) > 0
    values = [float(x) for line in lines[10+n:] for x in line.split()]
    assert len(values) == math.prod(dims) and all(map(math.isfinite, values))
    species = [SYMBOLS.index(s)+1 for s, count in zip(symbols, counts) for _ in range(count)]
    atoms = [(z, [sum(frac[j]*lattice[j][i] for j in range(3)) for i in range(3)])
             for z, frac in zip(species, positions)]
    return lattice, atoms, dims, values


def validate(cube_path, vasp_path, vti_path, quantity):
    cube = READER.read_cube(cube_path)
    lattice, atoms, dims, values = read_vasp(vasp_path)
    assert dims == cube['dims']
    for row, expected in zip(lattice, cube['lattice_vectors_bohr']):
        for x, y in zip(row, expected):
            close(x, y*BOHR)
    expected_atoms = [(a['atomic_number'], [x*BOHR for x in a['position_bohr']]) for a in cube['atoms']]
    assert len(atoms) == len(expected_atoms)
    # Stable grouping by element must preserve every atom without duplicates/loss.
    for z, position in atoms:
        index = next(i for i, (ez, ep) in enumerate(expected_atoms)
                     if ez == z and all(math.isclose(x, y, abs_tol=1e-12) for x, y in zip(position, ep)))
        expected_atoms.pop(index)
    volume = abs(READER.determinant(cube['lattice_vectors_bohr']))
    scale = volume if quantity in ('density', 'partial_charge') else HARTREE_EV if quantity in ('hartree', 'local_potential') else 1.0
    nx, ny, nz = dims
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                close(values[(k*ny+j)*nx+i], cube['values'][(i*ny+j)*nz+k]*scale)
    if quantity in ('density', 'partial_charge'):
        close(math.fsum(values)/len(values), cube['integral'])
    root = ET.parse(vti_path).getroot()
    assert root.tag == 'VTKFile' and root.attrib['type'] == 'ImageData'
    image = root.find('ImageData')
    extent = [0, nx, 0, ny, 0, nz]
    assert list(map(int, image.attrib['WholeExtent'].split())) == extent
    assert list(map(float, image.attrib['Origin'].split())) == [0, 0, 0]
    spacing = list(map(float, image.attrib['Spacing'].split()))
    direction = list(map(float, image.attrib['Direction'].split()))
    assert len(spacing) == 3 and len(direction) == 9
    # Row-major transform columns times spacing reconstruct the actual voxel
    # vectors, including non-orthogonal and negative-determinant geometry.
    for cart in range(3):
        for axis in range(3):
            close(direction[3*cart+axis]*spacing[axis]*dims[axis], cube['lattice_vectors_bohr'][axis][cart]*BOHR)
    piece = image.find('Piece')
    assert list(map(int, piece.attrib['Extent'].split())) == extent
    point_data = piece.find('PointData')
    array = point_data.find('DataArray')
    assert array.attrib['type'] == 'Float64' and array.attrib['format'] == 'ascii'
    assert point_data.attrib['Scalars'] == array.attrib['Name']
    units = {'density': 'electrons_per_bohr3', 'partial_charge': 'per_bohr3', 'elf': 'dimensionless',
             'hartree': 'hartree', 'local_potential': 'hartree'}
    assert array.attrib['Name'] == quantity + '_' + units[quantity]
    samples = list(map(float, array.text.split()))
    assert len(samples) == (nx+1)*(ny+1)*(nz+1)
    for k in range(nz+1):
        for j in range(ny+1):
            for i in range(nx+1):
                close(samples[(k*(ny+1)+j)*(nx+1)+i], cube['values'][((i%nx)*ny+(j%ny))*nz+(k%nz)])
    return dict(quantity=quantity, dims=dims, cube_integral=cube['integral'],
                vasp_mean=math.fsum(values)/len(values), minimum=cube['minimum'], maximum=cube['maximum'])


if __name__ == '__main__':
    print(json.dumps(validate(*sys.argv[1:]), allow_nan=False))
