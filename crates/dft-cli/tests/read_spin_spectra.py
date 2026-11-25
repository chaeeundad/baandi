"""R4b independent stdlib reader: explicit spin axes, portable tables and VASP layouts.

DOSCAR's five-column definition: https://vasp.at/wiki/DOSCAR .
EIGENVAL ISPIN=2 rows: band index, epsilon_up/down (eV), occupation_up/down.
CSV/dat/Grace/gnuplot/SVG bands use channel-major columns; total DOS is spin-summed.
"""
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

root, command = Path(sys.argv[1]), sys.argv[2]
meta = json.loads((root / (command + '.json')).read_text())
nb = meta['nbands']
nk = len(meta['kpoints'])

def close(a, b):
    assert abs(a-b) < 1e-10, (a,b)

def rows(name):
    return [[float(v) for v in line.split()] for line in (root / name).read_text().splitlines()
            if line.strip() and not line.startswith(('#','@','&'))]

with (root / (command + '.csv')).open() as stream:
    table = list(csv.DictReader(line for line in stream if not line.startswith('#')))
assert meta['band_data']['shape'] == [2,nk,nb]
assert meta['spin']['channel_order'] == ['up','down']
ET.parse(root / (command + '.svg'))
if command == 'bands':
    assert len(table) == nk
    assert list(table[0]) == ['k_index','distance_inv_bohr','kx','ky','kz'] + [f'band_{b+1}' for b in range(2*nb)]
    lines = (root/'EIGENVAL').read_text().splitlines()
    assert lines[0].split()[2:] == ['1','2']
    assert list(map(int,lines[5].split())) == [int(meta['electrons']),nk,nb]
    cursor = 6
    dat = rows('bands.dat')
    expected_blocks = []
    for k, row in enumerate(table):
        assert lines[cursor] == ''
        point = list(map(float,lines[cursor+1].split()))
        assert point[3] == 0
        for a,b in zip(point[:3],meta['kpoints'][k]['fractional']): close(a,b)
        cursor += 2
        for s in range(2):
            for b in range(nb):
                energy = (meta['band_data']['eigenvalues_ha'][s][k][b]-meta['energy_reference_hartree'])*meta['hartree_to_ev']
                close(float(row[f'band_{s*nb+b+1}']),energy)
                close(dat[k][1+s*nb+b],energy)
        for b in range(nb):
            values = list(map(float,lines[cursor].split())); cursor += 1
            assert len(values) == 5 and values[0] == b+1
            for s in range(2):
                close(values[1+s],float(row[f'band_{s*nb+b+1}']))
                assert values[3+s] == (1 if b < meta['spin']['occupied_bands'][s] else 0)
    assert cursor == len(lines)
    boundaries = [0]+meta['path_break_indices']+[nk]
    for b in range(2*nb):
        for start,end in zip(boundaries,boundaries[1:]):
            expected_blocks.extend([[float(r['distance_inv_bohr']),float(r[f'band_{b+1}'])] for r in table[start:end]])
    for name in ['bands.dat.gnu','bands.agr']:
        data = rows(name)
        assert len(data)==len(expected_blocks)
        for a,b in zip(data,expected_blocks):
            assert len(a)==2
            for x,y in zip(a,b): close(x,y)
    ns = {'s':'http://www.w3.org/2000/svg'}
    assert len(ET.parse(root/'bands.svg').findall('s:polyline',ns)) == 2*nb*(len(boundaries)-1)
else:
    lines = (root/'DOSCAR').read_text().splitlines()
    assert int(lines[5].split()[2]) == len(table)
    assert float(lines[5].split()[3]) == 0
    dat, agr = rows('dos.dat'), rows('dos.agr')
    for i,(line,row) in enumerate(zip(lines[6:],table)):
        values = list(map(float,line.split()))
        assert len(values)==5
        close(values[0],float(row['energy_eV']))
        close(values[1]+values[2],float(row['dos_states_per_eV_per_cell']))
        close(values[3]+values[4],float(row['integrated_dos']))
        for a,b in zip(dat[i],row.values()): close(a,float(b))
        for a,b in zip(agr[i],list(row.values())[:2]): close(a,float(b))
    close(float(table[-1]['integrated_dos']),2*nb)
