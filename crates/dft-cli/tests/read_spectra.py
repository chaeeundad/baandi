"""Independent standard-library schema reader; no production writer imports.

Checks spin-unpolarized EIGENVAL (six-line header, k blocks), total-only
DOSCAR (six-line header and NEDOS triples), Grace XY datasets and gnuplot
blocks against CSV. These are format checks, not third-party application tests.
"""
import builtins
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

def zip(*items, strict=False):
    if strict:
        items = [list(item) for item in items]
        assert len({len(item) for item in items}) <= 1
    return builtins.zip(*items)

root = Path(sys.argv[1])
def table(name):
    with (root / name).open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))
def numbers(line):
    return [float(x) for x in line.split()]
def data(name):
    return [numbers(line) for line in (root / name).read_text().splitlines() if line.strip() and not line.startswith(('#', '@', '&'))]
def close(a, b):
    assert abs(a-b) < 1e-11, (a, b)
bands = table('bands.csv')
dos = table('dos.csv')
meta = json.loads((root / 'bands.json').read_text())
nb = meta['nbands']
nk = len(bands)
assert list(bands[0]) == ['k_index','distance_inv_bohr','kx','ky','kz'] + [f'band_{b+1}' for b in range(nb)]
assert list(dos[0]) == ['energy_eV','dos_states_per_eV_per_cell','integrated_dos']
assert '1/Bohr' in (root / 'bands.csv').read_text().splitlines()[0]
assert 'states/eV/cell' in (root / 'dos.csv').read_text().splitlines()[0]
assert len(table('bands_labels.csv')) == len(meta['path_labels'])
for row, expected in zip(data('bands.dat'), bands, strict=True):
    assert len(row) == nb+1
    close(row[0], float(expected['distance_inv_bohr']))
    for b in range(nb): close(row[b+1], float(expected[f'band_{b+1}']))
for name in ['dos.dat','DOSCAR']:
    rows = data(name) if name == 'dos.dat' else [numbers(s) for s in (root/name).read_text().splitlines()[6:]]
    for row, expected in zip(rows, dos, strict=True):
        assert len(row)==3
        for a,b in zip(row, expected.values(), strict=True): close(a,float(b))
header = (root/'DOSCAR').read_text().splitlines()
assert numbers(header[0])[2] == 0
assert int(numbers(header[5])[2]) == len(dos)
assert numbers(header[5])[3] == 0 # VBM placeholder, not Fermi energy
lines = (root/'EIGENVAL').read_text().splitlines()
assert list(map(int, lines[0].split()[2:])) == [1, 1]
ne, n, count = map(int, lines[5].split())
assert ne == int(meta['electrons']) and n == nk and count == nb
cursor = 6
for expected in bands:
    assert lines[cursor] == ''; cursor += 1
    point = numbers(lines[cursor]); cursor += 1
    assert len(point)==4 and point[3]==0
    for a,key in zip(point[:3], ['kx','ky','kz']): close(a,float(expected[key]))
    for b in range(nb):
        row=numbers(lines[cursor]);cursor+=1
        assert len(row)==3 and row[0]==b+1
        close(row[1],float(expected[f'band_{b+1}']))
        assert row[2] == (2 if b < ne//2 else 0)
assert cursor==len(lines)
blocks=[]; current=[]
for line in (root/'bands.dat.gnu').read_text().splitlines():
    if line.startswith('#'): continue
    if not line.strip():
        if current: blocks.append(current);current=[]
    else: current.append(numbers(line))
if current: blocks.append(current)
boundaries=[0]+meta['path_break_indices']+[nk]
assert len(blocks)==nb*(len(boundaries)-1)
index=0
for b in range(nb):
    for start,end in zip(boundaries,boundaries[1:]):
        block=blocks[index];index+=1
        assert len(block)==end-start
        for pair,row in zip(block,bands[start:end]):
            assert len(pair)==2
            close(pair[0],float(row['distance_inv_bohr']));close(pair[1],float(row[f'band_{b+1}']))
for name,expected in [('bands.agr',[r for block in blocks for r in block]),('dos.agr',[[float(r['energy_eV']),float(r['dos_states_per_eV_per_cell'])] for r in dos])]:
    text=(root/name).read_text()
    assert '@version' in text and '@type xy' in text and '&' in text
    for row,pair in zip(data(name),expected,strict=True):
        assert len(row)==2
        close(row[0],pair[0]);close(row[1],pair[1])
ns={'s':'http://www.w3.org/2000/svg'}
assert len(ET.parse(root/'bands.svg').findall('s:polyline',ns)) == len(blocks)
ET.parse(root/'dos.svg')
dmeta=json.loads((root/'dos.json').read_text())
close(sum(k['weight'] for k in dmeta['kpoints']),1)
close(float(dos[-1]['integrated_dos']),2*dmeta['nbands'])
