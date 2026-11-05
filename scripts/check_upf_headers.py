#!/usr/bin/env python3
"""Independent stdlib UPF metadata reader; never prints radial samples.

DFT_NC_DIR selects external files. --binary compares exactly with Rust inspection.
A successful header comparison is not numerical support or qualification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def inspect(path):
    data = path.read_bytes()
    root = ET.fromstring(data)
    h = root.find('PP_HEADER').attrib
    def number(key, cast):
        return cast(h[key]) if key in h else None
    def flag(value):
        if value.strip().lower() in ('t', '.t.', 'true', '.true.', '1'):
            return True
        if value.strip().lower() in ('f', '.f.', 'false', '.false.', '0'):
            return False
        raise ValueError('invalid core_correction flag')
    return hashlib.sha256(data).hexdigest(), {
        'element': h['element'].strip(), 'z_valence': float(h['z_valence']),
        'pseudo_type': h['pseudo_type'], 'relativistic': h['relativistic'],
        'functional': h['functional'], 'nlcc': flag(h['core_correction']),
        'number_of_wfc': int(h['number_of_wfc']), 'mesh_size': int(h['mesh_size']),
        'number_of_proj': int(h['number_of_proj']),
        'projector_channels': [{'index': int(b.attrib['index']),
                                'l': int(b.attrib['angular_momentum'])}
                               for b in root.find('PP_NONLOCAL')
                               if b.tag.startswith('PP_BETA.')],
        'cutoff_hints_ry': {'wavefunction': number('wfc_cutoff', float),
                            'density': number('rho_cutoff', float)},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--directory', default=os.environ.get('DFT_NC_DIR'))
    ap.add_argument('--manifest', type=Path, default=Path('docs/data/nc-datasets.json'))
    ap.add_argument('--binary', type=Path)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    if not args.directory:
        ap.error('set DFT_NC_DIR or --directory')
    manifest = json.loads(args.manifest.read_text())
    results = []
    for entry in manifest['datasets']:
        path = Path(args.directory) / entry['filename']
        sha, facts = inspect(path)
        if sha != entry['sha256'] or facts['element'] != entry['element']:
            raise ValueError(f'{path.name}: manifest identity mismatch')
        for key in ('z_valence', 'functional', 'relativistic', 'pseudo_type', 'nlcc'):
            if entry[key] != facts[key]:
                raise ValueError(f'{path.name}: manifest mismatch in {key}')
        counts = {}
        for channel in facts['projector_channels']:
            key = str(channel['l'])
            counts[key] = counts.get(key, 0) + 1
        if counts != entry['projector_counts']:
            raise ValueError(f'{path.name}: manifest projector counts mismatch')
        record = {'element': facts['element'], 'manifest_hash_match': True,
                  'header_facts': facts, 'qualification_status': 'unqualified'}
        if args.binary:
            run = subprocess.run([str(args.binary.resolve()), 'inspect-upf', str(path),
                                  '--json', '--manifest', str(args.manifest)],
                                 capture_output=True, text=True, check=True)
            rust = json.loads(run.stdout)
            for key, expected in [('sha256', sha), ('header_facts', facts),
                                  ('license_reference', entry['license']),
                                  ('qualification_status', entry['qualification_status'])]:
                if rust[key] != expected:
                    raise ValueError(f'{path.name}: Rust/Python mismatch in {key}')
            record.update(rust_python_equal=True, blockers=rust['blockers'],
                          supported_by_engine=rust['supported_by_engine'])
        results.append(record)
    report = {'schema_version': 1, 'scope': 'Header identity only; no radial data, SCF or qualification',
              'command': 'DFT_NC_DIR=<external directory> python3 -I scripts/check_upf_headers.py --binary target/debug/dft-cli --output tests/reference/datasets/pseudodojo-headers.json',
              'datasets': results}
    output = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
