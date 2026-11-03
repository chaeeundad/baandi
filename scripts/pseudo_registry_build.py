#!/usr/bin/env python3
"""Maintainer-only registry builder; metadata template supplies legal provenance.

Inputs are a directory of UPF files (case insensitive suffix) or explicit URLs.
Network use requires --allow-network. No archive extraction or redirects.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

_spec = importlib.util.spec_from_file_location('pseudo_registry_check', Path(__file__).with_name('pseudo_registry_check.py'))
checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(checker)
MAX_BYTES = 64 * 1024 * 1024


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('redirect refused: ' + newurl)


def download(url):
    if urllib.parse.urlsplit(url).scheme != 'https':
        raise ValueError('HTTPS required')
    with urllib.request.build_opener(NoRedirect).open(url, timeout=30) as response:
        if response.status != 200:
            raise ValueError('HTTP status is not 200')
        data = response.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError('file exceeds 64 MiB limit')
    return data


def header_facts(data):
    # Header-only facts do not establish support; the strict Rust reader does.
    if b'<!DOCTYPE' in data.upper() or b'<!ENTITY' in data.upper():
        raise ValueError('DTD/entities forbidden')
    root = ET.fromstring(data)
    headers = root.findall('PP_HEADER')
    if root.tag != 'UPF' or len(headers) != 1:
        raise ValueError('expected one UPF/PP_HEADER')
    h = headers[0].attrib
    def flag(key):
        value = h.get(key, '').strip().lower()
        if value in ('t', 'true', '.true.', '1'): return True
        if value in ('f', 'false', '.false.', '0'): return False
        return None
    counts = {}
    for beta in root.findall('PP_NONLOCAL/*'):
        if beta.tag.startswith('PP_BETA.'):
            l = str(int(beta.attrib['angular_momentum']))
            counts[l] = counts.get(l, 0) + 1
    version = root.attrib.get('version')
    return dict(element=h.get('element', '').strip(), valence_electrons=float(h['z_valence']),
                xc=h.get('functional'), relativistic=h.get('relativistic'),
                pseudo_type=h.get('pseudo_type'), nlcc=flag('core_correction'),
                spin_orbit=flag('has_so'), paw=flag('is_paw'), ultrasoft=flag('is_ultrasoft'),
                projector_counts=counts, formats=['UPF2'] if version in ('2.0', '2.0.0', '2.0.1') else ['UPF'])


def build_entry(template, data, cli):
    e = json.loads(json.dumps(template))
    e.update(header_facts(data))
    # Temporary bytes stay outside the repository and are removed even on error.
    with tempfile.TemporaryDirectory(prefix='baandi-inspect-') as tmp:
        path = Path(tmp) / 'input.upf'
        path.write_bytes(data)
        result = subprocess.run([str(cli), 'inspect-upf', str(path), '--json'], capture_output=True, text=True, check=True)
    inspection = json.loads(result.stdout)
    if type(inspection.get('supported_by_engine')) is not bool or not isinstance(inspection.get('support_reasons'), list):
        raise ValueError('invalid inspect-upf JSON')
    # Accepted files use the strict reader's facts, including original XC label.
    for key in ('element', 'valence_electrons', 'xc', 'relativistic', 'pseudo_type', 'nlcc', 'spin_orbit', 'paw', 'ultrasoft', 'projector_counts', 'formats'):
        if key in inspection:
            e[key] = inspection[key]
    e.update(sha256=hashlib.sha256(data).hexdigest(), size_bytes=len(data),
             date_checked=datetime.date.today().isoformat(), status='verified')
    reasons = list(inspection['support_reasons']) + checker.compatibility_reasons(e)
    if not inspection['supported_by_engine'] and not reasons:
        reasons.append('Strict reader did not establish support')
    e['support_reasons'] = list(dict.fromkeys(reasons))
    e['supported_by_engine'] = not e['support_reasons']
    checker.validate(dict(schema_version=1, entries=[e]))
    return e


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs', nargs='+', help='One directory, local files, or HTTPS URLs')
    p.add_argument('--metadata', type=Path, required=True, help='JSON map from filename or exact URL to complete entry template')
    p.add_argument('--cli', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--allow-network', action='store_true')
    args = p.parse_args()
    try:
        metadata = json.loads(args.metadata.read_text())
        inputs = []
        for value in args.inputs:
            path = Path(value)
            if path.is_dir():
                inputs.extend(str(f) for f in sorted(path.iterdir()) if f.is_file() and f.suffix.lower() == '.upf')
            else:
                inputs.append(value)
        if not inputs:
            raise ValueError('no UPF files found')
        entries = []
        for value in inputs:
            remote = urllib.parse.urlsplit(value).scheme in ('http', 'https')
            template = metadata[value if remote else Path(value).name]
            if not template['source']['file_url']:
                raise ValueError('exact original file URL required in metadata')
            if remote:
                if not args.allow_network:
                    raise ValueError('URL input requires --allow-network')
                if template['license']['verified'] is not True or not template['license']['identifier']:
                    raise ValueError('license unverified: read source page before downloading')
                if template['source']['file_url'] != value:
                    raise ValueError('metadata URL does not match requested URL')
                data = download(value)
            else:
                with Path(value).open('rb') as stream:
                    data = stream.read(MAX_BYTES + 1)
                if len(data) > MAX_BYTES:
                    raise ValueError('file exceeds 64 MiB limit')
            entries.append(build_entry(template, data, args.cli.resolve()))
        registry = checker.validate(dict(schema_version=1, entries=entries))
        # Refuse replacement: maintainer reviews a separate candidate registry.
        with args.output.open('x') as stream:
            json.dump(registry, stream, indent=2)
            stream.write('\n')
        print(f'Wrote {len(entries)} inspected entries to {args.output}')
    except (ValueError, KeyError, OSError, ET.ParseError, subprocess.CalledProcessError) as exc:
        p.exit(1, f'Build failed: {exc}\n')


if __name__ == '__main__':
    main()
