#!/usr/bin/env python3
"""Restore byte-identical benchmark evidence, including original gzip envelopes."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import zlib

HERE = Path(__file__).resolve().parent
FORMAT = 'ic-audit-bundle-v1'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def reconstruct(data, entry):
    if entry['codec'] == 'gzip-deflate9':
        compressor = zlib.compressobj(9, zlib.DEFLATED, -15)
        data = (bytes.fromhex(entry['header']) + compressor.compress(data)
                + compressor.flush() + bytes.fromhex(entry['trailer']))
    elif entry['codec'] != 'identity':
        raise ValueError('Unknown evidence codec')
    if len(data) != entry['bytes'] or digest(data) != entry['sha256']:
        raise ValueError('Evidence reconstruction failed: ' + entry['path']
                         + '; recorded zlib=' + entry.get('zlib', 'none')
                         + ', current zlib=' + zlib.ZLIB_RUNTIME_VERSION)
    return data


def restore(archive, expected_sha256, destination):
    if digest(archive.read_bytes()) != expected_sha256:
        raise ValueError('Archive checksum mismatch: ' + str(archive))
    destination = destination.resolve()
    seen = set()
    with tarfile.open(archive, 'r|xz') as tar:
        first = next(iter(tar))
        if first.name != 'index.json' or not first.isfile():
            raise ValueError('Missing evidence index')
        index = json.load(tar.extractfile(first))
        if index['format'] != FORMAT:
            raise ValueError('Unsupported evidence format')
        entries = {entry['path']: entry for entry in index['files']}
        if len(entries) != len(index['files']):
            raise ValueError('Duplicate evidence paths')
        for member in tar:
            if member.name == 'index.json':
                continue
            if not member.isfile() or not member.name.startswith('payload/'):
                raise ValueError('Unexpected archive entry')
            name = member.name[len('payload/'):]
            if name not in entries or name in seen:
                raise ValueError('Unexpected or repeated evidence path')
            relative = Path(name)
            if relative.is_absolute() or '..' in relative.parts or relative.parts[0] != 'research':
                raise ValueError('Unsafe evidence path')
            target = destination / relative
            if not target.resolve().is_relative_to(destination):
                raise ValueError('Evidence path leaves destination')
            data = reconstruct(tar.extractfile(member).read(), entries[name])
            if target.exists():
                if not target.is_file() or target.read_bytes() != data:
                    raise ValueError('Refusing to overwrite different file: ' + str(target))
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open('xb') as output:
                    output.write(data)
                target.chmod(entries[name]['mode'] & 0o777)
            seen.add(name)
        if seen != set(entries):
            raise ValueError('Incomplete evidence archive')
    return {'archive': archive.name, 'files': len(seen), 'status': 'VERIFIED'}


def restore_bundle(row, source, destination):
    if 'parts' not in row:
        return restore(source / row['archive'], row['sha256'], destination)
    with tempfile.TemporaryDirectory(prefix='ic-evidence-') as temporary:
        archive = Path(temporary) / 'bundle.tar.xz'
        with archive.open('wb') as output:
            for part in row['parts']:
                data = (source / part['file']).read_bytes()
                if len(data) != part['bytes'] or digest(data) != part['sha256']:
                    raise ValueError('Archive part checksum mismatch: ' + part['file'])
                output.write(data)
        if archive.stat().st_size != row['bytes']:
            raise ValueError('Incomplete multipart evidence archive')
        result = restore(archive, row['sha256'], destination)
        return {**result, 'archive': row['archive'], 'parts': len(row['parts'])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, default=HERE.parents[1])
    parser.add_argument('--only', help='Restore one named bundle')
    args = parser.parse_args()
    manifest = json.loads((HERE / 'bundles.json').read_text())
    selected = [row for row in manifest['archives'] if args.only is None or row['name'] == args.only]
    if not selected:
        parser.error('No matching evidence bundle')
    for row in selected:
        print(json.dumps(restore_bundle(row, HERE, args.destination)), flush=True)


if __name__ == '__main__':
    main()
