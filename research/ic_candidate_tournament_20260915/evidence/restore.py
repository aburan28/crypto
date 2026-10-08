#!/usr/bin/env python3
"""Restore hash-checked frozen evidence without changing existing files.

Requires Python 3.11+ and zstd. Gzip streams use zlib's level-9 raw DEFLATE;
their original headers, trailers and SHA-256 identities are preserved in PAX.
An incompatible zlib encoder fails closed rather than changing any receipt.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
import zlib

HERE = Path(__file__).resolve().parent


def restore(entry, destination):
    archive = HERE / entry['file']
    with archive.open('rb') as stream:
        actual = hashlib.file_digest(stream, 'sha256').hexdigest()
    if actual != entry['sha256']:
        raise ValueError(f'archive hash mismatch: {archive.name}')
    process = subprocess.Popen(['zstd', '-d', '-q', '-c', str(archive)], stdout=subprocess.PIPE)
    count = total = 0
    try:
        with tarfile.open(fileobj=process.stdout, mode='r|') as tar:
            for member in tar:
                relative = PurePosixPath(member.name)
                if not member.isfile() or relative.is_absolute() or '..' in relative.parts:
                    raise ValueError(f'unsafe archive member: {member.name}')
                target = destination / member.name
                if not target.resolve().is_relative_to(destination) or target.is_symlink():
                    raise ValueError(f'unsafe destination: {target}')
                content = tar.extractfile(member).read()
                metadata = member.pax_headers
                if 'ic.gzip.sha256' in metadata:
                    encoder = zlib.compressobj(9, zlib.DEFLATED, -15)
                    content = (bytes.fromhex(metadata['ic.gzip.header']) + encoder.compress(content)
                               + encoder.flush() + bytes.fromhex(metadata['ic.gzip.trailer']))
                    if hashlib.sha256(content).hexdigest() != metadata['ic.gzip.sha256']:
                        raise ValueError(f'gzip reconstruction differs with zlib {zlib.ZLIB_RUNTIME_VERSION}: {member.name}')
                if target.exists():
                    if target.read_bytes() != content:
                        raise ValueError(f'refusing to overwrite changed file: {target}')
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with target.open('xb') as stream:
                        stream.write(content)
                    target.chmod(0o755 if member.mode & 0o111 else 0o644)
                count += 1
                total += len(content)
        if process.wait() != 0:
            raise ValueError(f'zstd failed: {archive.name}')
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
        process.wait()
    if count != entry['files'] or total != entry['uncompressed_file_bytes']:
        raise ValueError(f'archive member count or length mismatch: {archive.name}')
    return {'archive': archive.name, 'status': 'RESTORED', 'files': count, 'bytes': total}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', help='One archive name, without .tar.zst; default: all')
    parser.add_argument('--out', type=Path, default=HERE.parent)
    args = parser.parse_args()
    manifest = json.loads((HERE / 'manifest.json').read_text())
    entries = [e for e in manifest['archives'] if args.archive is None or e['file'] == args.archive + '.tar.zst']
    if not entries:
        parser.error('unknown archive')
    destination = args.out.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    for entry in entries:
        print(json.dumps(restore(entry, destination)), flush=True)


if __name__ == '__main__':
    main()
