#!/usr/bin/env python3
"""Pack a finished round into the hash-checked evidence format `restore.py` reads.

Regular files are stored unchanged. Each gzip Callgrind profile is stored
expanded so zstd can share content across profiles; its 10-byte header,
8-byte trailer and SHA-256 go into PAX fields, and the reconstruction is
checked here before the member is written. Build caches, Python bytecode and
operation locks are excluded, as in the committed archives.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXCLUDED_DIRS = {'build', '__pycache__', 'preflight-next-build', 'round6-preflight-build'}
EXCLUDED_FILES = {'operation.lock'}


def members(paths):
    for top in paths:
        for path in sorted((ROOT / top).rglob('*')):
            if not path.is_file():
                continue
            relative = path.relative_to(ROOT)
            if set(relative.parts[:-1]) & EXCLUDED_DIRS or relative.name in EXCLUDED_FILES:
                continue
            yield relative, path


def pack(name, paths):
    archive = HERE / (name + '.tar.zst')
    if archive.exists():
        raise SystemExit(f'refusing to overwrite {archive}')
    process = subprocess.Popen(['zstd', '-19', '--long=27', '-q', '-o', str(archive)], stdin=subprocess.PIPE)
    count = total = 0
    with tarfile.open(fileobj=process.stdin, mode='w|', format=tarfile.PAX_FORMAT) as tar:
        for relative, path in members(paths):
            data = path.read_bytes()
            info = tarfile.TarInfo(str(relative))
            info.mode = 0o755 if path.stat().st_mode & 0o111 else 0o644
            content = data
            if relative.name.startswith('callgrind.out') and relative.suffix == '.gz':
                expanded = gzip.decompress(data)
                encoder = zlib.compressobj(9, zlib.DEFLATED, -15)
                rebuilt = data[:10] + encoder.compress(expanded) + encoder.flush() + data[-8:]
                if rebuilt != data:
                    raise SystemExit(f'profile is not reproducible with zlib {zlib.ZLIB_RUNTIME_VERSION}: {relative}')
                info.pax_headers = {'ic.gzip.sha256': hashlib.sha256(data).hexdigest(),
                                    'ic.gzip.header': data[:10].hex(), 'ic.gzip.trailer': data[-8:].hex()}
                content = expanded
            info.size = len(content)
            tar.addfile(info, io.BytesIO(content))
            count += 1
            total += len(data)
    process.stdin.close()
    if process.wait() != 0:
        raise SystemExit('zstd failed')
    with archive.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    entry = {'file': archive.name, 'sha256': sha, 'bytes': archive.stat().st_size,
             'uncompressed_file_bytes': total, 'files': count}
    manifest_path = HERE / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    if any(e['file'] == entry['file'] for e in manifest['archives']):
        raise SystemExit('manifest already lists this archive')
    manifest['archives'].append(entry)
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(entry))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--name', required=True, help='archive name without .tar.zst')
    parser.add_argument('paths', nargs='+', help='directories relative to the research directory')
    args = parser.parse_args()
    pack(args.name, args.paths)


if __name__ == '__main__':
    main()
