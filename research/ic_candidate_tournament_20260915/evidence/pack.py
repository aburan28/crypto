#!/usr/bin/env python3
"""Pack a finished round into the hash-checked archive format `restore.py`
reads, and record it in `manifest.json`.

Every gzip profile is stored expanded, with its original 10-byte header,
8-byte trailer and SHA-256 in PAX headers, after checking that zlib's raw
level-9 DEFLATE reproduces the original bytes exactly. Build caches, Python
bytecode and the operation lock are excluded. Member order, modes and
timestamps are normalised so the archive depends only on the files.

    python3 evidence/pack.py --round round-0006
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
EXCLUDED_DIRS = {'build', '__pycache__', 'preflight-next-build'}
EXCLUDED_FILES = {'operation.lock'}


def members(round_dir):
    for path in sorted(round_dir.rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if any(part in EXCLUDED_DIRS for part in relative.parts) or path.name in EXCLUDED_FILES:
            continue
        yield relative, path


def expand_gzip(content):
    """Expanded bytes plus PAX fields, or None when not losslessly repackable."""
    if len(content) < 18 or content[:2] != b'\x1f\x8b' or content[2] != 8 or content[3] != 0:
        return None
    header, body, trailer = content[:10], content[10:-8], content[-8:]
    expanded = gzip.decompress(content)
    encoder = zlib.compressobj(9, zlib.DEFLATED, -15)
    if encoder.compress(expanded) + encoder.flush() != body:
        return None
    return expanded, {'ic.gzip.header': header.hex(), 'ic.gzip.trailer': trailer.hex(),
                      'ic.gzip.sha256': hashlib.sha256(content).hexdigest()}


def pack(round_name, out):
    round_dir = ROOT / 'runs' / round_name
    if not round_dir.is_dir():
        raise SystemExit(f'no such round: {round_dir}')
    count = total = 0
    process = subprocess.Popen(['zstd', '-19', '-q', '-T2', '-o', str(out), '-f'], stdin=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=process.stdin, mode='w|', format=tarfile.PAX_FORMAT) as tar:
            for relative, path in members(round_dir):
                content = path.read_bytes()
                info = tarfile.TarInfo(str(PurePosixPath(relative)))
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                info.mode = 0o755 if path.stat().st_mode & 0o111 else 0o600
                stored = content
                if path.suffix == '.gz' and 'callgrind' in path.name:
                    expanded = expand_gzip(content)
                    if expanded is None:
                        raise SystemExit(f'profile is not losslessly repackable: {path}')
                    stored, info.pax_headers = expanded
                info.size = len(stored)
                tar.addfile(info, io.BytesIO(stored))
                count += 1
                total += len(content)
    finally:
        process.stdin.close()
        if process.wait() != 0:
            raise SystemExit('zstd failed')
    with out.open('rb') as stream:
        sha = hashlib.file_digest(stream, 'sha256').hexdigest()
    entry = {'file': out.name, 'sha256': sha, 'bytes': out.stat().st_size,
             'uncompressed_file_bytes': total, 'files': count}
    manifest_path = HERE / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest['archives'] = [e for e in manifest['archives'] if e['file'] != out.name] + [entry]
    manifest_path.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(entry))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--round', required=True, help='round directory name under runs/')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    out = args.out or HERE / f'{args.round}.tar.zst'
    pack(args.round, out.resolve())


if __name__ == '__main__':
    main()
