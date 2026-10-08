#!/usr/bin/env python3
"""Rebuild the measured source snapshot into a new output directory."""
import argparse
import hashlib
import json
import pathlib
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--output', required=True, type=pathlib.Path)
args = parser.parse_args()
evidence = pathlib.Path(__file__).resolve().parent
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=False)
archive = evidence / 'measured-source.tar.gz'
manifest = json.loads((evidence / 'measured-source.json').read_text())
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archiveSha256']
source = out / 'source'
source.mkdir()
with tarfile.open(archive) as tf:
    tf.extractall(source, filter='data')
build = json.loads((evidence / 'build-config.json').read_text())
(out / 'build-config.json').write_text(json.dumps(build, indent=2))
for mode in range(4):
    with (out / f'build-{mode}.log').open('w') as log:
        subprocess.run(build + [f'PACKED_INLINE_POLY={mode}'], cwd=source,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    (source / 'ecc2k130').rename(out / f'client-{mode}')
