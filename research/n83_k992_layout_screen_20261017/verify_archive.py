#!/usr/bin/env python3
"""Check sealed K992 layout archive membership and exact member bytes."""

import hashlib
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
archive = HERE / 'raw.tar.gz'
manifest_blob = (HERE / 'manifest.json').read_bytes()
manifest = json.loads(manifest_blob)
if manifest['schema'] != 'n83.local-k992-layout-archive/v1':
    raise ValueError('archive schema differs')
expected = manifest['files']
seen = set()
with tarfile.open(archive, 'r:gz') as saved:
    for item in saved:
        name = item.name
        path = Path(name)
        if (name in seen or path.is_absolute() or '..' in path.parts or
                not item.isfile() or name not in expected and name != 'manifest.json'):
            raise ValueError(f'unsafe or unlisted member: {name}')
        seen.add(name)
        blob = saved.extractfile(item).read()
        if name == 'manifest.json':
            if blob != manifest_blob:
                raise ValueError('manifest member differs')
        elif (hashlib.sha256(blob).hexdigest() != expected[name] or
              name in ('replay.json', 'screen.svg', 'screen.png', 'screen.pdf')
              and blob != (HERE / name).read_bytes()):
            raise ValueError(f'archive member differs: {name}')
if seen != set(expected) | {'manifest.json'}:
    raise ValueError('archive member inventory incomplete')
print(json.dumps({'status': 'PASS_EXACT_ARCHIVE_MEMBERS',
                  'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
                  'members': len(seen)}, sort_keys=True))
