#!/usr/bin/env python3
"""Retain the finished native interface-control build and outputs; no execution."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import stat
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from oracle import require  # noqa: E402


def publish(build, controls, out):
    build, controls, out = map(Path, (build, controls, out))
    require(not out.exists(), 'publication output already exists')
    result = json.loads((controls/'result.json').read_text())
    require(result['status'] == 'PASS_ACTUAL_PREPARED_REPORT_INTERFACE_CONTROLS',
            'native interface controls did not pass')
    files = {}
    for prefix, root in (('build', build), ('controls', controls)):
        for path in sorted(root.rglob('*')):
            require(not path.is_symlink(), 'symlink in original evidence')
            if path.is_file():
                files[prefix+'/'+path.relative_to(root).as_posix()] = path
    out.mkdir(parents=True)
    inventory = []
    archive = out/'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for name, path in sorted(files.items()):
                value = path.read_bytes()
                mode = stat.S_IMODE(path.stat().st_mode)
                item = tarfile.TarInfo(name)
                item.size, item.mode, item.mtime = len(value), mode, 0
                tar.addfile(item, io.BytesIO(value))
                inventory.append(dict(role=name, bytes=len(value), mode=mode,
                                      sha256=hashlib.sha256(value).hexdigest()))
    require(all(hashlib.sha256(files[row['role']].read_bytes()).hexdigest() == row['sha256']
                for row in inventory), 'original evidence changed during capture')
    receipt = dict(schema_version=1, archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),
                   archive_bytes=archive.stat().st_size, inventory=inventory,
                   publication_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   native_solvers_executed_by_publication=0,
                   full_registry_dependency_sources_retained=False,
                   source_bound_scientific_runtime_admitted=False, online_speedup=None,
                   scope='repeatable native-build/interface controls; versioned full retained-input adapter remains pending')
    with (out/'receipt.json').open('x') as stream:
        stream.write(json.dumps(receipt, sort_keys=True, separators=(',', ':'))+'\n')
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build', 'controls', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    receipt = publish(args.build, args.controls, args.out)
    print(json.dumps({key:value for key,value in receipt.items() if key != 'inventory'}, sort_keys=True))
