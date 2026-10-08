#!/usr/bin/env python3
"""Retain and restore these original incomplete controls without rerunning them."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import stat
import sys
import tarfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from identity import write_immutable  # noqa: E402
from oracle import require  # noqa: E402
from sat_runtime_execution_v3 import make_read_only  # noqa: E402


def digest(data):
    return hashlib.sha256(data).hexdigest()


def pack(root, diagnosis, out):
    root, diagnosis, out = map(Path, (root, diagnosis, out))
    require(not out.exists(), 'publication output already exists')
    files = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'publication cannot contain symlinks')
        if path.is_file():
            files[path.relative_to(root).as_posix()] = path
    files['diagnosis/diagnosis.json'] = diagnosis
    require(all(name in files for name in (
        'f5-registration/execution-claim.json', 'sat-registration/execution-claim.json',
        'f5-execution/process.json', 'sat-execution/process.json',
        'f5-audit/transport.json', 'sat-audit/transport.json')), 'publication lost a terminal control or claim')
    out.mkdir(parents=True)
    inventory = []
    archive = out/'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for name, path in sorted(files.items()):
                data = path.read_bytes()
                mode = stat.S_IMODE(path.stat().st_mode)
                item = tarfile.TarInfo(name)
                item.size, item.mode, item.mtime = len(data), mode, 0
                tar.addfile(item, io.BytesIO(data))
                inventory.append(dict(role=name, bytes=len(data), sha256=digest(data), mode=mode))
    # Native jobs are terminal. Reject a publication if any source/output byte
    # changed while copying it; preserve the rejected archive rather than rerun.
    require(all(path.read_bytes() == (root/name).read_bytes()
                for name, path in files.items() if name != 'diagnosis/diagnosis.json'),
            'publication input moved during capture')
    require(all(digest(files[row['role']].read_bytes()) == row['sha256'] for row in inventory),
            'publication bytes changed during capture')
    data = archive.read_bytes()
    receipt = dict(schema_version=1, status='RETAINED_ORIGINAL_INCOMPLETE_CONTROLS',
                   archive_sha256=digest(data), archive_bytes=len(data), inventory=inventory,
                   native_solvers_executed_by_publication=0, complete_ic_admitted=False,
                   promotion_eligible=False, online_speedup=None)
    write_immutable(out/'receipt.json', receipt)
    return receipt


def restore(bundle, out):
    bundle, out = map(Path, (bundle, out))
    require(not out.exists(), 'restoration output already exists')
    receipt = json.loads((bundle/'receipt.json').read_text())
    data = (bundle/'evidence.tar.gz').read_bytes()
    require(digest(data) == receipt['archive_sha256'] and len(data) == receipt['archive_bytes'],
            'published archive digest or byte count changed')
    expected = {row['role']: row for row in receipt['inventory']}
    require(len(expected) == len(receipt['inventory']), 'duplicate publication inventory role')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for member in tar:
            name = PurePosixPath(member.name)
            require(member.isfile() and not name.is_absolute() and '..' not in name.parts
                    and member.name not in files and member.name in expected,
                    'unsafe or unexpected publication member')
            value = tar.extractfile(member).read()
            row = expected[member.name]
            require(len(value) == row['bytes'] and digest(value) == row['sha256']
                    and member.mode == row['mode'], 'published member differs from inventory')
            files[member.name] = value
    require(set(files) == set(expected), 'publication lost a file')
    out.mkdir(parents=True)
    for name, value in files.items():
        path = out/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)
        path.chmod(expected[name]['mode'])
    for family in ('f5', 'sat'):
        for directory in ('extracted', 'asset-files'):
            make_read_only(out/(family+'-execution')/directory)
    return dict(status='RESTORED_EXACT_ORIGINAL_BYTES', members=len(files),
                archive_sha256=receipt['archive_sha256'], native_solvers_executed=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    create = sub.add_parser('pack')
    create.add_argument('--root', type=Path, required=True)
    create.add_argument('--diagnosis', type=Path, required=True)
    create.add_argument('--out', type=Path, required=True)
    extract = sub.add_parser('restore')
    extract.add_argument('--bundle', type=Path, required=True)
    extract.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = pack(args.root, args.diagnosis, args.out) if args.command == 'pack' else restore(args.bundle, args.out)
    print(json.dumps({key:value for key,value in result.items() if key != 'inventory'}, sort_keys=True))
