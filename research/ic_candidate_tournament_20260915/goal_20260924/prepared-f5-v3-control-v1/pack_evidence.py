#!/usr/bin/env python3
"""Capture/restore this consumed toy control; never execute a native solver."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import stat
import tarfile

EXECUTION = '5c0a1473b7d2637a7e97666db68262d038c628f937f3b0ea0338fd83f949f98d'


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_control(root):
    spec = json.loads((root/'execution/execution.json').read_text())
    canonical = json.dumps(spec, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()
    require(digest(canonical) == EXECUTION, 'not the externally registered control')
    claim = json.loads((root/'execution/execution-claim.json').read_text())
    require(claim['execution_sha256'] == EXECUTION
            and claim['status'] == 'CONSUMED_BEFORE_LAUNCH', 'original claim missing')
    admission = json.loads((root/'audit/admission.json').read_text())
    transport = json.loads((root/'audit/transport.json').read_text())
    require(admission['status'] == 'ADMITTED_COMPLETE_PREPARED_F5_CONTROL'
            and admission['scalar_verified'] and admission['target_attempt_count'] == 3
            and transport['status'] == 'PASS_FROZEN_PREPARED_TRANSPORT'
            and transport['execution_sha256'] == EXECUTION
            and transport['native_solvers_executed'] == 0, 'terminal record differs')
    require(sum(admission['online_phases_ns'].values()) == admission['online_wall_ns'],
            'exclusive phase sum differs')


def pack(root, out):
    require(not out.exists(), 'publication already exists')
    verify_control(root)
    files = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'symlink in original evidence')
        if path.is_file():
            files[path.relative_to(root).as_posix()] = path
    inventory = []
    out.mkdir(parents=True)
    archive = out/'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as zipped:
        with tarfile.open(fileobj=zipped, mode='w|') as tar:
            for name, path in sorted(files.items()):
                data = path.read_bytes()
                mode = stat.S_IMODE(path.stat().st_mode)
                require(not mode & ~0o777, 'special mode in original evidence')
                member = tarfile.TarInfo(name)
                member.size, member.mode, member.mtime = len(data), mode, 0
                tar.addfile(member, io.BytesIO(data))
                inventory.append(dict(role=name, bytes=len(data), sha256=digest(data), mode=mode))
    require(all(digest(files[r['role']].read_bytes()) == r['sha256'] for r in inventory),
            'evidence changed during capture; retain rejected archive')
    receipt = dict(schema_version=1, status='RETAINED_ORIGINAL_COMPLETE_DISCLOSED_CONTROL',
                   execution_sha256=EXECUTION, archive_sha256=digest(archive.read_bytes()),
                   archive_bytes=archive.stat().st_size, inventory=inventory,
                   native_solvers_executed_by_publication=0, promotion_eligible=False,
                   fresh_paired_qualification=False, online_speedup=None)
    with (out/'receipt.json').open('x') as stream:
        json.dump(receipt, stream, sort_keys=True, separators=(',', ':'))
        stream.write('\n')
    return {k: v for k, v in receipt.items() if k != 'inventory'}


def restore(bundle, out):
    require(not out.exists(), 'restoration already exists')
    receipt = json.loads((bundle/'receipt.json').read_text())
    data = (bundle/'evidence.tar.gz').read_bytes()
    require(receipt['execution_sha256'] == EXECUTION and len(data) == receipt['archive_bytes']
            and digest(data) == receipt['archive_sha256'], 'archive seal differs')
    expected = {r['role']: r for r in receipt['inventory']}
    require(len(expected) == len(receipt['inventory']), 'duplicate inventory member')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for member in tar:
            name = PurePosixPath(member.name)
            require(member.isfile() and not name.is_absolute() and '..' not in name.parts
                    and member.name not in files and member.name in expected, 'unsafe archive member')
            value, row = tar.extractfile(member).read(), expected[member.name]
            require(len(value) == row['bytes'] and digest(value) == row['sha256']
                    and member.mode == row['mode'] and not member.mode & ~0o777,
                    'member bytes or mode differ')
            files[member.name] = value
    require(set(files) == set(expected), 'archive inventory incomplete')
    out.mkdir(parents=True)
    for name, value in files.items():
        path = out/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(value)
        path.chmod(expected[name]['mode'])
    verify_control(out)
    return dict(status='RESTORED_EXACT_ORIGINAL_BYTES', members=len(files),
                archive_sha256=receipt['archive_sha256'], native_solvers_executed=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    create = sub.add_parser('pack')
    create.add_argument('--root', type=Path, required=True)
    create.add_argument('--out', type=Path, required=True)
    extract = sub.add_parser('restore')
    extract.add_argument('--bundle', type=Path, required=True)
    extract.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = pack(args.root, args.out) if args.command == 'pack' else restore(args.bundle, args.out)
    print(json.dumps(result, sort_keys=True))
