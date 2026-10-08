"""Retain terminal F4/F5 v1 controls and replay them without native execution.

Controller failures use publish_sat_runtime_failure_v3 instead. A native
timeout with complete controller gates is retained here with unknown math.
"""
import argparse
import gzip
import io
import json
from pathlib import Path, PurePosixPath
import tarfile

from audit_f5_runtime_v1 import audit
from identity import sha256, write_immutable
from oracle import require
from publish_static_sat_v3_control import collect, digest, inventory
from sat_runtime_execution_v3 import read


def publish(registration, execution, audit_file, output, expected_execution_sha256):
    registration, execution, output = map(Path, (registration, execution, output))
    require(not output.exists(), 'F5 publication exists; never overwrite evidence')
    spec = read(registration / 'execution.json')
    require(sha256(spec) == expected_execution_sha256,
            'F5 publication differs from externally frozen invocation')
    result = audit(execution, spec)
    require(result == read(audit_file), 'F5 independent receipt differs from fresh replay')
    files = collect(registration, 'registration')
    files.update(collect(execution, 'execution'))
    files['independent-audit.json'] = (Path(audit_file).read_bytes(), 0o444)
    files['publisher.py'] = (Path(__file__).read_bytes(), 0o444)
    output.mkdir(parents=True)
    archive = output / 'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for role, (data, mode) in sorted(files.items()):
                item = tarfile.TarInfo(role)
                item.size, item.mtime, item.mode = len(data), 0, mode
                tar.addfile(item, io.BytesIO(data))
    receipt = dict(schema_version=1, execution_sha256=expected_execution_sha256,
        archive_sha256=digest(archive.read_bytes()), archive_bytes=archive.stat().st_size,
        inventory=inventory(files), result=result, promotion_eligible=False, online_speedup=None)
    write_immutable(output / 'receipt.json', receipt)
    write_immutable(output / 'AUDIT.json', result)
    return result


def replay(bundle, output, expected_execution_sha256):
    bundle, output = map(Path, (bundle, output))
    receipt = read(bundle / 'receipt.json')
    require(receipt['execution_sha256'] == expected_execution_sha256,
            'F5 evidence differs from externally frozen invocation')
    data = (bundle / 'evidence.tar.gz').read_bytes()
    require(digest(data) == receipt['archive_sha256'] and len(data) == receipt['archive_bytes'],
            'F5 published archive changed')
    require(not output.exists(), 'F5 replay output exists')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            role = PurePosixPath(item.name)
            require(item.isfile() and not role.is_absolute() and '..' not in role.parts
                    and role.as_posix() == item.name and item.name not in files
                    and 0 <= item.mode <= 0o777 and not item.mode & 0o022,
                    'unsafe or duplicate F5 archive member')
            files[item.name] = (tar.extractfile(item).read(), item.mode)
    require(inventory(files) == receipt['inventory'], 'F5 published inventory changed')
    output.mkdir(parents=True)
    for role, (data, mode) in files.items():
        path = output / role
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(mode)
    spec = read(output / 'registration/execution.json')
    require(sha256(spec) == expected_execution_sha256, 'F5 archived invocation hash changed')
    result = audit(output / 'execution', spec)
    require(result == receipt['result'] == read(output / 'independent-audit.json')
            == read(bundle / 'AUDIT.json') and receipt['promotion_eligible'] is False
            and receipt['online_speedup'] is None,
            'F5 transported source/mathematical assessment differs or claims promotion')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    packing = sub.add_parser('publish')
    packing.add_argument('--registration', type=Path, required=True)
    packing.add_argument('--execution', type=Path, required=True)
    packing.add_argument('--audit', type=Path, required=True)
    transport = sub.add_parser('replay')
    transport.add_argument('--bundle', type=Path, required=True)
    for command in (packing, transport):
        command.add_argument('--out', type=Path, required=True)
        command.add_argument('--expected-execution-sha256', required=True)
    args = parser.parse_args()
    result = (publish(args.registration, args.execution, args.audit, args.out,
                      args.expected_execution_sha256) if args.command == 'publish'
              else replay(args.bundle, args.out, args.expected_execution_sha256))
    print(json.dumps(result, sort_keys=True))
