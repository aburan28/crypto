"""Transport failed/partial registered executions without scientific admission.

This preserves bytes and operational outcomes. It does not audit relations,
rank, scalar recovery, native completion, or a successful IC timing interval.
"""
import argparse
import gzip
import io
import json
from pathlib import Path, PurePosixPath
import tarfile

from identity import sha256, write_immutable
from oracle import InvalidEvidence, require
from publish_static_sat_v3_control import collect, digest, inventory
from sat_runtime_bundle import verified_files
from sat_runtime_execution_v3 import audit_execution, read, validate_spec
from static_sat_assets_v3 import verified_assets


def assess(registration, execution, expected_execution_sha256):
    registration, execution = map(Path, (registration, execution))
    spec = read(registration/'execution.json')
    require(sha256(spec) == expected_execution_sha256,
            'failure evidence differs from externally frozen invocation')
    validate_spec(spec)
    seal = dict(schema_version=3, execution_sha256=sha256(spec),
                binding=spec['binding'], registration_stage='before-execution')
    require(read(registration/'registration-seal.json') == seal
            and read(execution/'execution.json') == spec
            and read(execution/'registration-seal.json') == seal,
            'failure evidence registration or retained invocation changed')
    verified_files(registration/'runtime', spec['runtime_manifest'], spec['runtime_seal'])
    if spec.get('asset_seal') is not None:
        verified_assets(registration/'assets', spec['asset_manifest'], spec['asset_seal'])
    process = read(execution/'process.json') if (execution/'process.json').is_file() else None
    if process is not None:
        require(process['execution_sha256'] == sha256(spec)
                and process['binding'] == spec['binding']
                and type(process['exit_code']) is int
                and type(process['timed_out']) is bool
                and (process['timed_out'] or process['exit_code'] != 0)
                and process['python_flags'] == ['-I', '-S', '-B']
                and type(process['process_wall_ns']) is int
                and process['process_wall_ns'] >= 0
                and process['watchdog_seconds'] == spec['runtime_watchdog_seconds']
                and process['stdout_sha256'] == digest((execution/'stdout.txt').read_bytes())
                and process['stderr_sha256'] == digest((execution/'stderr.txt').read_bytes())
                and process['promotion_eligible'] is False
                and process['online_speedup'] is None,
                'failure process receipt changed or claims a successful execution')
    source_audit = None
    source_error = None
    try:
        source_audit = audit_execution(execution, spec)
    except (InvalidEvidence, OSError, ValueError, KeyError, TypeError) as error:
        # Exception messages can contain extraction paths. Preserve those in
        # raw stderr; this relocation-stable assessment records only the type.
        source_error = type(error).__name__
    require(source_audit is None or source_audit['entrypoint_succeeded'] is False,
            'successful source execution belongs to the scientific audit path')
    mathematical_seal = spec['arguments'].get('seal', {}) if isinstance(spec['arguments'], dict) else {}
    return dict(schema_version=3, status=(
        'PROCESS_RECORD_MISSING' if process is None else
        'CONTROLLER_TIMEOUT' if process['timed_out'] else 'EXECUTION_FAILURE'),
        execution_sha256=sha256(spec), binding=spec['binding'],
        candidate_id=mathematical_seal.get('candidate_id'),
        workload_id=mathematical_seal.get('workload_id'), run_id=mathematical_seal.get('run_id'),
        process=process, python_source_audit=source_audit, source_audit_error_type=source_error,
        claim_scope='raw failed/partial execution transport; no mathematical or timing admission',
        complete_ic_source_bound=False, mathematical_audit=None,
        verified_target_count=None, final_rank=None, online_wall_ns=None,
        online_phases_ns=None, online_speedup=None, promotion_eligible=False,
        headline_online_admissible=False)


def publish(registration, execution, output, expected_execution_sha256):
    output = Path(output)
    require(not output.exists(), 'failure publication exists; never overwrite evidence')
    result = assess(registration, execution, expected_execution_sha256)
    files = collect(registration, 'registration')
    files.update(collect(execution, 'execution'))
    files['publisher.py'] = (Path(__file__).read_bytes(), 0o444)
    output.mkdir(parents=True)
    archive = output/'evidence.tar.gz'
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as tar:
            for role, (data, mode) in sorted(files.items()):
                item = tarfile.TarInfo(role)
                item.size, item.mtime, item.mode = len(data), 0, mode
                tar.addfile(item, io.BytesIO(data))
    receipt = dict(schema_version=3, execution_sha256=expected_execution_sha256,
                   archive_sha256=digest(archive.read_bytes()), archive_bytes=archive.stat().st_size,
                   inventory=inventory(files), result=result,
                   promotion_eligible=False, online_speedup=None)
    write_immutable(output/'receipt.json', receipt)
    write_immutable(output/'FAILURE.json', result)
    return result


def replay(bundle, output, expected_execution_sha256):
    bundle, output = map(Path, (bundle, output))
    receipt = read(bundle/'receipt.json')
    require(receipt['execution_sha256'] == expected_execution_sha256,
            'failure evidence differs from externally frozen invocation')
    data = (bundle/'evidence.tar.gz').read_bytes()
    require(digest(data) == receipt['archive_sha256'] and len(data) == receipt['archive_bytes'],
            'published failure archive changed')
    require(not output.exists(), 'failure replay output exists')
    files = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
        for item in tar:
            role = PurePosixPath(item.name)
            require(item.isfile() and not role.is_absolute() and '..' not in role.parts
                    and role.as_posix() == item.name and item.name not in files
                    and 0 <= item.mode <= 0o777 and not item.mode & 0o022,
                    'unsafe or duplicate failure archive member')
            files[item.name] = (tar.extractfile(item).read(), item.mode)
    require(inventory(files) == receipt['inventory'], 'published failure inventory changed')
    output.mkdir(parents=True)
    for role, (data, mode) in files.items():
        path = output/role
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(mode)
    result = assess(output/'registration', output/'execution', expected_execution_sha256)
    require(result == receipt['result'] == read(bundle/'FAILURE.json')
            and receipt['promotion_eligible'] is False and receipt['online_speedup'] is None,
            'transported failure assessment differs or claims promotion')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    packing = sub.add_parser('publish')
    packing.add_argument('--registration', type=Path, required=True)
    packing.add_argument('--execution', type=Path, required=True)
    transport = sub.add_parser('replay')
    transport.add_argument('--bundle', type=Path, required=True)
    for command in (packing, transport):
        command.add_argument('--out', type=Path, required=True)
        command.add_argument('--expected-execution-sha256', required=True)
    args = parser.parse_args()
    result = (publish(args.registration, args.execution, args.out, args.expected_execution_sha256)
              if args.command == 'publish' else replay(args.bundle, args.out, args.expected_execution_sha256))
    print(json.dumps(result, sort_keys=True))
