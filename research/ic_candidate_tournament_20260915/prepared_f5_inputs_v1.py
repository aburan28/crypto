"""Retain a newly built prepared n17 worker and its complete Rust sources.

Assets retain code and mathematical state; the Rust sources include disclosed
unit controls excluded from the release worker. The historical preparation
certificates remain invocation data. Staging checks the controlled build and
dependency trees; it neither compiles nor executes a solver. A separately
sealed invocation imports its certificate and disclosed public point.
"""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import platform
import tarfile

import generic_build
from identity import sha256
from oracle import require
from static_sat_assets_v3 import freeze_assets

WORKER = 'examples/ic_tournament_worker.rs'
MATHEMATICS = ('research/ic_candidate_tournament_20260915/goal_20260924/'
               'prepared-f5-runtime-v1/mathematics.json')
WORKER_SOURCE_SHA256 = '3e82ceaa7e1da6ef295b161d73c6f64472ef901a85b46b523977713431072816'
MATHEMATICS_SHA256 = '54e7f33679c9cc86d09fc883c9e2b32710540b0f2dba8fe30602e5b5d8cce325'
ROLES = {'bin/worker', 'build/build-record.json', 'build/build-policy.json',
         'build/build-exit.json', 'build/build.log', 'rust/source-manifest.json',
         'rust/root-source.tar.gz', 'rust/dependency-source.tar.gz'}
BUILD_ARGUMENTS = ['cargo', 'build', '--locked', '--offline', '--release',
                   '--no-default-features', '--example', 'ic_tournament_worker']
PLATFORMS = {
    ('aarch64', 'macos'): ('Darwin', 'arm64', 'physical-macos-arm64'),
    ('x86_64', 'linux'): ('Linux', 'x86_64', 'development-linux-x86_64'),
}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def archive_members(data, expected):
    """Check every retained file without extraction or executing build scripts."""
    observed = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        for item in archive:
            name = PurePosixPath(item.name)
            require(item.isfile() and item.name not in observed and item.name
                    and name.as_posix() == item.name and not name.is_absolute()
                    and '..' not in name.parts, 'unsafe or duplicate prepared Rust source member')
            observed[item.name] = digest(archive.extractfile(item).read())
    require(observed == expected, 'prepared Rust source bytes differ from build manifest')


def dependency_members(source):
    expected = {}
    for item in source['dependencies']:
        require(item['source'].startswith('registry+') and item['files'],
                'prepared worker dependency is not a retained registry source')
        for name, value in item['files'].items():
            role = item['package']+'-'+item['version']+'/'+name
            require(role not in expected, 'ambiguous prepared dependency source role')
            expected[role] = value
    return expected


def native_admission(files, *, check_host=False):
    require(set(files) == ROLES, 'prepared worker assets have missing or extra roles')
    record = json.loads(files['build/build-record.json'])
    source = json.loads(files['rust/source-manifest.json'])
    identity = generic_build.verify_build_record(record, source)
    policy = record['build']
    require(policy['arguments'] == BUILD_ARGUMENTS
            and policy['flags'] == dict(rustflags='', incremental=False, features=[], profile='release')
            and json.loads(files['build/build-policy.json']) == policy
            and json.loads(files['build/build-exit.json']) == {'exit_code':0}
            and digest(files['bin/worker']) == record['worker_sha256']
            and record['builder_sha256'] == generic_build.digest(generic_build.__file__),
            'prepared worker executable, controlled recipe or builder differs')
    require(source['root_files'].get(WORKER) == WORKER_SOURCE_SHA256
            and source['root_files'].get(MATHEMATICS) == MATHEMATICS_SHA256
            and not any(name.endswith('/f5-preparation.json') or name.endswith('/sat-preparation.json')
                        for name in source['root_files']),
            'prepared worker source lacks the admitted mode or hashes preparation history')
    archive_members(files['rust/root-source.tar.gz'], source['root_files'])
    archive_members(files['rust/dependency-source.tar.gz'], dependency_members(source))
    key = (policy['target_arch'], policy['target_os'])
    require(key in PLATFORMS, 'prepared worker platform needs a separately validated adapter')
    system, machine, host_class = PLATFORMS[key]
    if check_host:
        require((platform.system(), platform.machine()) == (system, machine),
                'prepared native worker requires its own build platform')
    native = dict(source_manifest_sha256=sha256(source), worker_sha256=record['worker_sha256'],
        build_identity=identity, build_record_sha256=sha256(record),
        platform=dict(system=system, machine=machine), host_class=host_class,
        scope='controlled local build and retained sources; no remote attestation or speedup')
    return record, source, native


def source_archive(files):
    raw = io.BytesIO()
    with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='') as compressed:
        with tarfile.open(fileobj=compressed, mode='w|') as archive:
            for name, data in sorted(files.items()):
                item = tarfile.TarInfo(name)
                item.size, item.mode, item.mtime = len(data), 0o444, 0
                archive.addfile(item, io.BytesIO(data))
    return raw.getvalue()


def stage(build, root, output):
    """Seal already-built bytes once; verify the complete dependency closure."""
    build, root, output = Path(build), Path(root).resolve(), Path(output)
    require(not output.exists(), 'prepared worker asset snapshot already exists')
    source = json.loads((build/'source-manifest.json').read_text())
    env = generic_build.controlled_environment(root)
    metadata = json.loads(generic_build.command(
        ['cargo', 'metadata', '--locked', '--offline', '--format-version', '1'], root, env))
    require(generic_build.source_manifest(root, metadata) == source,
            'source or dependency trees changed after prepared worker build')
    packages = {(p['name'],p['version'],p['source']):Path(p['manifest_path']).parent
                for p in metadata['packages'] if Path(p['manifest_path']).parent != root}
    dependencies = {}
    for item in source['dependencies']:
        directory = packages[(item['package'],item['version'],item['source'])]
        for name in item['files']:
            dependencies[item['package']+'-'+item['version']+'/'+name] = (directory/name).read_bytes()
    files = {'bin/worker':(build/'worker').read_bytes(),
        'rust/source-manifest.json':(build/'source-manifest.json').read_bytes(),
        'rust/root-source.tar.gz':(build/'root-source.tar.gz').read_bytes(),
        'rust/dependency-source.tar.gz':source_archive(dependencies)}
    for name in ('build-record.json','build-policy.json','build-exit.json','build.log'):
        files['build/'+name] = (build/name).read_bytes()
    record, _, native = native_admission(files, check_host=True)
    manifest, seal = freeze_assets(files, {'bin/worker'}, output)
    return dict(native=native, manifest_sha256=sha256(manifest), seal=seal,
                build_sha256=record['build_sha256'], solver_executed=False, fresh_targets_generated=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('build','root','out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(stage(args.build,args.root,args.out), sort_keys=True))


if __name__ == '__main__':
    main()
