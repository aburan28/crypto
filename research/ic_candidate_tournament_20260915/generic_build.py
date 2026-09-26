"""Controlled, source-bound build of the bounded generic research worker.

This is an auditable local build receipt, not remote attestation or a hermetic
toolchain claim. Reject Cargo configuration overrides, hash actual dependency
trees as well as the lockfile, and check inputs again after compilation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

from identity import sha256, write_immutable
from oracle import require

ROOT = Path(__file__).resolve().parents[2]
WORKER = 'examples/ic_tournament_worker.rs'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def controlled_environment(root):
    require(os.name == 'posix', 'builder requires POSIX host')
    env = {key: os.environ[key] for key in ('PATH', 'HOME', 'CARGO_HOME', 'RUSTUP_HOME', 'TMPDIR')
           if key in os.environ}
    cargo_home = Path(env.get('CARGO_HOME', str(Path.home() / '.cargo')))
    configs = [parent / '.cargo' / name for parent in (root, *root.parents)
               for name in ('config', 'config.toml')]
    configs += [cargo_home / name for name in ('config', 'config.toml')]
    require(not any(path.exists() for path in configs), 'Cargo config override needs a separate build adapter')
    env.update(LC_ALL='C', CARGO_TERM_COLOR='never', CARGO_INCREMENTAL='0', RUSTFLAGS='')
    for name in ('cc', 'ar'):
        tool = shutil.which(name, path=env['PATH'])
        require(tool is not None, 'missing native build tool ' + name)
        env[name.upper()] = tool
    return env


def command(args, root, env):
    return subprocess.check_output(args, cwd=root, env=env, text=True, stderr=subprocess.STDOUT).strip()


def source_manifest(root, metadata):
    require(not (root / 'build.rs').exists(), 'root build script needs explicit admission')
    names = [Path('Cargo.toml'), Path('Cargo.lock'), Path(WORKER)]
    names += sorted(path.relative_to(root) for path in (root / 'src').rglob('*') if path.is_file())
    # Library modules also include files outside src (for example the pinned
    # boundary calibration). Hash literal include inputs, including test-only
    # ones, so the manifest is a conservative superset of compilation inputs.
    included = set()
    for name in names:
        if name.suffix == '.rs':
            text = (root / name).read_text()
            for relative in re.findall(r'include(?:_str|_bytes)?!\(\s*"([^"]+)"', text):
                path = (root / name).parent / relative
                require(path.resolve().is_relative_to(root), 'include input escapes source root')
                included.add(path.resolve().relative_to(root))
    names = sorted(set(names) | included)
    require(not any((root / path).is_symlink() for path in names), 'symlinked root build input')
    root_files = {path.as_posix(): digest(root / path) for path in names}
    dependencies = []
    for package in metadata['packages']:
        directory = Path(package['manifest_path']).parent
        if directory == root:
            continue
        require(package['source'] is not None and package['source'].startswith('registry+'),
                'path/git dependency needs explicit build admission')
        # Include native headers/build scripts/data, not only Rust modules.
        files = sorted(path for path in directory.rglob('*') if path.is_file()
                       and path.name not in {'.cargo-ok', '.cargo_vcs_info.json'})
        require(not any(path.is_symlink() for path in files), 'symlinked dependency input')
        dependencies.append(dict(package=package['name'], version=package['version'],
                                 source=package['source'],
                                 files={path.relative_to(directory).as_posix(): digest(path) for path in files}))
    dependencies.sort(key=lambda p: (p['package'], p['version']))
    return dict(schema_version=1, root_files=root_files, dependencies=dependencies)


def verify_build_record(record, source):
    require(record.get('schema_version') == 1 and source.get('schema_version') == 1,
            'unknown build record schema')
    require(record['source_manifest_sha256'] == sha256(source), 'source manifest digest mismatch')
    require(record['build_sha256'] == sha256(record['build']), 'build policy digest mismatch')
    require(record['build']['source_manifest_sha256'] == record['source_manifest_sha256'],
            'build policy refers to different source')
    require(record['identity'] == dict(schema_version=1,
                                       source_manifest_sha256=record['source_manifest_sha256'],
                                       build_sha256=record['build_sha256'],
                                       target_arch=record['build']['target_arch'],
                                       target_os=record['build']['target_os']), 'embedded build identity mismatch')
    for value in (record['worker_sha256'], record['builder_sha256']):
        require(type(value) is str and len(value) == 64
                and all(ch in '0123456789abcdef' for ch in value), 'invalid build artifact digest')
    return record['identity']


def verify_binding(report, record, source, *, executable):
    identity = verify_build_record(record, source)
    require(report.get('generic_build') == identity, 'worker reports a different build identity')
    require(digest(executable) == record['worker_sha256'], 'worker executable digest mismatch')
    return dict(schema_version=1, source_manifest_sha256=record['source_manifest_sha256'],
                build_sha256=record['build_sha256'], worker_sha256=record['worker_sha256'],
                scope='controlled local build and retained executable, not remote attestation')


def build(root, out):
    root, out = root.resolve(), out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = controlled_environment(root)
    rustc = command(['rustc', '-Vv'], root, env)
    cargo = command(['cargo', '-Vv'], root, env)
    cfg = command(['rustc', '--print', 'cfg'], root, env).splitlines()
    target = {name: next(json.loads(line.split('=', 1)[1]) for line in cfg if line.startswith(name + '='))
              for name in ('target_arch', 'target_os')}
    metadata = json.loads(command(['cargo', 'metadata', '--locked', '--offline', '--format-version', '1'], root, env))
    source = source_manifest(root, metadata)
    write_immutable(out / 'source-manifest.json', source)
    with tarfile.open(out / 'root-source.tar.gz', 'x:gz') as archive:
        for name in source['root_files']:
            archive.add(root / name, arcname=name, recursive=False)
    native_tools = {}
    for name in ('CC', 'AR'):
        cp = subprocess.run([env[name], '--version'], cwd=root, env=env,
                            capture_output=True, text=True, check=False)
        native_tools[name] = dict(binary_sha256=digest(env[name]), exit_code=cp.returncode,
                                  version_stdout=cp.stdout, version_stderr=cp.stderr)
    args = ['cargo', 'build', '--locked', '--offline', '--release', '--no-default-features',
            '--example', 'ic_tournament_worker']
    policy = dict(schema_version=1, source_manifest_sha256=sha256(source),
                  rustc=rustc, cargo=cargo, native_tools=native_tools, arguments=args,
                  **target, flags=dict(rustflags='', incremental=False, features=[], profile='release'),
                  scope='local toolchain build; OS and system libraries are not hermetic')
    build_hash = sha256(policy)
    env.update(IC_GENERIC_SOURCE_MANIFEST_SHA256=sha256(source), IC_GENERIC_BUILD_SHA256=build_hash)
    write_immutable(out / 'build-policy.json', policy)
    with (out / 'build.log').open('x') as log:
        cp = subprocess.run(args, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    write_immutable(out / 'build-exit.json', dict(exit_code=cp.returncode))
    require(cp.returncode == 0, 'worker build failed; log retained')
    require(source_manifest(root, metadata) == source, 'source changed during build')
    executable = Path(metadata['target_directory']) / 'release/examples/ic_tournament_worker'
    shutil.copy2(executable, out / 'worker')
    identity = json.loads(command([str(out / 'worker'), '--build-identity'], root, env))
    record = dict(schema_version=1, source_manifest_sha256=sha256(source), build_sha256=build_hash,
                  build=policy, identity=identity, worker_sha256=digest(out / 'worker'),
                  builder_sha256=digest(__file__))
    verify_build_record(record, source)
    write_immutable(out / 'build-record.json', record)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = build(args.root, args.out)
    print(json.dumps({key: result[key] for key in ('build_sha256', 'source_manifest_sha256', 'worker_sha256')}))
