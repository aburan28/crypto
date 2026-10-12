#!/usr/bin/env python3
"""Build the accepted pairinv source for a separately labelled macOS arm."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

from campaign_rules import IC_SOURCE
from identity import sha256
from oracle import require
from producer.prepare import verify_source
from tournament import write

TARGET = 'aarch64-apple-darwin'
OVERRIDES = ('build.target="aarch64-apple-darwin"',
             'build.rustflags=[]')


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dependency_manifest(metadata, root):
    packages = []
    for package in metadata['packages']:
        directory = Path(package['manifest_path']).parent
        if directory == root:
            continue
        require(package['source'] is not None
                and package['source'].startswith('registry+'),
                'local/git dependency needs separate build admission')
        files = sorted(path for path in directory.rglob('*') if path.is_file()
                       and path.name not in {'.cargo-ok', '.cargo_vcs_info.json'})
        require(not any(path.is_symlink() for path in files),
                'symlinked dependency source')
        packages.append(dict(name=package['name'], version=package['version'],
                             source=package['source'], files={
                                 path.relative_to(directory).as_posix(): digest(path)
                                 for path in files}))
    packages.sort(key=lambda row: (row['name'], row['version'], row['source']))
    return dict(schema_version=1, packages=packages)


def build(prepared, out):
    require(platform.system() == 'Darwin' and platform.machine() == 'arm64',
            'local reference build requires physical macOS arm64')
    prepared, out = Path(prepared).resolve(), Path(out).resolve()
    receipt = json.loads((prepared/'preparation.json').read_text())
    require(receipt['reference'] == 'pairinv'
            and receipt['instrumented'] is True
            and receipt['source_manifest_sha256'] == IC_SOURCE,
            'source is not the accepted instrumented pairinv reference')
    source = prepared/'source'
    manifest = verify_source(prepared, IC_SOURCE)
    require((source/'.cargo/config.toml').is_file()
            and 'x86_64-unknown-linux-musl'
                in (source/'.cargo/config.toml').read_text(),
            'expected archived Linux build policy changed')
    require(not out.exists(), 'local pairinv build output already exists')
    out.mkdir(parents=True)
    environment = {key: os.environ[key] for key in
                   ('PATH', 'HOME', 'CARGO_HOME', 'RUSTUP_HOME', 'TMPDIR')
                   if key in os.environ}
    environment.update(IC_SOURCE_MANIFEST_SHA256=IC_SOURCE,
                       CARGO_TARGET_DIR=str(out/'target'),
                       CARGO_INCREMENTAL='0', RUSTFLAGS='',
                       LC_ALL='C', CARGO_TERM_COLOR='never')
    prefix = ['cargo', '--config', OVERRIDES[0], '--config', OVERRIDES[1]]
    def output(command):
        process = subprocess.run(command, cwd=source, env=environment,
                                 capture_output=True, text=True, check=False)
        require(process.returncode == 0,
                f'local pairinv preflight failed: {command}: {process.stderr}')
        return process.stdout.strip()
    rustc = output(['rustc', '-Vv'])
    cargo = output(['cargo', '-Vv'])
    metadata = json.loads(output(prefix + [
        'metadata', '--locked', '--offline', '--format-version', '1']))
    deps = dependency_manifest(metadata, source)
    write(out/'dependency-manifest.json', deps, exclusive=True)
    arguments = prefix + ['build', '--locked', '--offline', '--release',
                          '--target', TARGET, '--example',
                          'ic_tournament_worker', '--jobs', '2']
    policy = dict(schema_version=1, accepted_source_manifest_sha256=IC_SOURCE,
                  prepared_source_receipt_sha256=sha256(receipt),
                  compiler=rustc, cargo=cargo, target=TARGET,
                  cargo_config_sha256=manifest['.cargo/config.toml'],
                  cargo_lock_sha256=manifest['Cargo.lock'],
                  dependency_manifest_sha256=sha256(deps),
                  cargo_overrides=list(OVERRIDES), arguments=arguments,
                  flags=dict(rustflags='', incremental=False,
                             release=True, features='default'),
                  scope='local macOS build of accepted source; not Linux instruction qualification')
    write(out/'build-policy.json', policy, exclusive=True)
    with (out/'build.log').open('x') as log:
        process = subprocess.run(arguments, cwd=source, env=environment,
                                 stdout=log, stderr=subprocess.STDOUT,
                                 check=False)
    write(out/'build-exit.json', dict(exit_code=process.returncode),
          exclusive=True)
    require(process.returncode == 0, 'local pairinv build failed; log retained')
    require(verify_source(prepared, IC_SOURCE) == manifest,
            'accepted source changed during local build')
    require(dependency_manifest(metadata, source) == deps,
            'registry dependency contents changed during local build')
    executable = out/'target'/TARGET/'release/examples/ic_tournament_worker'
    require(executable.is_file(), 'local pairinv worker was not produced')
    shutil.copy2(executable, out/'worker')
    build_record = dict(schema_version=1, accepted_source_manifest_sha256=IC_SOURCE,
                        source_manifest=manifest,
                        prepared_source_receipt=receipt,
                        build_policy=policy,
                        build_policy_sha256=sha256(policy),
                        worker_sha256=digest(out/'worker'),
                        build_log_sha256=digest(out/'build.log'),
                        dependency_manifest_sha256=sha256(deps),
                        builder_sha256=digest(__file__))
    write(out/'build-record.json', build_record, exclusive=True)
    return build_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    record = build(args.prepared, args.out)
    print(json.dumps({key: record[key] for key in
                      ('accepted_source_manifest_sha256', 'worker_sha256',
                       'build_policy_sha256')}), flush=True)


if __name__ == '__main__':
    main()
