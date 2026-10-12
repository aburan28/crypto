"""Build a diagnostic overlay on exact historical Rust sources; run once.

This never writes to the historical worktree or executes a historical worker.
All failures and partial output remain in the new output directory.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[2]
sys.path.insert(0, str(PACKAGE))
from f5_job_schema_control import audit
from generic_build import command, controlled_environment, digest, source_manifest
from identity import sha256, write_immutable
from oracle import require


def supervised(args, cwd, env, out, name, limit, stdin=None):
    with (out/(name+'.stdout')).open('xb') as stdout, (out/(name+'.stderr')).open('xb') as stderr:
        child = subprocess.Popen(args, cwd=cwd, env=env, stdout=stdout, stderr=stderr, stdin=stdin, start_new_session=True)
        timed_out = False
        try:
            child.wait(timeout=limit)
        except subprocess.TimeoutExpired:
            timed_out = True
        finally:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
    write_immutable(out/(name+'-process.json'), dict(argv=args, exit_code=child.returncode,
        timeout=timed_out, watchdog_seconds=limit, owned_group_cleanup=True,
        stdout_sha256=digest(out/(name+'.stdout')), stderr_sha256=digest(out/(name+'.stderr'))))
    require(not timed_out and child.returncode == 0, name+' failed; terminal/partial output retained')


def tool_receipt(path, cwd, env):
    # A version flag is optional metadata, not the build command. In particular
    # BSD ar does not implement --version. Retain its actual failure verbatim.
    args = [path, '--version']
    cp = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True,
                        check=False, timeout=10)
    return dict(path=path, binary_sha256=digest(path), argv=args, exit_code=cp.returncode,
                stdout=cp.stdout, stderr=cp.stderr, version_probe_succeeded=cp.returncode == 0)


def run(original, out, protocol_file):
    require(not out.exists(), 'control output exists; do not retry or overwrite')
    protocol = json.loads(protocol_file.read_text())
    require(protocol['schema_version'] == 1, 'unknown Job schema control registration')
    require(not (HERE/'TERMINAL.json').exists(), 'Job schema registration closed; use retained replay')
    require(protocol['executions'] == 1 and protocol['build_seconds'] == 1800
            and protocol['native_seconds'] == 30, 'registered schema control bounds changed')
    out.mkdir(parents=True)
    analysis_files = [Path(__file__), HERE/'export.rs', HERE/'inputs.json', PACKAGE/'f5_job_schema_control.py', HERE/'kernel_append.rs',
                      PACKAGE/'generic_build.py', PACKAGE/'identity.py', PACKAGE/'oracle.py',
                      protocol_file, HERE/'PROTOCOL.md']
    frozen = {str(path): digest(path) for path in analysis_files}
    require(digest(HERE/'inputs.json') == protocol['input_sha256'], 'disclosed inputs changed')
    env = controlled_environment(original)
    env.update(RAYON_NUM_THREADS='1', CARGO_BUILD_JOBS='1', PYTHONDONTWRITEBYTECODE='1')
    metadata = json.loads(command(['cargo', 'metadata', '--locked', '--offline', '--format-version', '1'], original, env))
    source = source_manifest(original, metadata)
    require(sha256(source) == protocol['historical_source_manifest_sha256'],
            'historical source or dependency tree differs; do not substitute current main')
    root = out/'source'
    root.mkdir()
    for name in source['root_files']:
        destination = root/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original/name, destination)
    overlay = root/'examples/f5_job_schema_control.rs'
    shutil.copyfile(HERE/'export.rs', overlay)
    kernel_name = 'examples/ic_tournament_worker.rs'
    kernel = root/kernel_name
    base_kernel = kernel.read_bytes()
    append = (HERE/'kernel_append.rs').read_bytes()
    require(digest(HERE/'kernel_append.rs') == protocol['kernel_append_sha256'], 'diagnostic append changed')
    kernel.write_bytes(base_kernel+append)
    shutil.copyfile(original/kernel_name, out/'kernel-original.rs')
    expected_source = json.loads(json.dumps(source))
    expected_source['root_files'][kernel_name] = digest(kernel)
    write_immutable(out/'base-source-manifest.json', source)
    compiled_source = dict(expected_source, diagnostic_overlay={
        'examples/f5_job_schema_control.rs':digest(overlay),
        kernel_name:dict(base_sha256=source['root_files'][kernel_name],
                        append_sha256=digest(HERE/'kernel_append.rs'), compiled_sha256=digest(kernel))})
    write_immutable(out/'source-manifest.json', compiled_source)
    with tarfile.open(out/'root-source.tar.gz', 'x:gz') as archive:
        for name in [*source['root_files'], 'examples/f5_job_schema_control.rs']:
            archive.add(root/name, arcname=name, recursive=False)
    with tarfile.open(out/'dependency-source.tar.gz', 'x:gz') as archive:
        directories = {Path(item['manifest_path']).parent.name: Path(item['manifest_path']).parent
                       for item in metadata['packages'] if item['source'] is not None}
        for item in source['dependencies']:
            directory = directories[item['package']+'-'+item['version']]
            for name in item['files']:
                archive.add(directory/name, arcname=item['package']+'-'+item['version']+'/'+name, recursive=False)
    env['CARGO_TARGET_DIR'] = str(out/'target')
    args = ['cargo', 'build', '--locked', '--offline', '--release', '--no-default-features', '--example', 'f5_job_schema_control']
    record = dict(schema_version=1, historical_source_manifest_sha256=sha256(source),
        compiled_source_manifest_sha256=sha256(compiled_source), analysis_sources=frozen,
        protocol_sha256=digest(protocol_file), protocol_document_sha256=digest(HERE/'PROTOCOL.md'),
        rustc=command(['rustc', '-Vv'], root, env), cargo=command(['cargo', '-Vv'], root, env),
        cc=tool_receipt(env['CC'], root, env), ar=tool_receipt(env['AR'], root, env),
        rustc_launcher=tool_receipt(shutil.which('rustc', path=env['PATH']), root, env),
        cargo_launcher=tool_receipt(shutil.which('cargo', path=env['PATH']), root, env),
        installed_rustc_sha256=digest(Path(command(['rustc', '--print', 'sysroot'], root, env))/'bin/rustc'),
        interpreter=dict(executable=sys.executable, sha256=digest(sys.executable), version=sys.version),
        environment=env, build_argv=args, input_sha256=protocol['input_sha256'],
        limits=dict(build_seconds=protocol['build_seconds'], native_seconds=protocol['native_seconds']),
        scope='local source-bound correctness diagnostic; not remote attestation or complete DLP admission')
    write_immutable(out/'preexecution.json', record)
    supervised(args, root, env, out, 'build', protocol['build_seconds'])
    staged_metadata = json.loads(command(['cargo', 'metadata', '--locked', '--offline', '--format-version', '1'], root, env))
    require(source_manifest(root, staged_metadata) == expected_source and digest(overlay) == frozen[str(HERE/'export.rs')],
            'compiled sources changed during build')
    executable = out/'target/release/examples/f5_job_schema_control'
    shutil.copyfile(executable, out/'diagnostic')
    (out/'diagnostic').chmod(0o755)
    write_immutable(out/'native-registration.json', dict(binary_sha256=digest(out/'diagnostic'),
        preexecution_sha256=digest(out/'preexecution.json'), inputs_sha256=protocol['input_sha256'],
        argv=[str(out/'diagnostic')], native_seconds=protocol['native_seconds']))
    with (HERE/'inputs.json').open('rb') as stdin:
        supervised([str(out/'diagnostic')], root, env, out, 'native', protocol['native_seconds'], stdin=stdin)
    exported = json.loads((out/'native.stdout').read_text())
    result = audit(json.loads((HERE/'inputs.json').read_text()), exported)
    require(source_manifest(root, staged_metadata) == expected_source and source_manifest(original, metadata) == source
            and digest(overlay) == frozen[str(HERE/'export.rs')]
            and all(digest(path) == expected for path, expected in frozen.items()),
            'control code or Rust input changed during execution')
    require(digest(out/'diagnostic') == json.loads((out/'native-registration.json').read_text())['binary_sha256'],
            'diagnostic executable changed during execution')
    result.update(preexecution_sha256=digest(out/'preexecution.json'),
        native_registration_sha256=digest(out/'native-registration.json'),
        raw_export_sha256=digest(out/'native.stdout'), analysis_sources=frozen,
        historical_source_manifest_sha256=sha256(source), compiled_source_manifest_sha256=sha256(compiled_source))
    write_immutable(out/'RESULT.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--historical-source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, default=HERE/'protocol.json')
    args = parser.parse_args()
    require(not args.out.exists(), 'control output exists; do not retry or overwrite')
    try:
        result = run(args.historical_source.resolve(), args.out.resolve(), args.protocol.resolve())
    except Exception as error:
        if args.out.is_dir() and not (args.out/'RESULT.json').exists():
            write_immutable(args.out/'RESULT.json', dict(schema_version=1, status='CONTROL_FAILURE',
                error_type=type(error).__name__, error=str(error), complete_dlp=False,
                candidate_id=None, measured_costs=None, online_speedup=None, promotion_eligible=False))
        raise
    print(json.dumps(result, sort_keys=True))
