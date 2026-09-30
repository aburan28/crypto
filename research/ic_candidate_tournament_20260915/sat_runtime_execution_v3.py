#!/usr/bin/env python3
"""Execute a future SAT controller from retained Python, with source gates.

This proves an execution's Python binding, not an IC solve or a speedup.
Candidate/method/workload registration and mathematical auditing are separate.
Historical SAT controllers may be imported for controls, never remeasured here.
"""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import sysconfig
import time
import traceback

# -I deliberately removes the script directory. Add only this extracted,
# source-bound module directory, never a caller's checkout or PYTHONPATH.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from identity import sha256, write_immutable  # noqa: E402
from oracle import require  # noqa: E402
from sat_runtime_bundle import (  # noqa: E402
    DIRECTORY, check_loaded_modules, extract, freeze, source_manifest,
    verified_files,
)

HERE = Path(__file__).resolve().parent
SOURCE_SUFFIXES = {'.py', '.pyc', '.so', '.dylib', '.pyd', '.zip'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def interpreter_record():
    """Bind the executable and this interpreter's complete stdlib code surface.

    Site packages are forbidden by the module gate. Existing stdlib bytecode
    is included because CPython may read it even when writes are disabled.
    The record retains hashes, not portable copies of a Python installation.
    """
    library = Path(sysconfig.get_path('stdlib')).resolve()
    sites = {Path(sysconfig.get_path(key)).resolve()
             for key in ('purelib', 'platlib')}
    files = {}
    for directory, children, names in os.walk(library):
        root = Path(directory)
        children[:] = sorted(name for name in children
                             if (root/name).resolve() not in sites
                             and name not in ('site-packages', 'dist-packages'))
        for name in sorted(names):
            path = root/name
            if path.suffix not in SOURCE_SUFFIXES:
                continue
            require(path.is_file() and path.resolve().is_relative_to(library),
                    'stdlib source escapes the recorded interpreter')
            data = path.read_bytes()
            files[path.relative_to(library).as_posix()] = dict(
                bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    require(files, 'Python standard-library inventory is empty')
    return dict(schema_version=1, implementation=sys.implementation.name,
                version=sys.version, system=platform.system(),
                machine=platform.machine(), executable_sha256=digest(sys.executable),
                stdlib_inventory=files,
                scope='local interpreter and stdlib content binding; not remote attestation')


def entry_role(module):
    require(type(module) is str and module.isidentifier()
            and not module.startswith('_'), 'invalid bound SAT module name')
    return (DIRECTORY/(module+'.py')).as_posix()


def binding(spec):
    """Algorithm identity fields; invocation arguments remain run data.

    The future registrar must put algorithm settings in its canonical method
    record and freeze the complete invocation in the workload/run registration.
    This utility cannot infer those mathematical or configuration contracts.
    """
    result = dict(python_manifest_sha256=spec['runtime_seal']['manifest_sha256'],
                python_archive_sha256=spec['runtime_seal']['archive_sha256'],
                interpreter_sha256=sha256(spec['interpreter']),
                module=spec['entrypoint']['module'],
                callable=spec['entrypoint']['callable'])
    if spec.get('asset_seal') is not None:
        result['asset_manifest_sha256'] = spec['asset_seal']['manifest_sha256']
        result['asset_archive_sha256'] = spec['asset_seal']['archive_sha256']
    return result


def validate_spec(spec):
    require(spec['schema_version'] == 3,
            'unsupported SAT runtime execution schema')
    role = entry_role(spec['entrypoint']['module'])
    action = spec['entrypoint']['callable']
    require(type(action) is str and action.isidentifier()
            and not action.startswith('_'), 'invalid bound SAT callable')
    require(role in {item['role'] for item in spec['runtime_manifest']['components']},
            'SAT entrypoint is outside the frozen surface')
    require(spec['binding'] == binding(spec),
            'SAT runtime execution binding changed')
    require(('asset_manifest' in spec) == ('asset_seal' in spec),
            'SAT asset manifest and seal must be registered together')


def register(repository, output, *, module, action, arguments, timeout_seconds,
             asset_snapshot=None, arguments_factory=None):
    """Create a new preexecution runtime registration, without executing it."""
    output = Path(output)
    require(not output.exists(), 'SAT runtime registration already exists')
    require(type(timeout_seconds) is int and timeout_seconds > 0,
            'SAT runtime watchdog must be registered as a positive integer')
    roles = {item['role'] for item in source_manifest(repository)['components']}
    require(entry_role(module) in roles,
            'SAT entrypoint source is missing')
    require((DIRECTORY/'sat_runtime_execution_v3.py').as_posix()
            in roles,
            'SAT isolated launcher source is missing')
    # Reject malformed entrypoints before creating any registration directory.
    require(type(action) is str and action.isidentifier()
            and not action.startswith('_'), 'invalid bound SAT callable')
    output.mkdir(parents=True)
    manifest, seal = freeze(repository, output/'runtime')
    spec = dict(schema_version=3, runtime_manifest=manifest, runtime_seal=seal,
                interpreter=interpreter_record(),
                entrypoint=dict(module=module, callable=action),
                arguments=arguments, runtime_watchdog_seconds=timeout_seconds,
                claim_scope='Python execution binding only')
    if asset_snapshot is not None:
        from static_sat_assets_v3 import verified_assets
        spec.update(asset_manifest=read(Path(asset_snapshot)/'manifest.json'),
                    asset_seal=read(Path(asset_snapshot)/'seal.json'))
        verified_assets(asset_snapshot, spec['asset_manifest'], spec['asset_seal'])
        shutil.copytree(asset_snapshot, output/'assets')
    spec['binding'] = binding(spec)
    if arguments_factory is not None:
        require(arguments is None, 'SAT invocation and factory are mutually exclusive')
        spec['arguments'] = arguments_factory(spec)
        require(source_manifest(repository) == manifest,
                'SAT source changed while mathematical registration was constructed')
    validate_spec(spec)
    write_immutable(output/'execution.json', spec)
    write_immutable(output/'registration-seal.json', dict(
        schema_version=3, execution_sha256=sha256(spec), binding=spec['binding'],
        registration_stage='before-execution'))
    return spec


def frozen_environment():
    # Isolation flags also disable user/site initialization. Keep a minimal
    # environment for the separately source-bound native toolchain/solver.
    result = {key: value for key, value in os.environ.items()
              if key in ('PATH', 'HOME', 'TMPDIR', 'LANG', 'LC_ALL')}
    result.update(PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1')
    return result


def make_read_only(root):
    for path in sorted(Path(root).rglob('*'), reverse=True):
        if path.is_dir():
            path.chmod(0o555)
    Path(root).chmod(0o555)


def execute(registration, output, *, expected_spec, timeout_seconds):
    """Launch once from a verified extraction and retain even failed output."""
    registration, output = map(lambda p: Path(p).resolve(), (registration, output))
    require(not output.exists(), 'SAT execution output exists; no retries')
    require(type(timeout_seconds) is int and timeout_seconds > 0,
            'SAT runtime watchdog must be a positive integer')
    require(os.name == 'posix', 'SAT process-group watchdog requires a POSIX host')
    spec = read(registration/'execution.json')
    require(spec == expected_spec,
            'SAT execution differs from independently sealed registration')
    require(timeout_seconds == spec['runtime_watchdog_seconds'],
            'SAT runtime watchdog differs from sealed registration')
    validate_spec(spec)
    registration_seal = read(registration/'registration-seal.json')
    require(registration_seal == dict(schema_version=3, execution_sha256=sha256(spec),
                binding=spec['binding'], registration_stage='before-execution'),
            'SAT complete invocation seal changed')
    verified_files(registration/'runtime', spec['runtime_manifest'], spec['runtime_seal'])
    if spec.get('asset_seal') is not None:
        from static_sat_assets_v3 import verified_assets
        verified_assets(registration/'assets', spec['asset_manifest'], spec['asset_seal'])
    require(interpreter_record() == spec['interpreter'],
            'SAT registered interpreter changed before execution')
    output.mkdir(parents=True)
    write_immutable(output/'execution.json', spec)
    write_immutable(output/'registration-seal.json', registration_seal)
    shutil.copytree(registration/'runtime', output/'runtime')
    root = extract(output/'runtime', output/'extracted',
                   spec['runtime_manifest'], spec['runtime_seal'])
    make_read_only(root)
    if spec.get('asset_seal') is not None:
        from static_sat_assets_v3 import extract_assets
        shutil.copytree(registration/'assets', output/'assets')
        asset_root = extract_assets(output/'assets', output/'asset-files',
                                   spec['asset_manifest'], spec['asset_seal'])
        make_read_only(asset_root)
    (output/'entry-output').mkdir()
    command = [sys.executable, '-I', '-S', '-B',
               str(root/DIRECTORY/'sat_runtime_execution_v3.py'),
               '--worker', str(output)]
    started = time.monotonic_ns()
    timed_out = False
    with (output/'stdout.txt').open('x') as stdout, (output/'stderr.txt').open('x') as stderr:
        process = subprocess.Popen(command, cwd=root, env=frozen_environment(),
                                   stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            returncode = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            # The future controller may have solver children. Kill the process
            # group, wait for the leader, and retain the incomplete source gate.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass  # The whole process group exited at the watchdog boundary.
            returncode = process.wait()
        finally:
            # A failed controller can leave an inherited-group meter alive.
            # Reap the entire owned group on every terminal path, including
            # normal leader exit. Native v3 tools never start another group.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    wall = time.monotonic_ns()-started
    receipt = dict(schema_version=3, execution_sha256=sha256(spec),
                   binding=spec['binding'], exit_code=returncode, timed_out=timed_out,
                   watchdog_seconds=timeout_seconds, process_wall_ns=wall,
                   python_flags=['-I', '-S', '-B'],
                   stdout_sha256=digest(output/'stdout.txt'),
                   stderr_sha256=digest(output/'stderr.txt'),
                   status='EXECUTED_SOURCE_CONTROL',
                   scope='whole runtime process; outside single-target online interval',
                   promotion_eligible=False, online_speedup=None)
    write_immutable(output/'process.json', receipt)
    return receipt


def worker(output):
    output = Path(output).resolve()
    root = HERE.parents[1]
    spec = read(output/'execution.json')
    validate_spec(spec)
    require(sys.flags.isolated == 1 and sys.flags.no_site == 1
            and sys.dont_write_bytecode, 'SAT Python process is not isolated')
    require(source_manifest(root) == spec['runtime_manifest']
            and interpreter_record() == spec['interpreter'],
            'SAT extracted sources or interpreter changed before import')
    if spec.get('asset_seal') is not None:
        from static_sat_assets_v3 import check_extracted_assets
        check_extracted_assets(output/'asset-files', spec['asset_manifest'])
    module = importlib.import_module(spec['entrypoint']['module'])
    action = getattr(module, spec['entrypoint']['callable'])
    require(callable(action), 'SAT registered entrypoint is not callable')
    before = check_loaded_modules(root, spec['runtime_manifest'])
    pre = dict(schema_version=3, execution_sha256=sha256(spec),
               binding=spec['binding'], loaded_modules=before,
               interpreter_sha256=sha256(spec['interpreter']),
               flags=dict(isolated=True, no_site=True, bytecode_writes=False))
    write_immutable(output/'before.json', pre)
    succeeded, error = False, None
    result = None
    try:
        result = action(spec['arguments'], output/'entry-output')
        succeeded = True
    except Exception as exception:
        error = dict(type=type(exception).__name__, message=str(exception))
        traceback.print_exc()
    # A thrown entrypoint still needs source evidence. A failing terminal
    # source/import gate leaves no valid after.json and cannot be admitted.
    after = check_loaded_modules(root, spec['runtime_manifest'])
    if spec.get('asset_seal') is not None:
        check_extracted_assets(output/'asset-files', spec['asset_manifest'])
    require(source_manifest(root) == spec['runtime_manifest']
            and interpreter_record() == spec['interpreter'],
            'SAT extracted sources or interpreter changed at termination')
    write_immutable(output/'after.json', dict(pre, loaded_modules=after,
                    entrypoint_succeeded=succeeded, error=error, result=result))
    return 0 if succeeded else 1


def audit_execution(output, expected_spec):
    """Check retained source gates against an independently sealed spec."""
    output = Path(output)
    validate_spec(expected_spec)
    require(read(output/'execution.json') == expected_spec,
            'retained SAT execution differs from sealed registration')
    require(read(output/'registration-seal.json') == dict(
                schema_version=3, execution_sha256=sha256(expected_spec),
                binding=expected_spec['binding'], registration_stage='before-execution'),
            'SAT retained invocation seal changed')
    files = verified_files(output/'runtime', expected_spec['runtime_manifest'],
                           expected_spec['runtime_seal'])
    require(source_manifest(output/'extracted') == expected_spec['runtime_manifest'],
            'retained SAT extraction differs from registered sources')
    if expected_spec.get('asset_seal') is not None:
        from static_sat_assets_v3 import check_extracted_assets, verified_assets
        verified_assets(output/'assets', expected_spec['asset_manifest'], expected_spec['asset_seal'])
        check_extracted_assets(output/'asset-files', expected_spec['asset_manifest'])
    process, before, after = map(read, (output/'process.json',
                                      output/'before.json', output/'after.json'))
    require(process['execution_sha256'] == sha256(expected_spec)
            and process['binding'] == expected_spec['binding']
            and process['timed_out'] is False
            and process['python_flags'] == ['-I', '-S', '-B']
            and type(process['process_wall_ns']) is int
            and process['process_wall_ns'] >= 0
            and type(process['watchdog_seconds']) is int
            and process['watchdog_seconds'] > 0
            and process['watchdog_seconds'] == expected_spec['runtime_watchdog_seconds']
            and process['stdout_sha256'] == digest(output/'stdout.txt')
            and process['stderr_sha256'] == digest(output/'stderr.txt')
            and process['promotion_eligible'] is False
            and process['online_speedup'] is None,
            'SAT process/source evidence is incomplete or changed')
    for gate in (before, after):
        require(gate['execution_sha256'] == sha256(expected_spec)
                and gate['binding'] == expected_spec['binding']
                and gate['interpreter_sha256'] == sha256(expected_spec['interpreter'])
                and gate['flags'] == dict(isolated=True, no_site=True, bytecode_writes=False)
                and gate['loaded_modules']
                and all(role in files for role in gate['loaded_modules'].values())
                and entry_role(expected_spec['entrypoint']['module'])
                    in gate['loaded_modules'].values(),
                'SAT imported-module gate differs from registered surface')
    require(before['loaded_modules'].items() <= after['loaded_modules'].items(),
            'SAT terminal gate drops previously loaded modules')
    require(type(after['entrypoint_succeeded']) is bool
            and process['exit_code'] == (0 if after['entrypoint_succeeded'] else 1)
            and (after['error'] is None) is after['entrypoint_succeeded'],
            'SAT entrypoint outcome differs from retained process')
    return dict(schema_version=3, status='AUDITED_PYTHON_EXECUTION_BINDING',
                binding=expected_spec['binding'], complete_source_gates=True,
                entrypoint_succeeded=after['entrypoint_succeeded'],
                loaded_module_count=len(after['loaded_modules']),
                claim_scope='Python binding only; mathematical admission required separately',
                promotion_eligible=False, online_speedup=None)


def import_probe(arguments, output):
    """Import-only control for real historical modules; never call their jobs."""
    require(type(arguments) is list and all(type(name) is str for name in arguments),
            'import probe requires module names')
    for name in arguments:
        importlib.import_module(name)
    return dict(imported=arguments, measured_solver_executed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', type=Path, required=True)
    arguments = parser.parse_args()
    try:
        sys.exit(worker(arguments.worker))
    except Exception:
        traceback.print_exc()
        sys.exit(2)
