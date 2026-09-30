"""Fresh, isolated meters for the two admitted non-forking SAT native tools.

Native tools inherit the enclosing controller process group. This is not a
general supervisor for arbitrary plugins that fork or create their own groups.
All wrapper/source checks are charged to the caller's PDP interval.
"""
import argparse
import json
import math
import os
from pathlib import Path
import platform
import resource
import signal
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from identity import sha256, write_immutable  # noqa: E402
from oracle import require  # noqa: E402
from sat_runtime_bundle import DIRECTORY, check_loaded_modules, source_manifest  # noqa: E402
from sat_runtime_execution_v3 import digest, interpreter_record, read, validate_spec  # noqa: E402
from static_sat_assets_v3 import check_extracted_assets  # noqa: E402
from tournament import write  # noqa: E402

HERE = Path(__file__).resolve().parent
THREAD_ENVIRONMENT = {name: '1' for name in (
    'RAYON_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}


def context(execution):
    execution = Path(execution).resolve()
    spec = read(execution/'execution.json')
    validate_spec(spec)
    require(spec.get('asset_seal') is not None, 'native SAT assets are not sealed')
    return execution, spec


def meter(execution, asset_role, arguments, directory, name, seconds):
    """Run once. A wrapper failure remains evidence and fails the controller."""
    execution, spec = context(execution)
    directory = Path(directory).resolve()
    require(directory.is_relative_to(execution/'entry-output')
            and directory.is_dir(), 'native SAT output escapes its execution')
    require(name.isidentifier() and type(seconds) is int and seconds > 0,
            'invalid native SAT meter name or watchdog')
    call = dict(schema_version=3, execution_sha256=sha256(spec),
                asset_role=asset_role, arguments=[str(x) for x in arguments],
                watchdog_seconds=seconds, name=name,
                directory=directory.relative_to(execution).as_posix())
    intent = directory/(name+'.intent.json')
    require(not intent.exists(), 'native SAT invocation already exists; no retries')
    write_immutable(intent, call)
    command = [sys.executable, '-I', '-S', '-B', str(HERE/'static_sat_native_v3.py'),
               '--execution', str(execution), '--intent', str(intent),
               '--expected-sha256', sha256(call)]
    with (directory/(name+'.wrapper.stdout')).open('x') as stdout, (
            directory/(name+'.wrapper.stderr')).open('x') as stderr:
        # The outer controller watchdog owns this group, including the meter.
        process = subprocess.run(command, stdout=stdout, stderr=stderr,
                                 start_new_session=False, check=False)
    require(process.returncode == 0, 'native SAT meter failed; retain partial evidence')
    return read(directory/(name+'.metrics.json'))


def worker(execution, intent, expected_sha256):
    execution, spec = context(execution)
    call = read(intent)
    require(sha256(call) == expected_sha256
            and call['execution_sha256'] == sha256(spec),
            'native SAT invocation differs from its sealed parent')
    directory = (execution/call['directory']).resolve()
    name = call['name']
    require(directory.is_relative_to(execution/'entry-output')
            and Path(intent).resolve() == directory/(name+'.intent.json')
            and name.isidentifier(), 'unsafe native SAT invocation output')
    require(type(call['watchdog_seconds']) is int and call['watchdog_seconds'] > 0
            and type(call['arguments']) is list
            and all(type(x) is str for x in call['arguments']),
            'invalid native SAT invocation arguments')
    root = HERE.parents[1]
    require(root == execution/'extracted' and sys.flags.isolated
            and sys.flags.no_site and sys.dont_write_bytecode,
            'native SAT meter is outside the isolated extraction')
    require(source_manifest(root) == spec['runtime_manifest']
            and interpreter_record() == spec['interpreter'],
            'native SAT meter source or interpreter changed')
    assets = check_extracted_assets(execution/'asset-files', spec['asset_manifest'])
    roles = {item['role'] for item in spec['asset_manifest']['components']
             if item['executable']}
    require(call['asset_role'] in roles, 'native SAT tool is not an executable asset')
    executable = execution/'asset-files'/call['asset_role']
    command = [str(executable)]+call['arguments']
    loaded = check_loaded_modules(root, spec['runtime_manifest'])
    gate = dict(schema_version=3, execution_sha256=sha256(spec),
                invocation_sha256=sha256(call), binding=spec['binding'],
                executable_sha256=digest(executable),
                interpreter_sha256=sha256(spec['interpreter']),
                loaded_modules=loaded,
                process_group=os.getpgrp(), wrapper_pid=os.getpid(),
                flags=dict(isolated=True, no_site=True, bytecode_writes=False))
    gate['execution_directory'] = str(execution)
    write_immutable(directory/(name+'.before.json'), gate)
    environment = dict(os.environ, **THREAD_ENVIRONMENT)
    timed_out = False
    process = None
    start = time.monotonic_ns()
    with (directory/(name+'.stdout')).open('xb') as stdout, (
            directory/(name+'.stderr')).open('xb') as stderr:
        try:
            process = subprocess.Popen(command, cwd=directory, env=environment,
                                       stdout=stdout, stderr=stderr,
                                       start_new_session=False)
            native_pid = process.pid
            try:
                native_group = os.getpgid(native_pid)
                require(native_group == os.getpgrp(),
                        'native SAT tool escaped the controller watchdog group')
            except ProcessLookupError:
                # A very short preflight may finish before getpgid. Popen's
                # inherited-group configuration is still recorded explicitly.
                native_group = None
            write_immutable(directory/(name+'.spawned.json'), dict(
                native_pid=native_pid, native_group_observed=native_group,
                controller_group=os.getpgrp(), invocation_sha256=sha256(call)))
            try:
                returncode = process.wait(timeout=call['watchdog_seconds'])
            except subprocess.TimeoutExpired:
                timed_out = True
                process.send_signal(signal.SIGKILL)
                returncode = process.wait()
        finally:
            if process is not None and process.poll() is None:
                process.kill()
                process.wait()
    stop = time.monotonic_ns()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    require(digest(executable) == gate['executable_sha256']
            and check_extracted_assets(execution/'asset-files', spec['asset_manifest']) == assets
            and source_manifest(root) == spec['runtime_manifest']
            and interpreter_record() == spec['interpreter'],
            'native SAT sources or assets changed during execution')
    after = check_loaded_modules(root, spec['runtime_manifest'])
    write_immutable(directory/(name+'.after.json'), dict(gate, loaded_modules=after))
    metrics = dict(wall_seconds=(stop-start)/1_000_000_000,
                   user_seconds=usage.ru_utime, system_seconds=usage.ru_stime,
                   total_core_seconds=usage.ru_utime+usage.ru_stime,
                   single_core_seconds=usage.ru_utime+usage.ru_stime,
                   peak_rss_bytes=int(usage.ru_maxrss)*(1 if platform.system() == 'Darwin' else 1024),
                   meter='fresh isolated wrapper RUSAGE_CHILDREN; native process only')
    receipt = dict(schema_version=3, command=command, returncode=returncode,
                   native_wall_ns=stop-start,
                   timed_out=timed_out, watchdog_seconds=call['watchdog_seconds'],
                   metrics=metrics, execution_sha256=sha256(spec),
                   invocation_sha256=sha256(call), executable_sha256=digest(executable),
                   environment_threads=THREAD_ENVIRONMENT,
                   process_group_policy='inherit-controller-group; no native fork or setsid',
                   native_pid=native_pid, native_group_observed=native_group,
                   controller_group=gate['process_group'],
                   stdout_sha256=digest(directory/(name+'.stdout')),
                   stderr_sha256=digest(directory/(name+'.stderr')),
                   before_sha256=digest(directory/(name+'.before.json')),
                   after_sha256=digest(directory/(name+'.after.json')))
    write(directory/(name+'.metrics.json'), receipt, exclusive=True)


def audit_meter(execution, directory, name, *, asset_role, arguments, seconds):
    """Check exact argv, native bytes, isolated child gates and output receipts."""
    execution, spec = context(execution)
    directory = Path(directory).resolve()
    call = read(directory/(name+'.intent.json'))
    expected = dict(schema_version=3, execution_sha256=sha256(spec),
                    asset_role=asset_role, arguments=[str(x) for x in arguments],
                    watchdog_seconds=seconds, name=name,
                    directory=directory.relative_to(execution).as_posix())
    require(call == expected, 'native SAT command or watchdog differs from registration')
    receipt = read(directory/(name+'.metrics.json'))
    spawned = read(directory/(name+'.spawned.json'))
    require(spawned == dict(native_pid=receipt['native_pid'],
                            native_group_observed=receipt['native_group_observed'],
                            controller_group=receipt['controller_group'],
                            invocation_sha256=sha256(call)),
            'native SAT launch record differs from terminal receipt')
    binary = execution/'asset-files'/asset_role
    require(receipt['command'][1:] == expected['arguments']
            and receipt['execution_sha256'] == sha256(spec)
            and receipt['invocation_sha256'] == sha256(call)
            and receipt['executable_sha256'] == digest(binary)
            and receipt['watchdog_seconds'] == seconds
            and type(receipt['timed_out']) is bool
            and type(receipt['returncode']) is int
            and type(receipt['native_wall_ns']) is int and receipt['native_wall_ns'] >= 0
            and receipt['environment_threads'] == THREAD_ENVIRONMENT
            and receipt['process_group_policy'] == 'inherit-controller-group; no native fork or setsid'
            and receipt['native_group_observed'] in (None, receipt['controller_group']),
            'native SAT process receipt changed')
    files = {item['role'] for item in spec['runtime_manifest']['components']}
    gates = []
    for stage in ('before', 'after'):
        path = directory/(name+'.'+stage+'.json')
        gate = read(path)
        require(receipt[stage+'_sha256'] == digest(path)
                and gate['execution_sha256'] == sha256(spec)
                and gate['invocation_sha256'] == sha256(call)
                and gate['binding'] == spec['binding']
                and gate['executable_sha256'] == digest(binary)
                and gate['interpreter_sha256'] == sha256(spec['interpreter'])
                and gate['process_group'] == receipt['controller_group']
                and receipt['command'][0] == str(Path(gate['execution_directory'])/'asset-files'/asset_role)
                and gate['flags'] == dict(isolated=True, no_site=True, bytecode_writes=False)
                and (DIRECTORY/'static_sat_native_v3.py').as_posix()
                    in gate['loaded_modules'].values()
                and set(gate['loaded_modules'].values()) <= files,
                'native SAT child source gate missing or changed')
        gates.append(gate)
    require(gates[0]['loaded_modules'].items() <= gates[1]['loaded_modules'].items(),
            'native SAT terminal gate drops imports')
    for stream in ('stdout', 'stderr'):
        require(receipt[stream+'_sha256'] == digest(directory/(name+'.'+stream)),
                'native SAT output bytes changed')
    require(all(type(receipt['metrics'][key]) in (int, float)
                and math.isfinite(receipt['metrics'][key])
                and receipt['metrics'][key] >= 0 for key in (
                    'wall_seconds', 'user_seconds', 'system_seconds',
                    'total_core_seconds', 'single_core_seconds', 'peak_rss_bytes')),
            'native SAT resource metrics missing or negative')
    return receipt


def process_control(arguments, output):
    """Non-SAT process-group control for the bounded supervisor tests only."""
    receipt = meter(Path(output).parent, arguments['role'], arguments['argv'],
                    output, 'control', arguments['seconds'])
    return dict(native_returncode=receipt['returncode'], timed_out=receipt['timed_out'],
                scope='process supervisor control; no IC or SAT admission')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execution', type=Path, required=True)
    parser.add_argument('--intent', type=Path, required=True)
    parser.add_argument('--expected-sha256', required=True)
    args = parser.parse_args()
    worker(args.execution, args.intent, args.expected_sha256)
