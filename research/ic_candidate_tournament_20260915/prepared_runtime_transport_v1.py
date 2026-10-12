#!/usr/bin/env python3
"""Replay a prepared control audit from its preexecution-frozen sources.

This launches an isolated Python auditor, never a native solver. Admission is
delegated to the registered family's mathematical checker; transport alone
does not establish a solved target, fresh-target qualification or a speedup.
"""
import argparse
import importlib
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from identity import sha256, write_immutable  # noqa: E402
from oracle import InvalidEvidence, require  # noqa: E402
from sat_runtime_bundle import DIRECTORY, check_loaded_modules, source_manifest  # noqa: E402
from sat_runtime_execution_v3 import (  # noqa: E402
    audit_execution, digest, entry_role, frozen_environment, interpreter_record, read,
)

AUDITORS = {
    'prepared_sat_runtime_v1': {
        'ADMITTED_COMPLETE_PREPARED_SAT_CONTROL', 'ADMITTED_INCOMPLETE_PREPARED_SAT_CONTROL'},
    'prepared_f5_runtime_v2': {
        'ADMITTED_COMPLETE_PREPARED_F5_CONTROL', 'ADMITTED_INCOMPLETE_TARGET_PREPARED_F5_CONTROL',
        'AUDITED_NATIVE_TIMEOUT', 'AUDITED_NATIVE_ERROR'},
}
HELPER_ROLE = (DIRECTORY/'prepared_runtime_transport_v1.py').as_posix()
FLAGS = dict(isolated=True, no_site=True, bytecode_writes=False)
TIMEOUT_SECONDS = 180


def auditor_name(spec):
    name = spec['entrypoint']['module']
    require(name in AUDITORS and spec['entrypoint']['callable'] == 'run',
            'prepared transport requires a registered prepared SAT or F5 entrypoint')
    require(HELPER_ROLE in {item['role'] for item in spec['runtime_manifest']['components']},
            'prepared transport helper was not frozen before execution')
    return name


def validate_admission(result, spec):
    name = auditor_name(spec)
    require(type(result) is dict and result.get('status') in AUDITORS[name]
            and result.get('source_bound_execution_admitted') is True
            and result.get('promotion_eligible') is False
            and result.get('fresh_paired_qualification') is False
            and result.get('headline_online_admissible') is False
            and result.get('online_speedup') is None,
            'prepared admission status or development claim boundary changed')
    seal = spec['arguments']['seal']
    require(all(result.get(key) == seal[key] for key in
                ('candidate_id', 'workload_id', 'run_id')),
            'prepared admission differs from registered candidate/workload/run')
    complete = result['status'].startswith('ADMITTED_COMPLETE_')
    require(result.get('scalar_verified') is complete
            and ((type(result.get('online_wall_ns')) is int and result['online_wall_ns'] > 0)
                 if complete else result.get('online_wall_ns') is None),
            'prepared admission completeness differs from scalar or online cost')


def gate(spec, loaded):
    return dict(schema_version=1, execution_sha256=sha256(spec), binding=spec['binding'],
                interpreter_sha256=sha256(spec['interpreter']), flags=FLAGS,
                loaded_modules=loaded)


def validate_gates(before, after, spec):
    roles = {item['role'] for item in spec['runtime_manifest']['components']}
    for item in (before, after):
        require(item == gate(spec, item['loaded_modules'])
                and item['loaded_modules']
                and all(role in roles for role in item['loaded_modules'].values())
                and HELPER_ROLE in item['loaded_modules'].values()
                and entry_role(auditor_name(spec)) in item['loaded_modules'].values(),
                'prepared frozen audit source/interpreter gate differs')
    require(before['loaded_modules'].items() <= after['loaded_modules'].items(),
            'prepared frozen audit dropped previously imported modules')


def frozen_audit(execution, expected_sha256, out):
    execution, out = Path(execution).resolve(), Path(out).resolve()
    spec = read(execution/'execution.json')
    root = HERE.parents[1]
    name = auditor_name(spec)
    require(root == execution/'extracted' and sha256(spec) == expected_sha256
            and sys.flags.isolated == 1 and sys.flags.no_site == 1
            and sys.dont_write_bytecode
            and not out.is_relative_to(execution) and not execution.is_relative_to(out)
            and interpreter_record() == spec['interpreter']
            and source_manifest(root) == spec['runtime_manifest'],
            'prepared audit is outside its frozen source/interpreter or output tree')
    module = importlib.import_module(name)
    before = gate(spec, check_loaded_modules(root, spec['runtime_manifest']))
    write_immutable(out/'before.json', before)
    result = module.audit(execution, spec)
    validate_admission(result, spec)
    after = gate(spec, check_loaded_modules(root, spec['runtime_manifest']))
    require(source_manifest(root) == spec['runtime_manifest']
            and interpreter_record() == spec['interpreter'],
            'prepared audit source or interpreter changed at termination')
    validate_gates(before, after, spec)
    write_immutable(out/'after.json', after)
    write_immutable(out/'admission.json', result)


def transport(execution, expected_sha256, out):
    execution, out = Path(execution).resolve(), Path(out).resolve()
    spec = read(execution/'execution.json')
    auditor_name(spec)
    require(sha256(spec) == expected_sha256 and interpreter_record() == spec['interpreter']
            and not out.exists() and not out.is_relative_to(execution)
            and not execution.is_relative_to(out),
            'prepared transport seal, interpreter or one-use output differs')
    audit_execution(execution, spec)
    out.mkdir(parents=True)
    script = execution/'extracted'/HELPER_ROLE
    command = [sys.executable, '-I', '-S', '-B', str(script), '_frozen-audit',
               '--execution', str(execution), '--expected-execution-sha256', expected_sha256,
               '--out', str(out)]
    timed_out, exit_code, error = False, None, None
    with (out/'stdout.txt').open('x') as stdout, (out/'stderr.txt').open('x') as stderr:
        try:
            result = subprocess.run(command, cwd=execution/'extracted', env=frozen_environment(),
                                    stdout=stdout, stderr=stderr, timeout=TIMEOUT_SECONDS, check=False)
            exit_code = result.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
        except OSError as exception:
            error = dict(type=type(exception).__name__, message=str(exception))
    receipt = dict(schema_version=1, status='REJECTED_PREPARED_TRANSPORT',
                   exit_code=exit_code, timed_out=timed_out, error=error,
                   execution_sha256=expected_sha256, binding=spec['binding'],
                   interpreter_sha256=sha256(spec['interpreter']), python_flags=['-I', '-S', '-B'],
                   watchdog_seconds=TIMEOUT_SECONDS, stdout_sha256=digest(out/'stdout.txt'),
                   stderr_sha256=digest(out/'stderr.txt'), native_solvers_executed=0,
                   promotion_eligible=False, online_speedup=None,
                   scope='independent Python/math audit; no new solver attempt or performance claim')
    if exit_code == 0:
        try:
            validate_gates(read(out/'before.json'), read(out/'after.json'), spec)
            validate_admission(read(out/'admission.json'), spec)
        except (InvalidEvidence, KeyError, TypeError, ValueError, OSError) as exception:
            receipt['error'] = dict(type=type(exception).__name__, message=str(exception))
        else:
            receipt.update(status='PASS_FROZEN_PREPARED_TRANSPORT',
                           admission_sha256=digest(out/'admission.json'),
                           before_sha256=digest(out/'before.json'), after_sha256=digest(out/'after.json'))
    write_immutable(out/'transport.json', receipt)
    require(receipt['status'] == 'PASS_FROZEN_PREPARED_TRANSPORT',
            'prepared frozen transport rejected; retain failure, no native retry')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('audit', '_frozen-audit'):
        child = sub.add_parser(command)
        child.add_argument('--execution', type=Path, required=True)
        child.add_argument('--expected-execution-sha256', required=True)
        child.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'audit':
        import json
        print(json.dumps(transport(args.execution, args.expected_execution_sha256, args.out), sort_keys=True))
    else:
        frozen_audit(args.execution, args.expected_execution_sha256, args.out)


if __name__ == '__main__':
    main()
