#!/usr/bin/env python3
"""Actual native-to-Python prepared-report controls on the disclosed n17 fixture.

These repeatable interface tests are not a registered scientific IC execution,
natural-yield estimate, full Python runtime attestation or performance result.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import generic_build
from identity import canonical, write_immutable
from oracle import InvalidEvidence, require
from prepared_target_v1 import CERTIFICATE_SEALS, audit_native_target, native_job
from sat_runtime_execution_v3 import frozen_environment
from static_sat_native_v3 import THREAD_ENVIRONMENT

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
POINT = [52411, 72106]


def run_controls(build, out):
    build, out = Path(build).resolve(), Path(out).resolve()
    require(not out.exists(), 'prepared report control output already exists')
    record = json.loads((build/'build-record.json').read_text())
    source = json.loads((build/'source-manifest.json').read_text())
    worker = build/'worker'
    generic_build.verify_build_record(record, source)

    def source_gate():
        require(generic_build.digest(worker) == record['worker_sha256']
                and all(generic_build.digest(ROOT/name) == value
                        for name, value in source['root_files'].items()),
                'native control executable or root source changed')

    source_gate()
    document = json.loads((HERE/'goal_20260924/prepared-ic-state-v1/f5-preparation.json').read_text())
    out.mkdir(parents=True)
    write_immutable(out/'build-record.json', record)
    write_immutable(out/'source-manifest.json', source)
    write_immutable(out/'preparation.json', document)
    env = dict(frozen_environment(), **THREAD_ENVIRONMENT)
    controls = []
    for cap, expected in ((1, 'INCOMPLETE_TARGET'), (8, 'COMPLETE')):
        directory = out/f'cap-{cap}'
        directory.mkdir()
        job = native_job(document, CERTIFICATE_SEALS['f5'], point=POINT,
                         algorithm_seed=2026093032, max_attempts=cap)
        write_immutable(directory/'job.json', job)
        source_gate()
        with (directory/'stdout.json').open('xb') as stdout, (directory/'stderr.txt').open('xb') as stderr:
            process = subprocess.Popen([str(worker)], stdin=subprocess.PIPE, stdout=stdout,
                                       stderr=stderr, env=env, start_new_session=True)
            timed_out = False
            try:
                process.communicate(canonical(job), timeout=90)
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                if process.poll() is None:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                process.communicate()
        receipt = dict(exit_code=process.returncode, timed_out=timed_out, watchdog_seconds=90,
                       stdout_sha256=generic_build.digest(directory/'stdout.json'),
                       stderr_sha256=generic_build.digest(directory/'stderr.txt'),
                       environment_threads=THREAD_ENVIRONMENT,
                       process_group_policy='one standalone control session; kill whole group on timeout')
        write_immutable(directory/'process.json', receipt)
        source_gate()
        require(not timed_out and process.returncode == (0 if cap == 8 else 2),
                'actual native prepared control failed; retain output')
        report = json.loads((directory/'stdout.json').read_text())
        binding = generic_build.verify_binding(report, record, source, executable=worker)
        audit = audit_native_target(report, job, document, CERTIFICATE_SEALS['f5'])
        require(report.get('query_schema_version') == 1 and audit['status'] == expected
                and audit['ordinary_queries_executed'] == 0
                and audit['scalar_verified'] is (cap == 8),
                'actual native report failed its complete/incomplete interface contract')
        broken = copy.deepcopy(report)
        del broken['query_schema_version']
        try:
            audit_native_target(broken, job, document, CERTIFICATE_SEALS['f5'])
        except InvalidEvidence as error:
            require(str(error) == 'missing query schema', 'missing-header control rejected for another reason')
        else:
            raise InvalidEvidence('missing query schema was accepted')
        write_immutable(directory/'mathematical-audit.json', audit)
        controls.append(dict(cap=cap, mathematical_status=audit['status'],
                             scalar_verified=audit['scalar_verified'],
                             retained_attempts=audit['query_audit']['descent_queries'],
                             native_build_binding=binding, missing_header_rejected=True))
    result = dict(schema_version=1, status='PASS_ACTUAL_PREPARED_REPORT_INTERFACE_CONTROLS',
                  controls=controls, target_point=POINT, fresh_targets_generated=0,
                  script_sha256=generic_build.digest(__file__),
                  source_bound_scientific_runtime_admitted=False,
                  verified_scientific_online_wall_ns=None, promotion_eligible=False,
                  online_speedup=None,
                  scope='repeatable disclosed native-build/math interface controls; not a campaign or full Python runtime attestation')
    write_immutable(out/'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_controls(args.build, args.out), sort_keys=True))
