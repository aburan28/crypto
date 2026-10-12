"""Execute the sealed native one-target F4/F5 job from a frozen source snapshot."""
import argparse
import json
import os
from pathlib import Path
import platform

from f5_runtime_inputs_v1 import native_admission
from f5_runtime_registration_v2 import mathematical_registration
from identity import sha256, write_immutable
from oracle import require
from sat_runtime_execution_v3 import execute, read
from static_sat_assets_v3 import check_extracted_assets
from static_sat_native_v3 import meter


def run(arguments, output):
    output = Path(output)
    execution = output.parent
    spec = read(execution/'execution.json')
    assets = check_extracted_assets(execution/'asset-files', spec['asset_manifest'])
    fixture, inventory, curve, base, build, rust_source, native = native_admission(assets)
    require(arguments == spec['arguments']
            == mathematical_registration(arguments['panel'], spec, assets),
            'F5 supplied invocation differs from full mathematical registration')
    for key, value in arguments.items():
        write_immutable(output/(key+'.json'), value)
    write_immutable(output/'host.json', dict(system=platform.system(), machine=platform.machine(),
        os_release=platform.release(), cpu_count=os.cpu_count(), cpu_model=None,
        scope='uncontrolled local development correctness; no performance claim'))
    preflight = meter(execution, 'bin/worker', ['--build-identity'], output, 'build_identity', 10)
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and read(output/'build_identity.stdout') == native['build_identity'],
            'F5 native worker embeds a different accepted build identity')
    # Reserve ending gates below the whole-controller watchdog. Both caps are
    # fixed by registration; this does not extend or restart a native attempt.
    seconds = arguments['panel']['resources']['total_wall_limit_seconds']-30
    result = meter(execution, 'bin/worker', [], output, 'pipeline', seconds, stdin_argument='job')
    if result['timed_out']:
        write_immutable(output/'summary.json', dict(status='NATIVE_TIMEOUT',
            candidate_id=arguments['seal']['candidate_id'], workload_id=arguments['seal']['workload_id'],
            run_id=arguments['seal']['run_id'], complete_ic_admitted=False,
            final_rank=None, verified_target_count=None, online_wall_ns=None, online_speedup=None))
        return dict(status='NATIVE_TIMEOUT', run_id=arguments['seal']['run_id'])
    report = read(output/'pipeline.stdout')
    require((result['returncode'] == 0 and report['status'] == 'complete')
            or (result['returncode'] == 2 and report['status'] == 'incomplete'),
            'F5 native process result is inconsistent or failed')
    write_immutable(output/'summary.json', dict(status=report['status'].upper(),
        run_id=arguments['seal']['run_id'], independent_mathematical_audit='pending',
        complete_ic_admitted=False, online_speedup=None))
    return dict(status=report['status'].upper(), run_id=arguments['seal']['run_id'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--expected-execution-sha256', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    spec = read(args.registration/'execution.json')
    require(sha256(spec) == args.expected_execution_sha256,
            'F5 dispatch differs from externally recorded invocation hash')
    result = execute(args.registration, args.out, expected_spec=spec,
                     timeout_seconds=spec['runtime_watchdog_seconds'])
    print(json.dumps(result, sort_keys=True))
