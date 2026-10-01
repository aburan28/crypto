"""One-use prepared n17 F5 development runtime and independent transport.

The disclosed point is the only admitted target. Fresh qualification requires
a separate reviewed protocol. All old registrations remain consumed; no API
here extends them, samples targets, subtracts observer cost or asserts a win.
"""
import argparse
import copy
import json
from pathlib import Path
import platform
import subprocess
import sys

# The transport CLI executes this exact file from its read-only extraction
# under -I -S -B. Imports may use only that registered directory and stdlib.
HERE = Path(__file__).resolve().parent
if __package__ in (None, ''):
    sys.path.insert(0, str(HERE))

from f5_runtime_receipts_v2 import publication_process  # noqa: E402
from generic_admission import method_record  # noqa: E402
from generic_bases import verify_base  # noqa: E402
from generic_build import verify_binding  # noqa: E402
from generic_stages import DEFAULTS  # noqa: E402
from identity import candidate_manifest, run_id, sha256, workload_manifest, write_immutable  # noqa: E402
from oracle import require  # noqa: E402
from prepared_f5_inputs_v1 import native_admission  # noqa: E402
from prepared_target_v1 import audit_native_target, native_job  # noqa: E402
from sat_runtime_bundle import check_loaded_modules, source_manifest  # noqa: E402
from sat_runtime_execution_v3 import (audit_execution, digest, execute,  # noqa: E402
    frozen_environment, interpreter_record, read, register as register_runtime)
from static_sat_assets_v3 import check_extracted_assets, verified_assets  # noqa: E402
from static_sat_native_v3 import audit_meter, meter  # noqa: E402

SETTINGS = {'question', 'resources', 'target_input', 'algorithm_seed', 'max_attempts', 'run_number'}
QUESTION = 'prepared-development-source-control'
TARGET = dict(point=[52411,72106], seed=None,
    input_law='one-disclosed-public-point; fixture-construction-excluded',
    point_was_previously_supplied=True, known_scalar_supplied=False)


def mathematical_registration(panel, spec, files, document, certificate_sha256):
    build, _, native = native_admission(files)
    require(type(panel) is dict and set(panel) == SETTINGS and panel['question'] == QUESTION,
            'prepared F5 fresh qualification needs the separately reviewed paired protocol')
    for key in ('algorithm_seed','max_attempts','run_number'):
        require(type(panel[key]) is int and 0 <= panel[key] < 2**64,
                'prepared F5 setting is not an exact bounded integer: '+key)
    require(1 <= panel['max_attempts'] <= 8 and panel['target_input'] == TARGET
            and panel['target_input']['point_was_previously_supplied'] is True
            and panel['target_input']['known_scalar_supplied'] is False,
            'prepared F5 requires the disclosed control point without a scalar')
    resources = panel['resources']
    require(type(resources) is dict and resources == dict(host_class=native['host_class'],
                cpu_workers=1, target_count=1, memory_limit_bytes=None,
                total_wall_limit_seconds=resources.get('total_wall_limit_seconds'))
            and type(resources['cpu_workers']) is int and type(resources['target_count']) is int
            and type(resources['total_wall_limit_seconds']) is int
            and 31 <= resources['total_wall_limit_seconds'] <= 900,
            'prepared F5 requires a bounded one-worker control envelope')
    require(spec['entrypoint'] == dict(module='prepared_f5_runtime_v1',callable='run')
            and spec['runtime_watchdog_seconds'] == resources['total_wall_limit_seconds'],
            'prepared F5 entrypoint or controller watchdog differs')
    job = native_job(document, certificate_sha256, point=panel['target_input']['point'],
        algorithm_seed=panel['algorithm_seed'], max_attempts=panel['max_attempts'])
    fixture = copy.deepcopy(document['certificate']['inputs']['fixture'])
    inventory = dict(fixture=fixture, factor_base=job['prepared']['factor_base'], columns=29,
        column_logs=job['prepared']['columns'], effective_factor_base=job['factor_base'],
        effective_config=dict(copy.deepcopy(DEFAULTS), **job['config']))
    stages = {'base':verify_base(inventory, fixture, job)}
    collector = dict(strategy='Groebner', field_kernel='runtime-capability-gated',
        pair_table=False, query_rule='none', collection_window=None)
    method = method_record(job, inventory, stages, build, collector_plan=collector)
    source = dict(execution_binding=spec['binding'], native=native)
    code_binding = {k:v for k,v in spec['binding'].items() if not k.startswith('asset_')}
    code_binding['native_components'] = dict(
        rust_source_manifest_sha256=native['source_manifest_sha256'],
        worker_binary_sha256=native['worker_sha256'],
        compile_flags=build['build']['flags'])
    code_sha = sha256(code_binding)
    for stage in ('point_decomposition','relation_collection','relation_linear_algebra','target_descent'):
        method[stage]['source_sha256'] = code_sha
    method['relation_collection'].update(query_distribution='reusable-target-independent-preparation',
        query_rule='import-certified-ordinary-rows; no new ordinary queries',
        stop_rule='independent-full-rank-and-all-column-log-replay-before-online')
    method['implementation'] = dict(source_manifest_sha256=code_sha,
        components=[dict(role='full-frozen-Python-interpreter-Rust-source-and-worker-binding',sha256=code_sha)],
        flags=dict(algorithm_execution_binding=code_binding,
            mathematical_preparation_sha256=document['record_sha256'],
            preparation_policy='replay-all-ordinary-rows-and-logs-before-online-v1',
            reusable_template='fixed-zero-symbolic-specialization-before-online; no query',
            field_kernel='runtime-capability-gated; observed-and-audited-per-run',
            runtime_policy='default-environment-one-rayon-v1',
            online_interval='first-target-query-through-general-group-scalar-replay',
            observer_policy='exclusive-owner-thread-v1; all target attempts charged',
            stdin_policy='canonical-registered-math-only-preparation-and-public-point',
            native_watchdog_group='inherit-controller-group-no-native-fork-or-setsid'))
    candidate = candidate_manifest(fixture, inventory, method)
    target_fixture = dict(fixture, targets=job['public_targets'], target_seeds=[None])
    workload = workload_manifest(target_fixture, input_law=TARGET['input_law'],
        algorithm_seed=panel['algorithm_seed'], resource_envelope=resources, cache_policy='warm')
    workload['record'].update(question=QUESTION, point_was_previously_supplied=True)
    workload['record_sha256'] = sha256(workload['record'])
    workload['workload_id'] = workload['record_sha256'][:12]
    seal = dict(registration_stage='before-execution', panel_sha256=sha256(panel),
        candidate_id=candidate['candidate_id'], workload_id=workload['workload_id'],
        run_id=run_id(candidate['candidate_id'],workload['workload_id'],panel['run_number']),
        candidate_sha256=sha256(candidate), workload_sha256=sha256(workload),
        source_sha256=sha256(source), method_sha256=sha256(method),
        preparation_certificate_sha256=certificate_sha256)
    return dict(panel=copy.deepcopy(panel), job=job, fixture=target_fixture,
        source=source, method=method, candidate=candidate, workload=workload, seal=seal,
        preparation=copy.deepcopy(document), preparation_certificate_sha256=certificate_sha256)


def register(repository, assets, panel, document, certificate_sha256, output):
    assets = Path(assets)
    files = verified_assets(assets,read(assets/'manifest.json'),read(assets/'seal.json'))
    native_admission(files,check_host=True)
    return register_runtime(repository,output,module='prepared_f5_runtime_v1',action='run',arguments=None,
        timeout_seconds=panel['resources']['total_wall_limit_seconds'],asset_snapshot=assets,
        arguments_factory=lambda spec:mathematical_registration(panel,spec,files,document,certificate_sha256))


def summary_for(arguments, status):
    return dict(status=status, **{k:arguments['seal'][k] for k in ('candidate_id','workload_id','run_id')},
        independent_mathematical_audit='pending', source_bound_execution_admitted=False,
        fresh_paired_qualification=False, promotion_eligible=False, online_speedup=None)


def run(arguments, out):
    out = Path(out).resolve()
    execution = out.parent
    spec = read(execution/'execution.json')
    files = check_extracted_assets(execution/'asset-files',spec['asset_manifest'])
    _,_,native = native_admission(files,check_host=True)
    require(out == execution/'entry-output' and not (out/'summary.json').exists()
            and mathematical_registration(arguments['panel'],spec,files,arguments['preparation'],
                    arguments['preparation_certificate_sha256']) == arguments == spec['arguments'],
            'prepared F5 target output reused or mathematical invocation changed')
    for key,value in arguments.items():
        write_immutable(out/(key+'.json'),value)
    write_immutable(out/'host.json',dict(system=platform.system(),machine=platform.machine(),
        os_release=platform.release(), scope='observed development platform; no hardware speedup claim'))
    preflight = meter(execution,'bin/worker',['--build-identity'],out,'build_identity',10)
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and read(out/'build_identity.stdout') == native['build_identity'],
            'prepared F5 native preflight substituted build identity')
    process = meter(execution,'bin/worker',[],out,'pipeline',
                    spec['runtime_watchdog_seconds']-30,stdin_argument='job')
    if process['timed_out']:
        status = 'NATIVE_TIMEOUT'
    elif process['returncode'] not in (0,2):
        status = 'NATIVE_ERROR'
    else:
        report = read(out/'pipeline.stdout')
        require((process['returncode'] == 0 and report['status'] == 'complete')
                or (process['returncode'] == 2 and report['status'] == 'incomplete'),
                'prepared F5 native terminal outcome differs')
        status = report['status'].upper()
    write_immutable(out/'summary.json',summary_for(arguments,status))
    return dict(status=status,run_id=arguments['seal']['run_id'])


def audit(execution, expected_spec):
    """Independent source, stdin, all target attempts and scalar replay; no solver."""
    execution = Path(execution)
    source_audit = audit_execution(execution,expected_spec)
    require(source_audit['entrypoint_succeeded'],'partial prepared F5 entrypoint remains a failure')
    files = check_extracted_assets(execution/'asset-files',expected_spec['asset_manifest'])
    build,rust_source,native = native_admission(files)
    arguments = expected_spec['arguments']
    require(arguments == mathematical_registration(arguments['panel'],expected_spec,files,
        arguments['preparation'],arguments['preparation_certificate_sha256']),
        'audited prepared F5 mathematical registration differs')
    root = execution/'entry-output'
    require(all(read(root/(key+'.json')) == value for key,value in arguments.items()),
            'prepared F5 executed invocation records differ')
    host = read(root/'host.json')
    require({k:host[k] for k in ('system','machine')} == native['platform']
            and type(host['os_release']) is str and host['os_release']
            and host['scope'] == 'observed development platform; no hardware speedup claim',
            'prepared F5 retained execution platform differs')
    preflight = audit_meter(execution,root,'build_identity',asset_role='bin/worker',
                            arguments=['--build-identity'],seconds=10)
    require(preflight['returncode'] == 0 and not preflight['timed_out']
            and read(root/'build_identity.stdout') == native['build_identity'],
            'prepared F5 retained preflight differs')
    process = audit_meter(execution,root,'pipeline',asset_role='bin/worker',arguments=[],
        seconds=expected_spec['runtime_watchdog_seconds']-30,stdin_argument='job')
    common = dict(**{k:arguments['seal'][k] for k in ('candidate_id','workload_id','run_id')},
        source_audit=source_audit, native_process=publication_process(process),
        raw_native_process_sha256=digest(root/'pipeline.metrics.json'),
        source_bound_execution_admitted=True, headline_online_admissible=False,
        fresh_paired_qualification=False, promotion_eligible=False, online_speedup=None)
    if process['timed_out'] or process['returncode'] not in (0,2):
        status = 'NATIVE_TIMEOUT' if process['timed_out'] else 'NATIVE_ERROR'
        require(read(root/'summary.json') == summary_for(arguments,status),
                'prepared F5 failure summary changed')
        return dict(common,status='AUDITED_'+status,scalar_verified=False,
            online_wall_ns=None,online_phases_ns=None,target_attempt_count=None,
            ordinary_queries_executed=None, scope='native failure retained; no verified target result')
    report = read(root/'pipeline.stdout')
    require(report['fixture'] == arguments['fixture']
            and ((process['returncode'] == 0 and report['status'] == 'complete')
                 or (process['returncode'] == 2 and report['status'] == 'incomplete'))
            and read(root/'summary.json') == summary_for(arguments,report['status'].upper()),
            'prepared F5 supplied fixture, status or summary changed')
    binding = verify_binding(report,build,rust_source,executable=execution/'asset-files/bin/worker')
    mathematical = audit_native_target(report,arguments['job'],arguments['preparation'],
                                       arguments['preparation_certificate_sha256'])
    require(mathematical['online_attempt_wall_ns'] <= process['native_wall_ns'],
            'prepared F5 native online interval exceeds its entire process')
    return dict(common,status='ADMITTED_'+mathematical['status']+'_PREPARED_F5_CONTROL',
        native_binding=binding,mathematical=mathematical,
        recovered_scalar=mathematical['recovered_scalar'],scalar_verified=mathematical['scalar_verified'],
        ordinary_queries_executed=0,target_attempt_count=report['solutions'][0]['trials'],
        online_wall_ns=mathematical['online_wall_ns'],online_phases_ns=mathematical['online_phases_ns'],
        scope='disclosed source-bound control; fresh comparison and calibration pending')


def frozen_audit(execution, expected_sha256, out):
    """Replay only using the preexecution-frozen auditor and interpreter."""
    execution,out = Path(execution).resolve(),Path(out).resolve()
    spec = read(execution/'execution.json')
    root = HERE.parents[1]
    require(root == execution/'extracted' and sha256(spec) == expected_sha256
            and sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode
            and interpreter_record() == spec['interpreter']
            and source_manifest(root) == spec['runtime_manifest'],
            'prepared F5 transport is outside its frozen source/interpreter')
    before = check_loaded_modules(root,spec['runtime_manifest'])
    write_immutable(out/'before.json',dict(binding=spec['binding'],loaded_modules=before))
    result = audit(execution,spec)
    after = check_loaded_modules(root,spec['runtime_manifest'])
    require(before.items() <= after.items() and source_manifest(root) == spec['runtime_manifest'],
            'prepared F5 audit source changed or dropped imports')
    write_immutable(out/'after.json',dict(binding=spec['binding'],loaded_modules=after))
    write_immutable(out/'admission.json',result)


def transport(execution, expected_sha256, out):
    execution,out = Path(execution).resolve(),Path(out).resolve()
    spec = read(execution/'execution.json')
    require(sha256(spec) == expected_sha256 and interpreter_record() == spec['interpreter']
            and not out.exists() and not out.is_relative_to(execution),
            'prepared F5 transport seal, interpreter or one-use output differs')
    audit_execution(execution,spec)
    out.mkdir(parents=True)
    script = execution/'extracted'/ 'research/ic_candidate_tournament_20260915/prepared_f5_runtime_v1.py'
    command = [sys.executable,'-I','-S','-B',str(script),'_frozen-audit',
               '--execution',str(execution),'--expected-execution-sha256',expected_sha256,'--out',str(out)]
    timed_out = False
    with (out/'stdout.txt').open('x') as stdout,(out/'stderr.txt').open('x') as stderr:
        try:
            cp = subprocess.run(command,env=frozen_environment(),stdout=stdout,stderr=stderr,
                                timeout=180,check=False)
            exit_code = cp.returncode
        except subprocess.TimeoutExpired:
            timed_out,exit_code = True,None
    receipt = dict(status='PASS_FROZEN_PREPARED_F5_TRANSPORT' if exit_code == 0 else 'REJECTED_TRANSPORT',
        exit_code=exit_code,timed_out=timed_out,execution_sha256=expected_sha256,interpreter_sha256=sha256(spec['interpreter']),
        stdout_sha256=digest(out/'stdout.txt'),stderr_sha256=digest(out/'stderr.txt'),
        native_solvers_executed=0,promotion_eligible=False,online_speedup=None)
    if exit_code == 0:
        before,after = read(out/'before.json'),read(out/'after.json')
        require(before['binding'] == after['binding'] == spec['binding']
                and before['loaded_modules'].items() <= after['loaded_modules'].items(),
                'prepared F5 frozen audit terminal gate differs')
        receipt.update(admission_sha256=digest(out/'admission.json'),
                       before_sha256=digest(out/'before.json'),after_sha256=digest(out/'after.json'))
    write_immutable(out/'transport.json',receipt)
    require(exit_code == 0,'prepared F5 frozen transport rejected; retain failure, no native retry')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    freeze = sub.add_parser('register')
    for flag in ('repository','assets','panel','certificate','out'):
        freeze.add_argument('--'+flag,type=Path,required=True)
    freeze.add_argument('--expected-certificate-sha256',required=True)
    launch = sub.add_parser('execute')
    launch.add_argument('--registration',type=Path,required=True)
    launch.add_argument('--out',type=Path,required=True)
    launch.add_argument('--expected-execution-sha256',required=True)
    for command in ('audit','_frozen-audit'):
        check = sub.add_parser(command)
        check.add_argument('--execution',type=Path,required=True)
        check.add_argument('--expected-execution-sha256',required=True)
        check.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.command == 'register':
        spec = register(args.repository,args.assets,read(args.panel),read(args.certificate),
                        args.expected_certificate_sha256,args.out)
        print(json.dumps(dict(execution_sha256=sha256(spec),**spec['arguments']['seal']),sort_keys=True))
    elif args.command == 'execute':
        spec = read(args.registration/'execution.json')
        require(sha256(spec) == args.expected_execution_sha256,'prepared F5 externally frozen invocation differs')
        print(json.dumps(execute(args.registration,args.out,expected_spec=spec,
                                 timeout_seconds=spec['runtime_watchdog_seconds']),sort_keys=True))
    elif args.command == 'audit':
        print(json.dumps(transport(args.execution,args.expected_execution_sha256,args.out),sort_keys=True))
    else:
        frozen_audit(args.execution,args.expected_execution_sha256,args.out)


if __name__ == '__main__':
    main()
