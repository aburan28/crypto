"""Freeze the complete SAT method and one supplied target before execution."""
import argparse
import copy
import json
from pathlib import Path

from identity import candidate_manifest, run_id, sha256, workload_manifest
from oracle import require
from sat_runtime_execution_v3 import register as register_runtime
from static_sat_assets_v3 import verified_assets
from static_sat_inputs_v3 import native_admission
from static_sat_registration_v2 import method_record

SETTINGS = {'cms_conflict_budget', 'cms_timeout_seconds', 'export_timeout_seconds',
            'max_relation_queries', 'max_descent_queries', 'relation_query_seed',
            'descent_query_seed', 'export_nonce', 'target_input', 'resources',
            'question', 'run_number'}


def validate_panel(panel, curve):
    require(set(panel) == SETTINGS, 'SAT v3 settings missing or extra')
    for key in SETTINGS-{'target_input', 'resources', 'question'}:
        require(type(panel[key]) is int and panel[key] >= 0,
                'SAT v3 integer setting invalid: '+key)
    require(all(panel[key] > 0 for key in ('cms_conflict_budget', 'cms_timeout_seconds',
                                         'export_timeout_seconds', 'max_relation_queries',
                                         'max_descent_queries')),
            'SAT v3 budgets must be positive')
    require(all(panel[key] < 2**64 for key in (
                'relation_query_seed', 'descent_query_seed', 'export_nonce')),
            'SAT v3 seed outside the registered RNG domain')
    require(panel['question'] in ('development-source-control', 'fresh-paired-qualification'),
            'SAT v3 scientific question unspecified')
    target_input = panel['target_input']
    require(set(target_input) == {'point', 'seed', 'input_law',
                                  'point_was_previously_supplied', 'known_scalar_supplied'}
            and target_input['known_scalar_supplied'] is False
            and type(target_input['point_was_previously_supplied']) is bool
            and type(target_input['input_law']) is str and target_input['input_law'],
            'SAT v3 supplied target input malformed')
    target = curve.decode(target_input['point'])
    require(target is not None and curve.mul(target, curve.r) is None,
            'SAT v3 target leaves the prime subgroup')
    # This version has no paired-arm/calibration/exposure ledger adapter yet.
    # A point flag by itself cannot prove freshness or scientific qualification.
    require(panel['question'] == 'development-source-control',
            'fresh paired qualification requires the reviewed campaign adapter')
    resources = panel['resources']
    require(resources == dict(host_class='physical-macos-arm64', cpu_workers=1,
                              target_count=1, memory_limit_bytes=None,
                              total_wall_limit_seconds=resources.get('total_wall_limit_seconds'))
            and type(resources['total_wall_limit_seconds']) is int
            and resources['total_wall_limit_seconds'] > 0,
            'SAT v3 requires a one-worker registered development envelope')
    return target


def mathematical_registration(panel, spec, files):
    fixture, report, curve, base, native = native_admission(files)
    target = validate_panel(panel, curve)
    require(spec['entrypoint'] == dict(module='static_sat_pipeline_v3', callable='run')
            and spec['runtime_watchdog_seconds'] == panel['resources']['total_wall_limit_seconds'],
            'SAT mathematical registration differs from execution envelope')
    source = dict(schema_version=3, execution_binding=spec['binding'], native=native)
    legacy_parameters = dict(panel, cms_max_models_per_query=1,
                             cms_executable_sha256=native['cms_binary_sha256'],
                             cms_build_receipt_sha256=native['cms_build_receipt_sha256'],
                             exporter_source_sha256=native['exporter_source_sha256'],
                             source_encoding='wide-symmetrised-S4-circuit-XOR-DIMACS')
    method = method_record(legacy_parameters, sha256(source), fixture)
    method['implementation']['flags'].update(
        execution_binding=spec['binding'], exporter_binary_sha256=native['exporter_binary_sha256'],
        native_thread_environment='all-listed-pools-one-thread; CMS-one-thread',
        native_watchdog_group='inherit-controller-group-no-native-fork-or-setsid',
        online_observer_cost='all-target-native-wrapper-checks-charged-to-PDP',
        exporter_nonce_policy='run-registered-nonce')
    candidate = candidate_manifest(fixture, report, method)
    target_fixture = dict(fixture, targets=[list(target)],
                          target_seeds=[panel['target_input']['seed']])
    workload = workload_manifest(
        target_fixture, input_law=panel['target_input']['input_law'],
        algorithm_seed=panel['relation_query_seed'],
        resource_envelope=panel['resources'], cache_policy='warm')
    # Both RNG streams and exporter seed must identify a measured workload.
    workload['record'].update(algorithm_seeds=dict(
        collection=panel['relation_query_seed'], descent=panel['descent_query_seed'],
        exporter=panel['export_nonce']), question=panel['question'],
        point_was_previously_supplied=panel['target_input']['point_was_previously_supplied'])
    workload['record_sha256'] = sha256(workload['record'])
    workload['workload_id'] = workload['record_sha256'][:12]
    seal = dict(schema_version=3, panel_sha256=sha256(panel), source_sha256=sha256(source),
                method_sha256=sha256(method), candidate_sha256=sha256(candidate),
                workload_sha256=sha256(workload), candidate_id=candidate['candidate_id'],
                workload_id=workload['workload_id'],
                run_id=run_id(candidate['candidate_id'], workload['workload_id'], panel['run_number']),
                registration_stage='before-execution')
    return dict(panel=copy.deepcopy(panel), source=source, method=method,
                candidate=candidate, workload=workload, seal=seal)


def register(repository, assets, panel, output):
    assets = Path(assets)
    manifest = json.loads((assets/'manifest.json').read_text())
    seal = json.loads((assets/'seal.json').read_text())
    files = verified_assets(assets, manifest, seal)
    native_admission(files)
    return register_runtime(repository, output, module='static_sat_pipeline_v3', action='run',
                            arguments=None, arguments_factory=lambda spec:
                            mathematical_registration(panel, spec, files),
                            timeout_seconds=panel['resources']['total_wall_limit_seconds'],
                            asset_snapshot=assets)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    spec = register(args.repository, args.assets, json.loads(args.panel.read_text()), args.out)
    print(json.dumps(dict(execution_sha256=sha256(spec), **spec['arguments']['seal']), sort_keys=True))
