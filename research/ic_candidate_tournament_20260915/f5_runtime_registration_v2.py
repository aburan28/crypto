"""Freeze one supplied n17 F4/F5 development solve, including all controllers."""
import argparse
import copy
import json
from pathlib import Path

from f5_runtime_inputs_v1 import native_admission
from generic_admission import method_record
from generic_stages import DEFAULTS, effective_config
from generic_solver_feasibility import assess as assess_layout
from identity import candidate_manifest, run_id, sha256, workload_manifest
from oracle import require
from sat_runtime_execution_v3 import register as register_runtime
from static_sat_assets_v3 import verified_assets

SETTINGS = {'question', 'resources', 'target_input', 'algorithm_seed', 'config', 'run_number'}
CONFIG = {'solver', 'linear_algebra', 'summands', 'groebner_degree', 'node_budget',
          'conflict_budget', 'batch_trials', 'max_trials'}
INPUT_LAW = 'one-supplied-public-point; seed-is-provenance'


def validate_panel(panel, curve):
    require(set(panel) == SETTINGS and panel['question'] == 'development-source-control',
            'F5 v2 only admits a declared development control; paired adapter pending')
    require(type(panel['algorithm_seed']) is int and 0 <= panel['algorithm_seed'] < 2**64
            and type(panel['run_number']) is int and panel['run_number'] >= 0,
            'F5 v2 seed/run outside registered domain')
    cfg = panel['config']
    require(type(cfg) is dict and set(cfg) == CONFIG and cfg['solver'] in ('f4', 'f5')
            and cfg['linear_algebra'] == 'dense' and cfg['summands'] == 3
            and cfg['groebner_degree'] == 3,
            'F5 v2 requires the admitted three-summand degree-three dense pipeline')
    require(all(type(cfg[key]) is int and 0 < cfg[key] < 2**64 for key in (
        'node_budget', 'conflict_budget', 'batch_trials', 'max_trials')),
        'F5 v2 pipeline limits must fit positive native 64-bit integers')
    resources = panel['resources']
    require(type(resources) is dict and resources == dict(host_class='physical-macos-arm64',
            cpu_workers=1, target_count=1, memory_limit_bytes=None,
            total_wall_limit_seconds=resources.get('total_wall_limit_seconds'))
            and type(resources['cpu_workers']) is int
            and type(resources['target_count']) is int
            and type(resources['total_wall_limit_seconds']) is int
            and resources['total_wall_limit_seconds'] > 30,
            'F5 v2 requires a bounded one-worker development envelope')
    target_input = panel['target_input']
    require(type(target_input) is dict and set(target_input) == {
        'point', 'seed', 'input_law', 'point_was_previously_supplied', 'known_scalar_supplied'}
        and target_input['input_law'] == INPUT_LAW and target_input['seed'] is None
        and target_input['point_was_previously_supplied'] is True
        and target_input['known_scalar_supplied'] is False,
        'F5 v2 requires a disclosed supplied point without a scalar')
    target = curve.decode(target_input['point'])
    require(target is not None and curve.mul(target, curve.r) is None,
            'F5 v2 target is outside the admitted subgroup')
    return target


def mathematical_registration(panel, spec, files, *, check_host=True):
    fixture, inventory, curve, base, build, rust_source, native = native_admission(files, check_host=check_host)
    target = validate_panel(panel, curve)
    require(spec['entrypoint'] == dict(module='f5_runtime_pipeline_v2', callable='run')
            and spec['runtime_watchdog_seconds'] == panel['resources']['total_wall_limit_seconds'],
            'F5 mathematical registration differs from execution envelope')
    job = dict(mode='ic', degree=17, curve_a=1,
        public_targets=[[str(value) for value in target]], target_seeds=[],
        factor_base=dict(kind='standard_subspace', dimension=6), exclusive_phases=True,
        algorithm_seed=panel['algorithm_seed'], config=copy.deepcopy(panel['config']))
    layout = assess_layout(dict(cells=['n17a1'], candidates=[dict(id='declared-development-arm',
        config=dict(job['config'], factor_base=job['factor_base']))]))
    require(layout['status'] == 'PASS_STATIC_LAYOUT_ONLY', 'F5 registered encoder exceeds layout cap')
    declared = copy.deepcopy(inventory)
    declared['effective_config'] = dict(copy.deepcopy(DEFAULTS), **job['config'])
    effective_config(job, declared)
    native_method = method_record(job, declared, {'base':base}, build,
                                  collector_plan=inventory['collector_dispatch'])
    method = copy.deepcopy(native_method)
    source = dict(schema_version=2, execution_binding=spec['binding'], native=native)
    method['implementation']['components'].append(dict(role='complete-frozen-Python-interpreter-native-binding',
                                                       sha256=sha256(source)))
    method['implementation']['flags'].update(execution_binding=spec['binding'],
        stdin_policy='canonical-registered-job-UTF8; empty-input-seeds-for-unseeded-public-point; fixture-null-provenance',
        native_watchdog_group='inherit-controller-group-no-native-fork-or-setsid',
        native_thread_environment='all-listed-pools-one-thread',
        online_interval='native first target-dependent work through general scalar replay',
        external_audit='independent Python audit outside native online interval')
    candidate = candidate_manifest(fixture, declared, method)
    target_fixture = dict(fixture, targets=job['public_targets'], target_seeds=[None])
    workload = workload_manifest(target_fixture, input_law=INPUT_LAW,
        algorithm_seed=panel['algorithm_seed'], resource_envelope=panel['resources'], cache_policy='cold')
    seal = dict(schema_version=2, registration_stage='before-execution',
        panel_sha256=sha256(panel), source_sha256=sha256(source), method_sha256=sha256(method),
        candidate_sha256=sha256(candidate), workload_sha256=sha256(workload),
        candidate_id=candidate['candidate_id'], workload_id=workload['workload_id'],
        run_id=run_id(candidate['candidate_id'], workload['workload_id'], panel['run_number']))
    return dict(panel=copy.deepcopy(panel), job=job, fixture=target_fixture, static_layout=layout,
                source=source, native_method=native_method, method=method,
                candidate=candidate, workload=workload, seal=seal)


def register(repository, assets, panel, output):
    assets = Path(assets)
    manifest, seal = [json.loads((assets/name).read_text()) for name in ('manifest.json', 'seal.json')]
    files = verified_assets(assets, manifest, seal)
    native_admission(files)
    return register_runtime(repository, output, module='f5_runtime_pipeline_v2', action='run',
        arguments=None, arguments_factory=lambda spec: mathematical_registration(panel, spec, files),
        timeout_seconds=panel['resources']['total_wall_limit_seconds'], asset_snapshot=assets)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--assets', type=Path, required=True)
    parser.add_argument('--panel', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    spec = register(args.repository, args.assets, json.loads(args.panel.read_text()), args.out)
    print(json.dumps(dict(execution_sha256=sha256(spec), **spec['arguments']['seal']), sort_keys=True))
