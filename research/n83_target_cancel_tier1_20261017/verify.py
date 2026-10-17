#!/usr/bin/env python3
"""Independently replay one target cancellation N83 Tier 1 stratum."""

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module('n83_finegrid_run', HERE / 'run.py')
prior = module('n83_finegrid_prior_verify',
               ROOT / 'research/n83_e2e_final_holdout_20261016/verify.py')
group = prior.group


def read(path):
    return json.loads(path.read_text())


def fixture_points():
    runner.fixture()
    encoded = [json.loads(line) for line in (HERE / 'targets.jsonl').read_text().splitlines()]
    answers = [json.loads(line) for line in (HERE / 'answers.jsonl').read_text().splitlines()]
    if len(encoded) != 24 or len(answers) != 24 or len({tuple(x) for x in encoded}) != 24:
        raise ValueError('twenty-four distinct sealed public points required')
    points, scalars = [], []
    for index, (item, answer) in enumerate(zip(encoded, answers)):
        seed = hashlib.sha256(runner.generator.DOMAIN +
                              index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.point_from_json(item)
        if (answer['fixture_index'] != index or
                answer['seed_sha256'] != seed.hex() or
                int(answer['scalar']) != scalar or
                not group.on_curve(point) or
                group.multiply(point, group.R) is not None or
                group.multiply(group.G, scalar) != point):
            raise ValueError(f'fresh point {index} fails scalar replay')
        points.append(point)
        scalars.append(scalar)
    return encoded, points, scalars


def folder(output, event):
    point_output = output / 'points' / f"p{event['point_index']:02d}"
    return point_output, point_output / 'runs' / event['name']


def replay_ic(output, event, point, scalar, by_name, field, curve):
    point_output, run_folder = folder(output, event)
    name = event['candidate']
    _, k, layout, expected_id = by_name[name]
    try:
        result = prior.dynamic.verify_ic(point_output, name, k, layout, point, scalar)
        result['candidate_id'] = expected_id
        if result['status'] == 'PASS_INDEPENDENT_REPLAY':
            summary = read(run_folder / 'stdout.jsonl')
            manifest = read(HERE / f'candidate_{expected_id}.json')
            source = runner.SOURCE_SHA256
            if (summary['target_threads'] != 2 or
                    summary['target_cancel'] is not (name == 'parallel_cancel') or
                    summary['root_table_indexed5_id30_enabled'] is not False or
                    manifest['candidate_id'] != expected_id or
                    manifest['factor_base']['ordered_base_dump_sha256'] !=
                    runner.parent.digest(run_folder / 'base_dump.jsonl') or
                    manifest['factor_base']['source_base_hash_blake3'] !=
                    summary['base_hash'] or
                    manifest['implementation']['source_sha256']['ic'] != source or
                    manifest['implementation']['algorithm_flags']['target_threads'] != 2 or
                    manifest['implementation']['algorithm_flags']['target_cancel'] is not
                    (name == 'parallel_cancel')):
                raise ValueError('cancellation mode, exact base or candidate differs')
        return result
    except Exception as error:
        return {'arm': name, 'k': k, 'layout': layout,
                'candidate_id': expected_id, 'status': 'VERIFICATION_FAILED',
                'error': str(error)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('--stratum', type=int, choices=range(4), required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    hashes, lines = runner.fixture()
    _, points, scalars = fixture_points()
    arms, choice = runner.choose()
    events = runner.schedule(args.stratum, arms)
    host, source, preflight = (read(output / name) for name in
                               ('host.json', 'source.json', 'preflight.json'))
    if (host['schema'] != 'n83.target-cancel-tier1-host/v1' or
            host['git_sha'] is None or host['stratum'] != args.stratum or
            host['instance_id'] !=
            f"{host['parent_instance_id']}-cancel-stratum{args.stratum}" or
            preflight['status'] != 'PASS_CAPACITY_PREFLIGHT' or
            source['schema'] != 'n83.target-cancel-tier1-source/v1' or
            source['hashes'] != hashes or
            source['source_sha256'] != runner.SOURCE_SHA256 or
            source['runner_sha256'] != runner.parent.digest(HERE / 'run.py') or
            source['targets_sha256'] != runner.TARGETS_SHA256 or
            source['answers_sha256'] != runner.ANSWERS_SHA256 or
            source['fixture_manifest_sha256'] !=
            runner.parent.digest(HERE / 'manifest.json') or
            source['stratum'] != args.stratum or
            source['arms'] != [list(arm) for arm in arms] or
            source['choice'] != choice or source['events'] != events or
            source['rho_seed_rule'] !=
            '2026110300+2*point_index+fixed_candidate_index' or
            source['calibration_seed'] != int(runner.cold.RHO_SEED) or
            source['calibration_target_sha256'] !=
            runner.parent.digest(runner.cold.TARGETS) or
            source['calibration_answer_sha256'] !=
            runner.parent.digest(runner.cold.ANSWERS) or
            source['global_producer_cap_seconds'] !=
            runner.GLOBAL_PRODUCER_CAP_SECONDS):
        raise ValueError('source, host, fixture or selected holdout arms differ')
    binaries_match = (
        runner.parent.digest(output / 'executed_bins/koblitz_orbit_dlp_fast') ==
        source['binary_sha256'] and
        runner.parent.digest(output / 'executed_bins/n83_rho_optimization_eval') ==
        source['rho_binary_sha256'])
    topology = host.get('topology', {})
    topology_status = ('PASS' if topology.get('lscpu_json') and
                       len(topology.get('visible_cpu_topology', {})) >= 4 and
                       (output / 'host_telemetry_after.json').exists() else 'MISSING')
    finished, telemetry, target_checks, launch_checks = [], [], [], []
    for event in events:
        point_output, run_folder = folder(output, event)
        usage_path = run_folder / 'usage.json'
        if usage_path.exists():
            finished.append(datetime.fromisoformat(read(usage_path)['finished_utc']))
        telemetry.append(prior.prior.telemetry_receipt(point_output, event['name']))
        target = point_output / 'target.jsonl'
        target_checks.append(target.exists() and
                             target.read_text() == lines[event['point_index']])
        try:
            launch = read(run_folder / 'launch.json')
            command = launch['command']
            if event['kind'] in ('rho', 'calibration'):
                expected_target = (runner.cold.TARGETS if
                                   event['kind'] == 'calibration' else target)
                seed = (int(runner.cold.RHO_SEED) if
                        event['kind'] == 'calibration' else
                        runner.PAIR_RHO_SEED_BASE + len(arms) * event['point_index'] +
                        next(i for i, arm in enumerate(arms)
                             if arm[0] == event['candidate']))
                valid = (prior.suffix(command[1], 3) ==
                         prior.suffix(expected_target, 3) and
                         command[2] == str(seed) and
                         prior.suffix(command[3], 3) ==
                         prior.suffix(runner.cold.RHO_CONFIG, 3) and
                         prior.suffix(command[4], 5) ==
                         prior.suffix(run_folder / 'receipt.json', 5) and
                         launch['environment_flags'] == {} and
                         launch['cap_seconds'] == runner.RHO_CAP_SECONDS and
                         Path(command[0]).parts[-2:] ==
                         ('executed_bins', 'n83_rho_optimization_eval'))
            else:
                arm = next(item for item in arms if item[0] == event['candidate'])
                _, k, layout, _ = arm
                flags = runner.ic_environment(run_folder, layout, event['candidate'])
                flags = {key: value for key, value in flags.items()
                         if key.startswith('KIC_')}
                observed = dict(launch['environment_flags'])
                observed_dump = observed.pop('KIC_DUMP_BASE', '')
                expected_dump = flags.pop('KIC_DUMP_BASE')
                cap = runner.IC_CAP_SECONDS
                valid = (command[1] == f'construct:83:1:{k}' and
                         prior.suffix(command[2], 3) == prior.suffix(target, 3) and
                         command[3] == '29' and
                         prior.suffix(command[4], 5) ==
                         prior.suffix(run_folder / 'targets_result.jsonl', 5) and
                         observed == flags and
                         prior.suffix(observed_dump, 5) ==
                         prior.suffix(expected_dump, 5) and
                         launch['cap_seconds'] == cap and
                         Path(command[0]).parts[-2:] ==
                         ('executed_bins', 'koblitz_orbit_dlp_fast'))
            launch_checks.append(valid)
        except (KeyError, IndexError, OSError, ValueError, TypeError, StopIteration):
            launch_checks.append(False)
    order_status = ('PASS' if len(finished) == len(events) and
                    all(a < b for a, b in zip(finished, finished[1:])) else
                    'INCOMPLETE_OR_OUT_OF_ORDER')
    telemetry_status = ('PASS' if topology_status == 'PASS' and
                        all(item['status'] == 'PASS' and
                            item['steal_jiffies_delta'] == 0 for item in telemetry)
                        else 'INCOMPLETE_OR_STEAL')
    cal_encoded = json.loads(runner.cold.TARGETS.read_text())
    cal_answer = json.loads(runner.cold.ANSWERS.read_text())
    cal_point = group.point_from_json(cal_encoded)
    cal_scalar = int(cal_answer['scalar'])
    if (cal_answer['fixture_index'] != 0 or cal_answer['point'] != cal_encoded or
            group.multiply(group.G, cal_scalar) != cal_point or
            group.multiply(cal_point, group.R) is not None):
        raise ValueError('calibration point fails scalar replay')
    field, curve = prior.finegrid.panel.curve_record()
    by_name = {arm[0]: arm for arm in arms}
    point_rows = {index: {'point_index': index, 'arms': [],
                          'rho_by_candidate': {}, 'calibration_by_candidate': {}}
                  for index in range(args.stratum * 6, args.stratum * 6 + 6)}
    calibrations = {}
    for event in events:
        index, candidate, kind = (event[key] for key in
                                  ('point_index', 'candidate', 'kind'))
        candidate_id = by_name[candidate][3]
        if kind == 'ic':
            point_rows[index]['arms'].append(replay_ic(
                output, event, points[index], scalars[index], by_name,
                field, curve))
        elif kind == 'rho':
            seed = (runner.PAIR_RHO_SEED_BASE + len(arms) * index +
                    next(i for i, arm in enumerate(arms) if arm[0] == candidate))
            point_rows[index]['rho_by_candidate'][candidate_id] = prior.replay_rho(
                output, event, points[index], scalars[index], seed,
                output / 'points' / f'p{index:02d}' / 'target.jsonl', hashes)
        else:
            calibrations[(index, candidate, event['name'].startswith('cal_before_'))] = (
                prior.replay_rho(output, event, cal_point, cal_scalar,
                                 int(runner.cold.RHO_SEED), runner.cold.TARGETS,
                                 hashes))
    for index, point in point_rows.items():
        for name, _, _, candidate_id in arms:
            point['calibration_by_candidate'][candidate_id] = prior.calibration(
                calibrations.get((index, name, True)),
                calibrations.get((index, name, False)))
    ic = [row for point in point_rows.values() for row in point['arms']]
    rho = [row for point in point_rows.values() for row in
           point['rho_by_candidate'].values()]
    cal = list(calibrations.values())
    base_status = ('PASS' if all(len({row.get('base_hash') for row in ic
                                  if row.get('candidate_id') == candidate_id and
                                  row['status'] == 'PASS_INDEPENDENT_REPLAY'}) == 1
                                 for _, _, _, candidate_id in arms) else 'MISMATCH')
    replay_status = ('PASS_INDEPENDENT_REPLAY' if binaries_match and
                     all(target_checks) and all(launch_checks) and
                     len(ic) == 6 * len(arms) and len(rho) == 6 * len(arms) and
                     len(cal) == 12 * len(arms) and
                     all(item['status'] == 'PASS_INDEPENDENT_REPLAY'
                         for item in ic + rho + cal) and base_status == 'PASS'
                     else 'VERIFICATION_FAILED')
    resource_status = ('PASS' if replay_status == 'PASS_INDEPENDENT_REPLAY' and
                       all(item['resources']['status'] == 'PASS'
                           for item in ic + rho + cal) else
                       'FAILED_OR_INCOMPLETE')
    calibration_status = ('PASS_10_PERCENT_PER_PAIR' if all(
        len(point['calibration_by_candidate']) == len(arms) and
        all(item['status'] == 'PASS_10_PERCENT' for item in
            point['calibration_by_candidate'].values())
        for point in point_rows.values()) else 'NOISY_OR_INCOMPLETE')
    admitted = (replay_status == 'PASS_INDEPENDENT_REPLAY' and
                resource_status == 'PASS' and order_status == 'PASS' and
                telemetry_status == 'PASS' and
                calibration_status == 'PASS_10_PERCENT_PER_PAIR' and
                not (output / 'producer_stop.json').exists())
    validated = []
    for index in sorted(point_rows):
        point = point_rows[index]
        validated.append({
            'point_index': index,
            'arms': [{'candidate_id': item.get('candidate_id'),
                      'status': item['status'],
                      'resource_status': item.get('resources', {}).get('status'),
                      'one_target_online_ms': item.get('one_target_online_ms'),
                      'cold_process_ms': item.get('cold_process_ms_exploratory'),
                      'raw_replay': item} for item in point['arms']],
            'rho_by_candidate': {
                candidate_id: {'status': item['status'],
                               'resource_status': item.get('resources', {}).get('status'),
                               'online_ms': item.get('online_ms'),
                               'raw_replay': item}
                for candidate_id, item in point['rho_by_candidate'].items()},
            'calibration_by_candidate': point['calibration_by_candidate'],
        })
    analysis = {
        'schema': 'n83.target-cancel-tier1-stratum-analysis/v1',
        'stratum': args.stratum,
        'point_indices': list(sorted(point_rows)),
        'source_sha256': runner.SOURCE_SHA256,
        'targets_sha256': runner.TARGETS_SHA256,
        'answers_sha256': runner.ANSWERS_SHA256,
        'choice': choice,
        'status': 'STRATUM_ADMITTED' if admitted else 'BOUNDED_UNKNOWN',
        'tier_status': 'PASS_TIER1' if admitted else 'NOT_ADMITTED',
        'hardware_tier': 'cloud_instance_tier1' if admitted else
                         'cloud_instance_exploratory',
        'evidence_status': replay_status,
        'resource_status': resource_status,
        'aa_status': calibration_status,
        'base_status': base_status,
        'execution_order_status': order_status,
        'topology_status': topology_status,
        'telemetry_status': telemetry_status,
        'host_telemetry_by_event': telemetry,
        'producer_stop': (read(output / 'producer_stop.json') if
                          (output / 'producer_stop.json').exists() else None),
        'all_per_point_targets_match': all(target_checks),
        'all_launch_commands_match': all(launch_checks),
        'archived_executed_binaries_match': binaries_match,
        'points': validated,
    }
    (output / 'analysis.json').write_text(json.dumps(analysis, indent=2,
                                                   sort_keys=True) + '\n')
    print(json.dumps(analysis, sort_keys=True))
    return 0 if admitted else 2


if __name__ == '__main__':
    raise SystemExit(main())
