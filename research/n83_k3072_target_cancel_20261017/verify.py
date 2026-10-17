#!/usr/bin/env python3
"""Independently replay the K3072 cancellation comparison."""

import argparse
import hashlib
import importlib.util
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location(
    'n83_parallel_cancel_group',
    ROOT / 'research/n83_indexed5_20261013/analyze.py')
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)
spec = importlib.util.spec_from_file_location(
    'n83_parallel_cancel_runner', HERE / 'run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def read(path):
    return json.loads(path.read_text())


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'existing analysis differs: {path}')
    path.write_text(blob)


def fixture():
    targets = [group.point_from_json(json.loads(line)) for line in
               runner.TARGETS.read_text().splitlines()]
    answers = [json.loads(line) for line in
               runner.ANSWERS.read_text().splitlines()]
    if len(targets) != 128 or len(answers) != 128 or len(set(targets)) != 128:
        raise ValueError('fixture count or distinctness differs')
    for index, (point, answer) in enumerate(zip(targets, answers)):
        seed = hashlib.sha256(
            b'n83-root11-finegrid-batch-20261017-v1' +
            index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        if (answer['fixture_index'] != index or
                answer['seed_sha256'] != seed.hex() or
                int(answer['scalar']) != scalar or
                group.multiply(group.G, scalar) != point):
            raise ValueError(f'fixture {index} differs')
    return targets, answers


def replay_one(policy, folder, plan, targets, answers):
    receipt = read(folder / 'run_receipt.json')
    command = read(folder / 'command.json')
    threads, cancel = runner.POLICIES[policy]
    manifest = read(HERE / f'candidate_{policy}.json')
    ids = read(HERE / 'candidate_ids.json')
    candidate_id = ids[policy]['candidate_id']
    expected_environment = runner.dynamic.environment(folder, 'root11_dynamic')
    expected_environment['KIC_TARGET_THREADS'] = str(threads)
    expected_environment['KIC_TARGET_CANCEL'] = cancel
    canonical = json.dumps({key: value for key, value in manifest.items()
                            if key != 'candidate_id'}, sort_keys=True,
                           separators=(',', ':'), ensure_ascii=False).encode()
    expected_id = ('IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0h'
                   + hashlib.sha256(canonical).hexdigest()[:12])
    if (candidate_id != expected_id or manifest['candidate_id'] != expected_id or
            ids[policy]['manifest_sha256'] !=
            runner.digest(HERE / f'candidate_{policy}.json') or
            manifest['implementation']['source_sha256']['ic'] !=
            plan['sha256']['source'] or
            manifest['implementation']['algorithm_flags']['target_threads'] != threads or
            manifest['implementation']['algorithm_flags']['target_cancel'] !=
            (cancel == '1') or
            receipt['candidate_id'] != candidate_id or
            command['candidate_id'] != candidate_id or
            receipt['policy'] != policy or
            receipt['producer_status'] != 'EXIT_0' or
            receipt['exit_code'] != 0 or
            receipt['source_sha256'] != plan['sha256']['source'] or
            receipt['binary_sha256'] != plan['sha256']['binary'] or
            command['policy'] != policy or
            command['plan_sha256'] != runner.digest(HERE / 'plan.json') or
            {key: value for key, value in command['environment'].items()
             if key.startswith('KIC_')} !=
            {key: value for key, value in expected_environment.items()
             if key.startswith('KIC_')} or
            command['argv'] != [str((Path(plan['binary_path'])).resolve()),
                                'construct:83:1:3072', str(runner.TARGETS),
                                '29', str(folder / 'targets_result.jsonl')] or
            receipt['stdout_sha256'] != runner.digest(folder / 'stdout.jsonl') or
            receipt['stderr_sha256'] != runner.digest(folder / 'stderr.txt') or
            (folder / 'stderr.txt').read_text() or
            receipt['target_result_sha256'] !=
                runner.digest(folder / 'targets_result.jsonl') or
            receipt['base_dump_sha256'] !=
                runner.digest(folder / 'base_dump.jsonl')):
        raise ValueError('producer, source, command or raw digest differs')
    summary = read(folder / 'stdout.jsonl')
    if (summary != receipt['summary'] or summary['n'] != 83 or
            summary['a'] != 1 or summary['orbit_columns'] != 3072 or
            summary['factor_base_points'] != 166 * 3072 or
            summary['rank'] != 3072 or summary['targets'] != 128 or
            summary['targets_solved'] != 128 or summary['targets_failed'] != 0 or
            summary['target_threads'] != threads or
            summary['target_cancel'] != (cancel == '1') or
            summary['base_hash'] !=
            manifest['factor_base']['source_base_hash_blake3']):
        raise ValueError('summary, rank, target coverage or policy differs')
    base = group.replay_base(read(folder / 'base_dump.jsonl'), summary)
    records = [json.loads(line) for line in
               (folder / 'targets_result.jsonl').read_text().splitlines()]
    if len(base) != 166 * 3072 or len(records) != 128:
        raise ValueError('base or target record count differs')
    times, probes = [], []
    for index, (record, point, answer) in enumerate(zip(records, targets, answers)):
        indices = record['point_indices']
        if not isinstance(indices, list) or len(indices) != 4:
            raise ValueError(f'point {index}: incomplete four-point relation')
        relation = None
        for base_index in indices:
            if not isinstance(base_index, int) or not 0 <= base_index < len(base):
                raise ValueError(f'point {index}: base index out of range')
            relation = group.add(relation, base[base_index])
        recovered = int(record['recovered_scalar'])
        interval = float(record['target_ms'])
        phase = sum(float(record[key]) for key in
                    ('target_query_ms', 'target_pdp_ms',
                     'target_relation_check_ms', 'target_descent_ms',
                     'target_recovery_check_ms'))
        if (record['fixture_index'] != index or
                group.point_from_json(record['target']) != point or
                recovered != int(answer['scalar']) or
                group.multiply(group.G, recovered) != point or
                relation != point or record['group_verified'] is not True or
                not math.isfinite(interval) or interval <= 0 or
                abs(interval - phase) > 1.0 or
                not isinstance(record['probes'], int) or record['probes'] <= 0):
            raise ValueError(f'point {index}: scalar, relation or phase differs')
        times.append(interval)
        probes.append(record['probes'])
    return {'status': 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY',
            'candidate_id': candidate_id,
            'policy': policy, 'target_count': len(times),
            'base_hash': summary['base_hash'],
            'actual_usable_points': len(base),
            'effective_columns': summary['orbit_columns'],
            'online_geometric_mean_ms':
                math.exp(statistics.mean(math.log(value) for value in times)),
            'online_median_ms': statistics.median(times),
            'cold_process_ms': summary['timing_ms']['process_total'],
            'per_target_online_ms': times,
            'per_target_winning_worker_probes': probes,
            'raw_receipt_sha256': runner.digest(folder / 'run_receipt.json')}


def verify(root):
    plan = runner.frozen(Path(read(HERE / 'plan.json')['binary_path']))
    targets, answers = fixture()
    rows = []
    for policy in runner.POLICIES:
        folder = root / policy
        try:
            rows.append(replay_one(policy, folder, plan, targets, answers))
        except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
            rows.append({'status': 'VERIFICATION_FAILED_OR_CENSORED',
                         'policy': policy, 'error': str(error)})
    complete = (all(row['status'] ==
                    'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'
                    for row in rows) and
                len({row['base_hash'] for row in rows}) == 1)
    result = {'schema': 'n83.k3072-target-cancel-analysis/v1',
              'status': 'PASS_FULL_INDEPENDENT_REPLAY' if complete else
                        'BOUNDED_UNKNOWN',
              'plan_sha256': runner.digest(HERE / 'plan.json'),
              'source_sha256': plan['sha256']['source'],
              'binary_sha256': plan['sha256']['binary'],
              'target_sha256': plan['sha256']['targets'],
              'rows': rows}
    if complete:
        by_policy = {row['policy']: row for row in rows}
        controls = by_policy['serial']['per_target_online_ms']
        for name in ('parallel_cancel',):
            values = by_policy[name]['per_target_online_ms']
            differences = [math.log(a / b) for a, b in zip(controls, values)]
            mean = statistics.mean(differences)
            halfwidth = 1.979 * statistics.stdev(differences) / math.sqrt(128)
            result[name + '_vs_serial'] = {
                'paired_geometric_speedup': math.exp(mean),
                'descriptive_95_interval': [math.exp(mean - halfwidth),
                                            math.exp(mean + halfwidth)],
                'faster_point_count': sum(b < a for a, b in zip(controls, values))}
    write_once(root / 'analysis.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    args = parser.parse_args()
    result = verify(args.root.resolve())
    print(json.dumps({'status': result['status'],
                      'rows': [(row['policy'], row['status'])
                               for row in result['rows']]}, sort_keys=True))


if __name__ == '__main__':
    main()
