#!/usr/bin/env python3
"""Independently verify each high-K N83 base, relation, scalar and timing."""
import argparse
import hashlib
import importlib.util
import json
import math
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

group = module('n83_highk_group', ROOT / 'research/n83_indexed5_20261013/analyze.py')


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fixture():
    public_path = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl'
    answer_path = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_answers.jsonl'
    public = [group.point_from_json(json.loads(line)) for line in public_path.read_text().splitlines()]
    answers = [json.loads(line) for line in answer_path.read_text().splitlines()]
    if len(public) != 128 or len(answers) != 128 or len(set(public)) != 128:
        raise ValueError('fixture count or distinctness differs')
    for index, (point, answer) in enumerate(zip(public, answers)):
        seed = hashlib.sha256(b'n83-root11-finegrid-batch-20261017-v1' +
                              index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        if (answer['fixture_index'] != index or
                answer['seed_sha256'] != seed.hex() or
                int(answer['scalar']) != scalar or
                group.multiply(group.G, scalar) != point or
                group.multiply(point, group.R) is not None):
            raise ValueError(f'fixture {index} differs')
    return public, answers


def verify_one(folder, k, plan):
    import run
    if folder.name != f'k{k}' or k not in run.ARMS:
        raise ValueError('unexpected arm path')
    command = json.loads((folder / 'command.json').read_text())
    receipt = json.loads((folder / 'run_receipt.json').read_text())
    if (command['k'] != k or receipt['k'] != k or
            command['plan_sha256'] != digest(HERE / 'plan.json') or
            command['cap_seconds'] != plan['cap_seconds_per_arm'] or
            command['argv'] != [plan['binary_path'], f'construct:83:1:{k}',
                                str(run.FIXTURE), '29',
                                str(folder / 'targets_result.jsonl')] or
            {key: value for key, value in command['environment'].items()
             if key.startswith('KIC_')} !=
            {key: value for key, value in run.environment(folder).items()
             if key.startswith('KIC_')} or
            not command['environment'].get('PATH') or
            receipt['producer_status'] != 'EXIT_0' or receipt['exit_code'] != 0 or
            receipt['stdout_sha256'] != digest(folder / 'stdout.jsonl') or
            receipt['stderr_sha256'] != digest(folder / 'stderr.txt') or
            receipt['base_dump_sha256'] != digest(folder / 'base_dump.jsonl') or
            receipt['target_results_sha256'] != digest(folder / 'targets_result.jsonl')):
        raise ValueError('command, producer status or raw receipt differs')
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    if (summary != receipt['summary'] or summary['n'] != 83 or
            summary['a'] != 1 or summary['orbit_columns'] != k or
            summary['factor_base_points'] != 166*k or
            summary['rank'] != k or
            summary['root_table_indexed5_id30_enabled'] is not False or
            summary['root_only11_enabled'] is not True or
            summary['targets'] != 128 or
            summary['targets_solved'] != 128 or summary['targets_failed'] != 0):
        raise ValueError('summary, base size, rank or target coverage differs')
    base = group.replay_base(json.loads((folder / 'base_dump.jsonl').read_text()), summary)
    rows = [json.loads(line) for line in
            (folder / 'targets_result.jsonl').read_text().splitlines()]
    points, answers = fixture()
    if len(base) != 166*k or len(rows) != 128:
        raise ValueError('base or target row count differs')
    times, probes = [], []
    for index, (row, point, answer) in enumerate(zip(rows, points, answers)):
        relation = None
        indices = row['point_indices']
        if not isinstance(indices, list) or len(indices) != 4:
            raise ValueError(f'point {index}: incomplete four-point relation')
        for base_index in indices:
            if not isinstance(base_index, int) or not 0 <= base_index < len(base):
                raise ValueError(f'point {index}: base index out of range')
            relation = group.add(relation, base[base_index])
        recovered = int(row['recovered_scalar'])
        interval = float(row['target_ms'])
        phases = sum(float(row[key]) for key in
                     ('target_query_ms', 'target_pdp_ms',
                      'target_relation_check_ms', 'target_descent_ms',
                      'target_recovery_check_ms'))
        if (row['fixture_index'] != index or
                group.point_from_json(row['target']) != point or
                relation != point or recovered != int(answer['scalar']) or
                group.multiply(group.G, recovered) != point or
                row['group_verified'] is not True or
                not math.isfinite(interval) or interval <= 0 or
                abs(interval-phases) > 1.0 or
                not isinstance(row['probes'], int) or row['probes'] <= 0):
            raise ValueError(f'point {index}: relation, scalar or phase differs')
        times.append(interval)
        probes.append(row['probes'])
    return {'schema':'n83.root11-midk-screen-analysis/v1',
            'status':'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY',
            'k':k,'actual_usable_points':len(base),
            'effective_columns':summary['orbit_columns'],
            'base_hash':summary['base_hash'],
            'online_geometric_mean_ms':math.exp(statistics.mean(math.log(x) for x in times)),
            'online_median_ms':statistics.median(times),
            'mean_probes':statistics.mean(probes),
            'cold_process_ms':summary['timing_ms']['process_total'],
            'root_state_allocated_bytes':summary['root_state_allocated_bytes'],
            'root_table_allocated_bytes':summary['root_table_allocated_bytes'],
            'per_target_online_ms':times,
            'per_target_probes':probes,
            'raw_receipt_sha256':digest(folder / 'run_receipt.json')}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('folder',type=Path)
    parser.add_argument('--k',type=int,required=True)
    args=parser.parse_args()
    import run
    plan=run.frozen(Path(json.loads((HERE / 'plan.json').read_text())['binary_path']))
    result=verify_one(args.folder.resolve(),args.k,plan)
    print(json.dumps({'k':args.k,'status':result['status'],
                      'online_geometric_mean_ms':result['online_geometric_mean_ms']},sort_keys=True))


if __name__=='__main__':
    main()
