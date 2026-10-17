#!/usr/bin/env python3
"""Independent full scalar and four-point relation replay for the batch screen."""

import hashlib
import importlib.util
import json
import math
import statistics
from pathlib import Path

ROOT = Path('/Volumes/SSD990/crypto/worktrees/n83-final-online-holdout-20261017')
HERE = Path(__file__).resolve().parent
COUNT = 128
K_VALUES = (2752, 2816, 2880, 2944, 3072)
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY = Path('/Volumes/SSD990/crypto/artifacts/n83_k_intermediate_screen_20261017/target/release/examples/koblitz_orbit_dlp_fast')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


group = module('n83_batch_replay_group',
               ROOT / 'research/n83_indexed5_20261013/analyze.py')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def manifest(k):
    return read(HERE / f'candidate_k{k}.json')


def fixture():
    inputs = read(HERE / 'batch_inputs.json')
    plan = read(HERE / 'batch_plan.json')
    if (inputs['source_sha256'] != SOURCE_SHA256 or
            sha256(HERE / 'batch_targets.jsonl') != inputs['targets_sha256'] or
            sha256(HERE / 'batch_answers.jsonl') != inputs['answers_sha256'] or
            sha256(BINARY) !=
            inputs['binary_sha256'] or
            sha256(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            plan['source_sha256'] != SOURCE_SHA256 or
            plan['binary_sha256'] != inputs['binary_sha256'] or
            plan['targets_sha256'] != inputs['targets_sha256'] or
            plan['answers_sha256'] != inputs['answers_sha256'] or
            plan['batch_script_sha256'] != sha256(HERE / 'batch.py') or
            plan['k_values'] != list(K_VALUES)):
        raise ValueError('source, binary or fixture hash differs')
    targets = [group.point_from_json(json.loads(line)) for line in
               (HERE / 'batch_targets.jsonl').read_text().splitlines()]
    answers = [json.loads(line) for line in
               (HERE / 'batch_answers.jsonl').read_text().splitlines()]
    if len(targets) != COUNT or len(answers) != COUNT or len(set(targets)) != COUNT:
        raise ValueError('fixture coverage or uniqueness differs')
    for index, (point, answer) in enumerate(zip(targets, answers)):
        scalar = int(answer['scalar'])
        seed = hashlib.sha256(
            b'n83-root11-finegrid-batch-20261017-v1' +
            index.to_bytes(4, 'little')).digest()
        if (answer['fixture_index'] != index or
                answer['seed_sha256'] != seed.hex() or
                scalar != 1 + int.from_bytes(seed, 'little') % (group.R - 1) or
                group.multiply(group.G, scalar) != point):
            raise ValueError(f'fixture {index} scalar replay differs')
    return inputs, targets, answers


def replay(k, inputs, targets, answers):
    folder = HERE / f'batch_k{k}'
    original = read(folder / 'batch_result.json')
    summary = read(folder / 'stdout.jsonl')
    candidate = manifest(k)
    command = read(folder / 'command.json')
    binary = BINARY
    if (original['producer_status'] != 'EXIT_0' or original['exit_code'] != 0 or
            (folder / 'stderr.txt').read_text() or
            command['argv'] != [str(binary), f'construct:83:1:{k}',
                                str(HERE / 'batch_targets.jsonl'), '29',
                                str(folder / 'targets_result.jsonl')] or
            summary['orbit_columns'] != k or summary['rank'] != k or
            summary['factor_base_points'] != 166 * k or
            summary['targets'] != COUNT or summary['targets_solved'] != COUNT or
            summary['targets_failed'] != 0 or
            candidate['factor_base']['actual_usable_points'] != 166 * k or
            candidate['factor_base']['effective_columns'] != k or
            candidate['factor_base']['source_base_hash_blake3'] !=
            summary['base_hash']):
        raise ValueError(f'K{k}: producer, command, rank or candidate differs')
    base = group.replay_base(read(folder / 'base_dump.jsonl'), summary)
    records = [json.loads(line) for line in
               (folder / 'targets_result.jsonl').read_text().splitlines()]
    if len(base) != 166 * k or len(records) != COUNT:
        raise ValueError(f'K{k}: base or result coverage differs')
    times, probes = [], []
    for index, (record, point, answer) in enumerate(zip(records, targets, answers)):
        recovered = int(record['recovered_scalar'])
        indices = record['point_indices']
        if not isinstance(indices, list) or len(indices) != 4:
            raise ValueError(f'K{k} point {index}: incomplete relation')
        relation = None
        for base_index in indices:
            if not isinstance(base_index, int) or not 0 <= base_index < len(base):
                raise ValueError(f'K{k} point {index}: invalid base index')
            relation = group.add(relation, base[base_index])
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
            raise ValueError(f'K{k} point {index}: scalar, relation or phase differs')
        times.append(interval)
        probes.append(record['probes'])
    if (times != original['per_target_online_ms'] or
            probes != original['per_target_probes'] or
            not math.isclose(statistics.mean(times), original['online_mean_ms'])):
        raise ValueError(f'K{k}: original result rows differ')
    return {'status': 'PASS_FULL_INDEPENDENT_REPLAY', 'k': k,
            'candidate_id': candidate['candidate_id'],
            'target_count': COUNT,
            'summary_sha256': sha256(folder / 'stdout.jsonl'),
            'base_dump_sha256': sha256(folder / 'base_dump.jsonl'),
            'target_results_sha256': sha256(folder / 'targets_result.jsonl'),
            'batch_result_sha256': sha256(folder / 'batch_result.json'),
            'mean_online_ms': statistics.mean(times),
            'mean_probes': statistics.mean(probes)}


def main():
    inputs, targets, answers = fixture()
    rows = [replay(k, inputs, targets, answers) for k in K_VALUES]
    result = {'schema': 'n83.intermediate-k-batch-replay/v1',
              'status': 'PASS_FULL_INDEPENDENT_REPLAY',
              'targets_sha256': inputs['targets_sha256'], 'rows': rows}
    path = HERE / 'batch_full_replay.json'
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError('sealed replay differs')
    path.write_text(blob)
    print(json.dumps({'status': result['status'],
                      'verified_candidate_count': len(rows)}))


if __name__ == '__main__':
    main()
