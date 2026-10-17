#!/usr/bin/env python3
"""Replay bound current-source K3450/wide3 and K3500/root4 archives."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_COMMIT = '416489e09c3d8885b6c12f262935e6ba29373460'
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
BINARY_SHA256 = 'd71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906'
CHECKS = {
    'wide3': {
        'folder': 'k3450root11wide3', 'k': 3450, 'points': 572700,
        'base_hash': 'ae75c86565f5fa12fb6f30a5539f6f4638bcf1ca8b347a5cefa20f5073a9f572',
        'root_entries': 987887478, 'slot_bytes': 5,
        'allocated_bytes': 5368709120, 'bloom_profile': 'wide3',
        'candidate_id': 'IC1N83Ckb1a1fb572700PDP4rootRCguidedLAgaussTDdirectISO0hbde4a785a36e',
        'archive_sha256': '1495f924ed125d4ea7aef000f2265b60f750126053ebb537b73c0f90735aff19',
        'receipt': 'receipt_K3450_wide3.json', 'replay': 'replay_K3450_wide3.json',
    },
    'root4': {
        'folder': 'k3500root4standard', 'k': 3500, 'points': 581000,
        'base_hash': '3387d0c40335c49256ad2048542455572980bef501c3735fd8a26f42ac035e67',
        'root_entries': 1016729255, 'slot_bytes': 4,
        'allocated_bytes': 8589934592, 'bloom_profile': 'standard',
        'candidate_id': 'IC1N83Ckb1a1fb581000PDP4rootRCguidedLAgaussTDdirectISO0hb3a616c4e22f',
        'archive_sha256': '1f278f491e20e0780be9d225a7ce62e4345e9777763a7f7bf73377ac100ac3f5',
        'receipt': 'receipt_K3500_root4.json', 'replay': 'replay_K3500_root4.json',
    },
}

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value

def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()

def candidate(kind, summary, dump_path):
    if kind == 'wide3':
        builder = module('n83_wide3_archive_builder',
                         ROOT / 'research/n83_best_e2e_results_20261015/materialize.py')
        builder.SOURCE_COMMIT = SOURCE_COMMIT
        builder.SOURCE_SHA256 = SOURCE_SHA256
        builder.panel.SOURCE_COMMIT = SOURCE_COMMIT
        builder.panel.SOURCE_SHA256 = SOURCE_SHA256
        builder.BASE_HASHES[3450] = CHECKS[kind]['base_hash']
        builder.ROOT_KEYS[3450] = CHECKS[kind]['root_entries']
        field, curve = builder.panel.curve_record()
        return builder.candidate(3450, 'root11', 'wide3', summary,
                                 dump_path, field, curve)
    builder = module('n83_root4_archive_builder',
                     ROOT / 'research/n83_compact4_crosslayout_results_20261016/materialize.py')
    builder.SOURCE_COMMIT = SOURCE_COMMIT
    builder.SOURCE_SHA256 = SOURCE_SHA256
    builder.common.SOURCE_COMMIT = SOURCE_COMMIT
    builder.common.SOURCE_SHA256 = SOURCE_SHA256
    builder.common.panel.SOURCE_COMMIT = SOURCE_COMMIT
    builder.common.panel.SOURCE_SHA256 = SOURCE_SHA256
    builder.common.BASE_HASHES[3500] = CHECKS[kind]['base_hash']
    builder.common.ROOT_KEYS[3500] = CHECKS[kind]['root_entries']
    field, curve = builder.common.panel.curve_record()
    return builder.candidate(3500, 'root4', 'standard', summary,
                             dump_path, field, curve)

def replay(kind, archive, extracted):
    config = CHECKS[kind]
    if digest(archive) != config['archive_sha256']:
        raise ValueError('raw archive digest differs')
    folder = extracted / 'runs' / config['folder']
    receipt = json.loads((extracted / 'receipt.json').read_text())
    archived_replay = json.loads((extracted / 'replay.json').read_text())
    archived_candidate = json.loads((extracted / 'candidate.json').read_text())
    versioned_receipt = json.loads((HERE / config['receipt']).read_text())
    versioned_replay = json.loads((HERE / config['replay']).read_text())
    versioned_candidate = json.loads((HERE / f"candidate_{config['candidate_id']}.json").read_text())
    if (receipt != versioned_receipt or archived_replay != versioned_replay or
            archived_candidate != versioned_candidate):
        raise ValueError('versioned receipt, replay or candidate differs')
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    target = json.loads((folder / 'targets_result.jsonl').read_text())
    cold = module('n83_highk_archive_cold',
                  ROOT / 'research/n83_cold_factor_base_20261016/run.py')
    group = module('n83_highk_archive_group',
                   ROOT / 'research/n83_root11_dynamic_20261016/verify.py').group
    encoded = json.loads(cold.TARGETS.read_text())
    answer = json.loads(cold.ANSWERS.read_text())
    point = group.point_from_json(encoded)
    scalar = int(answer['scalar'])
    flags = receipt['environment_flags']
    expected_flags = ({'KIC_ROOT_ONLY11': '1', 'KIC_INDEXED_ROOT4': '0',
                       'KIC_INDEXED_ROOT5': '1', 'KIC_INDEXED5_FIXED30': '1',
                       'KIC_INDEXED4_FIXED31': '0', 'KIC_ROOT_BLOOM': '1',
                       'KIC_ROOT_BLOOM_PROFILE': 'wide3'} if kind == 'wide3' else
                      {'KIC_ROOT_ONLY11': '1', 'KIC_INDEXED_ROOT4': '1',
                       'KIC_INDEXED_ROOT5': '0', 'KIC_INDEXED5_FIXED30': '0',
                       'KIC_INDEXED4_FIXED31': '1', 'KIC_ROOT_BLOOM': '1',
                       'KIC_ROOT_BLOOM_PROFILE': 'standard'})
    if (digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            receipt['status'] != 'PASS_PRODUCER' or receipt['exit_code'] != 0 or
            receipt['source_sha256'] != SOURCE_SHA256 or
            receipt['binary_sha256'] != BINARY_SHA256 or
            receipt['target_sha256'] != digest(cold.TARGETS) or
            receipt['stdout_sha256'] != digest(folder / 'stdout.jsonl') or
            receipt['stderr_sha256'] != digest(folder / 'stderr.txt') or
            (folder / 'stderr.txt').read_bytes() or
            any(flags.get(key) != value for key, value in expected_flags.items()) or
            summary['n'] != 83 or summary['a'] != 1 or
            summary['factor_base_points'] != config['points'] or
            summary['rank'] != config['k'] or
            summary['root_table_entries'] != config['root_entries'] or
            summary['root_table_slot_bytes'] != config['slot_bytes'] or
            summary['root_table_allocated_bytes'] != config['allocated_bytes'] or
            summary['root_bloom_profile'] != config['bloom_profile'] or
            summary['root_bloom_bytes'] != 1073741824 or
            summary['base_hash'] != config['base_hash'] or
            summary['targets_solved'] != 1 or summary['targets_failed'] != 0 or
            target['fixture_index'] != 0 or target['target'] != encoded or
            int(target['recovered_scalar']) != scalar or
            group.multiply(group.G, scalar) != point):
        raise ValueError('source, producer, base, target or scalar differs')
    base = group.replay_base(json.loads((folder / 'base_dump.jsonl').read_text()),
                             summary)
    indices = target['point_indices']
    if (len(base) != config['points'] or len(indices) != 4 or
            any(not isinstance(i, int) or not 0 <= i < len(base)
                for i in indices)):
        raise ValueError('factor base or witness differs')
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError('four-point relation differs')
    keys = ('target_query_ms', 'target_pdp_ms', 'target_relation_check_ms',
            'target_descent_ms', 'target_recovery_check_ms')
    online = float(target['target_ms'])
    phases = [float(target[key]) for key in keys]
    if (not math.isfinite(online) or online <= 0 or
            any(not math.isfinite(value) or value < 0 for value in phases) or
            abs(math.fsum(phases) - online) > 1):
        raise ValueError('online phase accounting differs')
    computed = candidate(kind, summary, folder / 'base_dump.jsonl')
    if computed != versioned_candidate or computed['candidate_id'] != config['candidate_id']:
        raise ValueError('exact candidate identity differs')
    if (archived_replay['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY' or
            archived_replay['source_commit'] != SOURCE_COMMIT or
            archived_replay['source_sha256'] != SOURCE_SHA256 or
            archived_replay['binary_sha256'] != BINARY_SHA256 or
            archived_replay['receipt_sha256'] != digest(extracted / 'receipt.json') or
            archived_replay['base_hash'] != config['base_hash'] or
            archived_replay['root_table_entries'] != config['root_entries'] or
            archived_replay['actual_usable_points'] != len(base) or
            archived_replay['one_target_online_ms'] != online):
        raise ValueError('archived arithmetic replay summary differs')
    return {'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
            'candidate_id': config['candidate_id'],
            'actual_usable_points': len(base), 'one_target_online_ms': online}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('kind', choices=CHECKS)
    parser.add_argument('archive', type=Path)
    parser.add_argument('extracted', type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.kind, args.archive.resolve(),
                            args.extracted.resolve()), sort_keys=True))

if __name__ == '__main__':
    main()
