#!/usr/bin/env python3
"""Replay the archived K3072 base, relation, scalar, and candidate identity."""

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
ARCHIVE_SHA256 = 'ba747b054d24af88abfb6df66acbd2efb407b5c92fb96092b6f2f0b4aa7f65d1'
CANDIDATE_ID = 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module('k3072_finegrid_manifest',
                  ROOT / 'research/n83_e2e_finegrid_20261016/materialize.py')
verifier = module('k3072_independent_group',
                  ROOT / 'research/n83_root11_dynamic_20261016/verify.py')
group = verifier.group


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def replay(archive, extracted):
    if digest(archive) != ARCHIVE_SHA256:
        raise ValueError('raw archive digest differs')
    folder = extracted / 'runs/k3072root11'
    receipt = json.loads((extracted / 'receipt.json').read_text())
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    target = json.loads((folder / 'targets_result.jsonl').read_text())
    frozen_target = json.loads(finegrid.runner.cold.TARGETS.read_text())
    answer = json.loads(finegrid.runner.cold.ANSWERS.read_text())
    point = group.point_from_json(frozen_target)
    scalar = int(answer['scalar'])
    if (digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            receipt['status'] != 'PASS_PRODUCER' or receipt['exit_code'] != 0 or
            receipt['source_sha256'] != SOURCE_SHA256 or
            receipt['target_sha256'] != digest(finegrid.runner.cold.TARGETS) or
            receipt['stdout_sha256'] != digest(folder / 'stdout.jsonl') or
            receipt['stderr_sha256'] != digest(folder / 'stderr.txt') or
            (folder / 'stderr.txt').read_bytes() or
            summary['rank'] != 3072 or summary['factor_base_points'] != 509952 or
            summary['targets_solved'] != 1 or summary['targets_failed'] != 0 or
            target['fixture_index'] != 0 or target['target'] != frozen_target or
            int(target['recovered_scalar']) != scalar or
            group.multiply(group.G, scalar) != point):
        raise ValueError('source, producer, target or scalar differs')
    base = group.replay_base(json.loads((folder / 'base_dump.jsonl').read_text()),
                             summary)
    indices = target['point_indices']
    if (len(base) != 509952 or len(indices) != 4 or
            any(not isinstance(i, int) or not 0 <= i < len(base) for i in indices)):
        raise ValueError('base or witness differs')
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError('four-point relation differs')
    keys = ('target_query_ms', 'target_pdp_ms', 'target_relation_check_ms',
            'target_descent_ms', 'target_recovery_check_ms')
    online_ms = float(target['target_ms'])
    phases = [float(target[key]) for key in keys]
    if (not math.isfinite(online_ms) or online_ms <= 0 or
            any(not math.isfinite(x) or x < 0 for x in phases) or
            abs(math.fsum(phases) - online_ms) > 1):
        raise ValueError('online phases differ')
    finegrid.panel.SOURCE_COMMIT = SOURCE_COMMIT
    finegrid.panel.SOURCE_SHA256 = SOURCE_SHA256
    field, curve = finegrid.panel.curve_record()
    candidate = finegrid.dynamic.candidate(
        3072, 'root11_dynamic', summary, folder / 'base_dump.jsonl',
        field, curve)
    versioned = json.loads((HERE / f'candidate_{CANDIDATE_ID}.json').read_text())
    if candidate != versioned or candidate['candidate_id'] != CANDIDATE_ID:
        raise ValueError('exact candidate identity differs')
    replayed = json.loads((extracted / 'replay.json').read_text())
    if (replayed['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY' or
            replayed['candidate_id'] != CANDIDATE_ID or
            replayed['base_hash'] != summary['base_hash'] or
            replayed['actual_usable_points'] != len(base) or
            replayed['one_target_online_ms'] != online_ms):
        raise ValueError('archived replay summary differs')
    return {'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
            'candidate_id': CANDIDATE_ID, 'actual_usable_points': len(base),
            'one_target_online_ms': online_ms}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('archive', type=Path)
    parser.add_argument('extracted', type=Path,
                        help='directory made by tar -xzf archive')
    args = parser.parse_args()
    print(json.dumps(replay(args.archive.resolve(), args.extracted.resolve()),
                     sort_keys=True))


if __name__ == '__main__':
    main()
