#!/usr/bin/env python3
"""Replay archived K1248/K1280 bases, relations, scalars and candidate IDs."""
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
ARCHIVE_SHA256 = 'f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685'
CANDIDATES = {
    1248: 'IC1N83Ckb1a1fb207168PDP4rootRCguidedLAgaussTDdirectISO0hc9eda193920f',
    1280: 'IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac',
}


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finegrid = module('smallk_finegrid', ROOT / 'research/n83_e2e_finegrid_20261016/materialize.py')
verifier = module('smallk_group', ROOT / 'research/n83_root11_dynamic_20261016/verify.py')
group = verifier.group


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def replay_one(extracted, k, frozen, point, scalar, field, curve):
    folder = extracted / 'runs' / f'k{k}root11'
    receipt_path = extracted / f'receipt_K{k}.json'
    receipt = json.loads(receipt_path.read_text())
    if receipt != json.loads((HERE / f'receipt_K{k}.json').read_text()):
        raise ValueError(f'K{k}: producer receipt differs')
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    target = json.loads((folder / 'targets_result.jsonl').read_text())
    recorded_folder = Path(receipt['command'][4]).parent
    expected_flags = {key: value for key, value in
                      finegrid.runner.dynamic.environment(recorded_folder, 'root11_dynamic').items()
                      if key.startswith('KIC_')}
    if (digest(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256 or
            receipt['k'] != k or receipt['status'] != 'PASS_PRODUCER' or
            len(receipt['command']) != 5 or
            receipt['command'][1] != f'construct:83:1:{k}' or
            receipt['command'][3] != '29' or
            Path(receipt['command'][4]).name != 'targets_result.jsonl' or
            receipt['exit_code'] != 0 or receipt['source_sha256'] != SOURCE_SHA256 or
            receipt['binary_sha256'] != BINARY_SHA256 or
            receipt['target_sha256'] != digest(finegrid.runner.cold.TARGETS) or
            receipt['environment_flags'] != expected_flags or
            receipt['stdout_sha256'] != digest(folder / 'stdout.jsonl') or
            receipt['stderr_sha256'] != digest(folder / 'stderr.txt') or
            (folder / 'stderr.txt').read_bytes() or
            summary['n'] != 83 or summary['a'] != 1 or
            summary['rank'] != k or summary['orbit_columns'] != k or
            summary['factor_base_points'] != 166 * k or
            summary['targets_solved'] != 1 or summary['targets_failed'] != 0 or
            target['fixture_index'] != 0 or target['target'] != frozen or
            target['group_verified'] is not True or
            int(target['recovered_scalar']) != scalar or
            group.multiply(group.G, scalar) != point):
        raise ValueError(f'K{k}: source, producer, target or scalar differs')
    base = group.replay_base(json.loads((folder / 'base_dump.jsonl').read_text()), summary)
    indices = target['point_indices']
    if (len(base) != 166 * k or len(indices) != 4 or
            any(not isinstance(i, int) or not 0 <= i < len(base) for i in indices)):
        raise ValueError(f'K{k}: base or witness differs')
    relation = None
    for index in indices:
        relation = group.add(relation, base[index])
    if relation != point:
        raise ValueError(f'K{k}: four-point relation differs')
    keys = ('target_query_ms', 'target_pdp_ms', 'target_relation_check_ms',
            'target_descent_ms', 'target_recovery_check_ms')
    online = float(target['target_ms'])
    phases = [float(target[key]) for key in keys]
    if (not math.isfinite(online) or online <= 0 or
            any(not math.isfinite(value) or value < 0 for value in phases) or
            abs(math.fsum(phases) - online) > 1):
        raise ValueError(f'K{k}: online phase accounting differs')
    candidate = finegrid.dynamic.candidate(k, 'root11_dynamic', summary,
                                            folder / 'base_dump.jsonl', field, curve)
    expected = json.loads((HERE / f'candidate_{CANDIDATES[k]}.json').read_text())
    archived = json.loads((extracted / f'candidate_K{k}.json').read_text())
    if candidate != expected or candidate != archived or candidate['candidate_id'] != CANDIDATES[k]:
        raise ValueError(f'K{k}: candidate identity differs')
    replayed = json.loads((extracted / f'replay_K{k}.json').read_text())
    if replayed != json.loads((HERE / f'replay_K{k}.json').read_text()):
        raise ValueError(f'K{k}: replay receipt differs')
    if (replayed['status'] != 'PASS_INDEPENDENT_ARITHMETIC_REPLAY' or
            replayed['source_commit'] != SOURCE_COMMIT or
            replayed['source_sha256'] != SOURCE_SHA256 or
            replayed['binary_sha256'] != BINARY_SHA256 or
            replayed['candidate_id'] != CANDIDATES[k] or
            replayed['receipt_sha256'] != digest(receipt_path) or
            replayed['base_hash'] != summary['base_hash'] or
            replayed['actual_usable_points'] != len(base) or
            replayed['one_target_online_ms'] != online):
        raise ValueError(f'K{k}: replay summary differs')
    return {'k': k, 'status': 'PASS_INDEPENDENT_ARITHMETIC_REPLAY',
            'candidate_id': CANDIDATES[k], 'actual_usable_points': len(base),
            'one_target_online_ms': online}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('archive', type=Path)
    parser.add_argument('extracted', type=Path)
    args = parser.parse_args()
    if digest(args.archive) != ARCHIVE_SHA256:
        raise ValueError('raw archive digest differs')
    extracted = args.extracted.resolve()
    frozen = json.loads(finegrid.runner.cold.TARGETS.read_text())
    answer = json.loads(finegrid.runner.cold.ANSWERS.read_text())
    point = group.point_from_json(frozen)
    scalar = int(answer['scalar'])
    finegrid.panel.SOURCE_COMMIT = SOURCE_COMMIT
    finegrid.panel.SOURCE_SHA256 = SOURCE_SHA256
    field, curve = finegrid.panel.curve_record()
    for k in (1248, 1280):
        print(json.dumps(replay_one(extracted, k, frozen, point, scalar, field, curve),
                         sort_keys=True))


if __name__ == '__main__':
    main()
