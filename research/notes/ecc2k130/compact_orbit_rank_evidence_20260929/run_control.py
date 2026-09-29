#!/usr/bin/env python3
"""Run the fixed n13 trace and untraced controls without a timing claim."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import traceback

from verify_rank import Curve, Field, point, read_jsonl, verify

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def deterministic(value):
    if isinstance(value, dict):
        return {key: deterministic(item) for key, item in value.items()
                if '_ms' not in key and key != 'peak_rss_bytes'}
    if isinstance(value, list):
        return [deterministic(item) for item in value]
    return value


def run_one(binary: Path, out: Path, traced: bool) -> None:
    out.mkdir()
    env = os.environ.copy()
    env['KIC_DUMP_BASE'] = str(out / 'base.jsonl')
    env.pop('KIC_DUMP_RANK', None)
    if traced:
        env['KIC_DUMP_RANK'] = str(out / 'rank.jsonl')
    command = [str(binary), 'construct:13:0:2', str(HERE / 'targets_n13.txt'),
               '7', str(out / 'targets.jsonl')]
    completed = subprocess.run(command, cwd=ROOT, env=env, capture_output=True,
                               text=True, timeout=60, check=False)
    (out / 'summary.jsonl').write_text(completed.stdout)
    (out / 'stderr.txt').write_text(completed.stderr)
    (out / 'command.json').write_text(json.dumps({
        'argv': command, 'exit_code': completed.returncode, 'traced': traced,
        'binary_sha256': sha(binary), 'source_head': subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'host': platform.platform(),
    }, sort_keys=True, indent=2) + '\n')
    assert completed.returncode == 0, (traced, completed.stderr)
    assert len(read_jsonl(out / 'summary.jsonl')) == 1
    assert len(read_jsonl(out / 'targets.jsonl')) == 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    assert not args.out.exists(), 'refusing to overwrite control evidence'
    args.out.mkdir(parents=True)
    receipt = {'status': 'STOP', 'schema': 'compact-orbit-rank-control-v1'}
    try:
        binary = args.binary.resolve()
        assert binary.is_file()
        run_one(binary, args.out / 'traced', True)
        run_one(binary, args.out / 'control', False)
        left = args.out / 'traced'
        right = args.out / 'control'
        assert deterministic(read_jsonl(left / 'summary.jsonl')) == deterministic(
            read_jsonl(right / 'summary.jsonl'))
        assert deterministic(read_jsonl(left / 'targets.jsonl')) == deterministic(
            read_jsonl(right / 'targets.jsonl'))
        assert (left / 'base.jsonl').read_bytes() == (right / 'base.jsonl').read_bytes()
        replay = verify(left / 'rank.jsonl', left / 'base.jsonl', left / 'summary.jsonl')
        (args.out / 'independent_rank_replay.json').write_text(
            json.dumps(replay, sort_keys=True, indent=2) + '\n')
        assert replay['status'] == 'PASS' and replay['rank'] == 2
        target = read_jsonl(left / 'targets.jsonl')[0]
        assert target['published_fixture_scalar'] == 7
        assert target['recovered_scalar'] == 7
        assert target['group_verified'] is True
        base = read_jsonl(left / 'base.jsonl')[0]
        curve = Curve(Field(13, base['field_modulus_low_terms']), 0)
        gen = point(target['generator'])
        q = point(target['target'])
        assert curve.on_curve(q) and curve.mul(7, gen) == q
        chosen = [point(base['factor_base_point_coordinates'][index])
                  for index in target['point_indices']]
        relation_sum = None
        for member in chosen:
            assert curve.on_curve(member)
            relation_sum = curve.add(relation_sum, member)
        assert relation_sum == q
        receipt.update({'status': 'PASS', 'rank': replay['rank'],
                        'rank_attempts': replay['attempts'],
                        'rank_relations': replay['relations'],
                        'representative_logs_verified': replay['representative_logs_verified'],
                        'target_scalars_verified': 1,
                        'traced_and_control_deterministic_outputs_match': True})
    except BaseException as exc:
        receipt.update({'error_type': type(exc).__name__, 'error': str(exc),
                        'traceback': traceback.format_exc()})
        raise
    finally:
        receipt['files'] = {str(path.relative_to(args.out)): sha(path)
                            for path in sorted(args.out.rglob('*')) if path.is_file()}
        (args.out / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({key: value for key, value in receipt.items() if key != 'files'},
                     sort_keys=True))


if __name__ == '__main__':
    main()
