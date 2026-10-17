#!/usr/bin/env python3
"""Seal new N83 public points before high-K candidate selection."""
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
DOMAIN = b'n83-root11-highk-holdout24-20261017-v1'
COUNT = 24
EXCLUDED = (
    'n83_final_online_holdout_20261017/targets.jsonl',
    'n83_root11_finegrid_holdout_20261017/targets.jsonl',
    'n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl',
    'n83_intermediate_factor_base_holdout_20261017/targets.jsonl',
    'n83_smallk_challenge_20261017/targets.jsonl',
    'n83_k992_cold_followup_20261017/targets.jsonl',
    'n83_k3072_current_source_20261017/targets.jsonl',
    'n83_highk_tiebreak_20261017/targets.jsonl',
    'n83_highk_k2543_tiebreak_20261017/targets.jsonl',
    'n83_e2e_final_holdout_20261016/targets.jsonl',
    'n83_e2e_precision_20261016/targets.jsonl',
)
spec = importlib.util.spec_from_file_location(
    'n83_highk_holdout_group', ROOT/'research/n83_indexed5_20261013/analyze.py')
group = importlib.util.module_from_spec(spec)
spec.loader.exec_module(group)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_once(path, value):
    if path.exists() and path.read_bytes() != value:
        raise ValueError(f'sealed fixture differs: {path}')
    path.write_bytes(value)


def main():
    excluded_points, excluded_hashes = set(), {}
    for name in EXCLUDED:
        path = ROOT/'research'/name
        excluded_hashes[name] = sha256(path)
        for line in path.read_text().splitlines():
            excluded_points.add(group.point_from_json(json.loads(line)))
    targets, answers, points = [], [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN+index.to_bytes(4,'little')).digest()
        scalar = 1+int.from_bytes(seed,'little')%(group.R-1)
        point = group.multiply(group.G,scalar)
        if (point is None or not group.on_curve(point) or
                group.multiply(point,group.R) is not None):
            raise ValueError(f'point {index} outside subgroup')
        points.append(point)
        targets.append(json.dumps([str(v) for v in point],separators=(',',':'))+'\n')
        answers.append(json.dumps({'fixture_index':index,'scalar':str(scalar),
                                   'seed_sha256':seed.hex()},
                                  separators=(',',':'),sort_keys=True)+'\n')
    if len(set(points))!=COUNT or set(points)&excluded_points:
        raise ValueError('duplicate or previously used holdout point')
    target_path,answer_path=HERE/'targets.jsonl',HERE/'answers.jsonl'
    write_once(target_path,''.join(targets).encode())
    write_once(answer_path,''.join(answers).encode())
    manifest={'schema':'n83.root11-highk-holdout-fixture/v1',
              'curve_id':'EC1N83Ckb1a1haacb69ee9f37',
              'subgroup_order':group.R,'point_count':COUNT,
              'domain':DOMAIN.decode(),
              'scalar_rule':'1+little_endian(SHA256(domain||u32_le(index))) mod (r-1)',
              'targets_sha256':sha256(target_path),
              'answers_sha256':sha256(answer_path),
              'generator_sha256':sha256(Path(__file__)),
              'selector_sha256':sha256(HERE/'select.py'),
              'excluded_files_sha256':excluded_hashes,
              'status':'SEALED_UNMEASURED'}
    write_once(HERE/'manifest.json',
               (json.dumps(manifest,sort_keys=True,indent=2)+'\n').encode())
    print(json.dumps({'status':manifest['status'],
                      'targets_sha256':manifest['targets_sha256'],
                      'answers_sha256':manifest['answers_sha256'],
                      'nonoverlap_file_count':len(EXCLUDED)},sort_keys=True))


if __name__=='__main__':
    main()
