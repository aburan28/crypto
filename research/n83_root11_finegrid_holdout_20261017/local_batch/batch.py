#!/usr/bin/env python3
"""Secondary shared-setup N83 query-distribution diagnostic on 128 points."""

import hashlib
import importlib.util
import json
import math
import statistics
import subprocess
import time
from pathlib import Path

ROOT = Path('/Volumes/SSD990/crypto/worktrees/n83-final-online-holdout-20261017')
OUT = Path(__file__).resolve().parent
BIN = Path('/Volumes/SSD990/crypto/artifacts/n83_k_intermediate_screen_20261017/target/release/examples/koblitz_orbit_dlp_fast')
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
DOMAIN = b'n83-root11-finegrid-batch-20261017-v1'
COUNT = 128
K_VALUES = (2752, 2816, 2880, 2944, 3072)


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


dynamic = module('n83_batch_dynamic',
                 ROOT / 'research/n83_root11_dynamic_20261016/run.py')
group = module('n83_batch_group',
               ROOT / 'research/n83_cold_factor_base_20261016/verify.py').group


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_once(path, blob):
    if path.exists() and path.read_bytes() != blob:
        raise ValueError(f'existing receipt differs: {path}')
    path.write_bytes(blob)


def pretty(value):
    return (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()


def fixture():
    points, scalars, target_lines, answer_lines = [], [], [], []
    for index in range(COUNT):
        seed = hashlib.sha256(DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (group.R - 1)
        point = group.multiply(group.G, scalar)
        if point is None or not group.on_curve(point):
            raise ValueError('invalid generated subgroup point')
        points.append(point)
        scalars.append(scalar)
        target_lines.append(json.dumps([str(x) for x in point],
                                       separators=(',', ':')) + '\n')
        answer_lines.append(json.dumps({'fixture_index': index,
                                        'scalar': str(scalar),
                                        'seed_sha256': seed.hex()},
                                       sort_keys=True, separators=(',', ':')) + '\n')
    if len(set(points)) != COUNT:
        raise ValueError('duplicate target point')
    target = OUT / 'batch_targets.jsonl'
    answers = OUT / 'batch_answers.jsonl'
    write_once(target, ''.join(target_lines).encode())
    write_once(answers, ''.join(answer_lines).encode())
    return target, answers, points, scalars


def verify(folder, points, scalars):
    summary = json.loads((folder / 'stdout.jsonl').read_text())
    targets = [json.loads(line) for line in
               (folder / 'targets_result.jsonl').read_text().splitlines()]
    if (len(targets) != COUNT or summary['targets'] != COUNT or
            summary['targets_solved'] != COUNT or summary['targets_failed'] != 0):
        raise ValueError('query count or recovery status differs')
    base = group.replay_base(
        json.loads((folder / 'base_dump.jsonl').read_text()), summary)
    if len(base) != summary['factor_base_points']:
        raise ValueError('independent factor-base replay differs')
    times, probes = [], []
    for index, target in enumerate(targets):
        recovered = int(target['recovered_scalar'])
        indices = target['point_indices']
        relation = None
        for base_index in indices:
            if not 0 <= base_index < len(base):
                raise ValueError('factor-base index out of range')
            relation = group.add(relation, base[base_index])
        if (target['fixture_index'] != index or
                group.point_from_json(target['target']) != points[index] or
                recovered != scalars[index] or relation != points[index] or
                target['group_verified'] is not True):
            raise ValueError(f'target {index}: scalar or relation replay differs')
        timing = float(target['target_ms'])
        if not math.isfinite(timing) or timing <= 0:
            raise ValueError('nonpositive target time')
        times.append(timing)
        probes.append(target['probes'])
    return {'status': 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY',
            'target_count': COUNT, 'actual_factor_base_points': len(base),
            'folded_columns': summary['orbit_columns'],
            'online_mean_ms': statistics.mean(times),
            'online_median_ms': statistics.median(times),
            'online_geometric_mean_ms': math.exp(statistics.mean(map(math.log, times))),
            'mean_target_probes': statistics.mean(probes),
            'median_target_probes': statistics.median(probes),
            'root_state_allocated_bytes': summary['root_state_allocated_bytes'],
            'root_table_allocated_bytes': summary['root_table_allocated_bytes'],
            'cold_process_ms': summary['timing_ms']['process_total'],
            'per_target_online_ms': times,
            'per_target_probes': probes,
            'base_hash': summary['base_hash']}


def main():
    if sha256(ROOT / 'examples/koblitz_orbit_dlp_fast.rs') != SOURCE_SHA256:
        raise ValueError('current source differs')
    if not BIN.is_file():
        raise ValueError('offline release binary missing')
    target, answers, points, scalars = fixture()
    write_once(OUT / 'batch_inputs.json', pretty({
        'schema': 'n83.intermediate-k-batch-screen/v1',
        'workload_role': 'secondary_shared_setup_query_distribution_diagnostic',
        'hardware_tier': 'ordinary_local_exploratory',
        'source_sha256': SOURCE_SHA256, 'binary_sha256': sha256(BIN),
        'targets_sha256': sha256(target), 'answers_sha256': sha256(answers),
        'point_count': COUNT, 'scalar_rule':
            '1+little_endian(SHA256(domain||u32_le(index))) mod (r-1)',
        'domain': DOMAIN.decode(), 'k_values': K_VALUES,
        'cap_seconds_per_arm': 1200}))
    for k in K_VALUES:
        folder = OUT / f'batch_k{k}'
        folder.mkdir(exist_ok=True)
        command = [str(BIN), f'construct:83:1:{k}', str(target), '29',
                   str(folder / 'targets_result.jsonl')]
        env = dynamic.environment(folder, 'root11_dynamic')
        write_once(folder / 'command.json', pretty({
            'argv': command, 'environment': env, 'cap_seconds': 1200}))
        receipt = folder / 'batch_result.json'
        if receipt.exists():
            print(f'K{k}: existing batch receipt; not rerunning', flush=True)
            continue
        print(f'K{k}: starting 128-point secondary batch', flush=True)
        start = time.monotonic()
        try:
            process = subprocess.run(command, cwd=ROOT, env=env,
                                     capture_output=True, text=True,
                                     timeout=1200, check=False)
            stdout, stderr = process.stdout, process.stderr
            result = {'producer_status': 'EXIT_0' if process.returncode == 0
                      else 'PROCESS_FAILED', 'exit_code': process.returncode}
        except subprocess.TimeoutExpired as error:
            stdout = (error.stdout or b'').decode(errors='replace')
            stderr = (error.stderr or b'').decode(errors='replace')
            result = {'producer_status': 'TIMEOUT', 'exit_code': None}
        result['outer_wall_seconds'] = time.monotonic() - start
        (folder / 'stdout.jsonl').write_text(stdout)
        (folder / 'stderr.txt').write_text(stderr)
        if result['producer_status'] == 'EXIT_0':
            try:
                result.update(verify(folder, points, scalars))
            except Exception as error:
                result.update({'status': 'VERIFICATION_FAILED',
                               'verification_error': str(error)})
        else:
            result['status'] = 'CENSORED_PRODUCER_STATUS'
        write_once(receipt, pretty({'k': k, **result}))
        print(f'K{k}: {result["status"]}, {result["outer_wall_seconds"]:.1f} s',
              flush=True)


if __name__ == '__main__':
    main()
