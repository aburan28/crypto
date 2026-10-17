#!/usr/bin/env python3
"""Run one six-point N83 target cancellation Tier 1 stratum."""

import argparse
import importlib.util
import json
import os
import shutil
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8109fb33374ce2aa66c6a5a0b6e144d3b0143168ac756f9e0bd70af60a911829'
TARGETS_SHA256 = 'aa02c031dedd97a5f4bb214838e550e4597f5af44b556c3f24f7e5f3eefd99f5'
ANSWERS_SHA256 = 'c0134002e32ff898c8799f2024649ef286d0fbf881dc5a1a3eff5b589e4ce7d2'
PAIR_RHO_SEED_BASE = 2026110400
IC_CAP_SECONDS = 900
RHO_CAP_SECONDS = 120
GLOBAL_PRODUCER_CAP_SECONDS = 300 * 60


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


panel = module('n83_finegrid_panel',
               ROOT / 'research/n83_smallk_challenge_20261017/run.py')
generator = module('n83_finegrid_generator', HERE / 'generate.py')
ARMS = (('parallel_join', 992, 'root11_dynamic', 'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h33c9c3b2f933'),
        ('parallel_cancel', 992, 'root11_dynamic', 'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h58ec337f71f7'))

def choose():
    return ARMS, {'schema':'n83.target-cancel-tier1-choice/v1',
                  'policy':'two_worker_join_vs_cooperative_cancel',
                  'ordered_candidate_ids':[arm[3] for arm in ARMS]}
finegrid = panel.finegrid
finegrid.screen.SOURCE_SHA256 = SOURCE_SHA256
finegrid.screen.dynamic.SOURCE_SHA256 = SOURCE_SHA256
finegrid.screen.parent.EXPECTED[finegrid.screen.SOURCE] = SOURCE_SHA256
holdout = panel.holdout
dynamic = panel.dynamic
cold = panel.cold
parent = panel.parent


def ic_environment(folder, layout, name):
    flags = dynamic.environment(folder, layout)
    flags['KIC_TARGET_THREADS'] = '2'
    flags['KIC_INDEXED5_ID30'] = '0'
    flags['KIC_TARGET_CANCEL'] = '1' if name == 'parallel_cancel' else '0'
    return flags


def fixture():
    hashes = finegrid.frozen_inputs()
    manifest = json.loads((HERE / 'manifest.json').read_text())
    lines = (HERE / 'targets.jsonl').read_text().splitlines(keepends=True)
    answers = [json.loads(line) for line in
               (HERE / 'answers.jsonl').read_text().splitlines()]
    points = [generator.group.point_from_json(json.loads(line)) for line in lines]
    if (parent.digest(HERE / 'targets.jsonl') != TARGETS_SHA256 or
            parent.digest(HERE / 'answers.jsonl') != ANSWERS_SHA256 or
            manifest['schema'] != 'n83.target-cancel-tier1-fixture/v1' or
            manifest['targets_sha256'] != TARGETS_SHA256 or
            manifest['answers_sha256'] != ANSWERS_SHA256 or
            manifest['point_count'] != 24 or
            manifest['generator_sha256'] != parent.digest(HERE / 'generate.py') or
            len(lines) != 24 or len(answers) != 24 or
            len(set(points)) != 24 or any(not line.endswith('\n') for line in lines)):
        raise ValueError('sealed cancellation fixture differs')
    for index, (point, answer) in enumerate(zip(points, answers)):
        seed = generator.hashlib.sha256(
            generator.DOMAIN + index.to_bytes(4, 'little')).digest()
        scalar = 1 + int.from_bytes(seed, 'little') % (generator.group.R - 1)
        if (answer['fixture_index'] != index or
                answer['seed_sha256'] != seed.hex() or
                int(answer['scalar']) != scalar or
                generator.group.multiply(generator.group.G, scalar) != point):
            raise ValueError(f'fixture point {index} differs')
    return hashes, lines


def schedule(stratum, arms):
    if (stratum not in range(4) or len(arms) != 2 or
            len({arm[3] for arm in arms}) != len(arms) or
            len({arm[0] for arm in arms}) != len(arms)):
        raise ValueError('two distinct frozen arms and stratum 0..3 required')
    events = []
    for local in range(6):
        index = stratum * 6 + local
        ordered = arms[index % len(arms):] + arms[:index % len(arms)]
        for arm in ordered:
            name = arm[0]
            fixed = arms.index(arm)
            for repetition in range(3):
                events.append({'point_index': index, 'candidate': name,
                               'name': f'cal_before_{repetition}_{name}',
                               'kind': 'calibration'})
            pair = ('ic', 'rho') if (index + fixed) % 2 == 0 else ('rho', 'ic')
            for kind in pair:
                events.append({'point_index': index, 'candidate': name,
                               'name': name if kind == 'ic' else f'rho_{name}',
                               'kind': kind})
            for repetition in range(3):
                events.append({'point_index': index, 'candidate': name,
                               'name': f'cal_after_{repetition}_{name}',
                               'kind': 'calibration'})
    return events


def run(binary, rho_binary, output, stratum):
    hashes, lines = fixture()
    arms, choice = choose()
    if output.exists():
        raise ValueError('output directory already exists')
    output.mkdir(parents=True)
    events = schedule(stratum, arms)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host['schema'] = 'n83.target-cancel-tier1-host/v1'
    host['parent_instance_id'] = host['instance_id']
    host['instance_id'] = f"{host['instance_id']}-cancel-stratum{stratum}"
    host['stratum'] = stratum
    host['topology'] = holdout.topology_receipt(host['visible_affinity'] or [])
    host['cpu_telemetry_before'] = holdout.cpu_telemetry()
    parent.write_json(output / 'host.json', host)
    memory_max = (host['cgroup_before']['memory.max'] or '').strip()
    memory_too_small = memory_max.isdigit() and int(memory_max) < cold.RSS_CAP_BYTES
    if (host['mem_total_kib'] is None or
            host['mem_total_kib'] * 1024 < cold.RSS_CAP_BYTES + 256 * 1024**2 or
            host['disk_free_bytes'] < 3 * 1024**3 or
            len(host['visible_affinity'] or []) < 4 or memory_too_small):
        parent.write_json(output / 'preflight.json',
                          {'status': 'REJECTED_HOST_CAPACITY', 'host': host})
        return 2
    parent.write_json(output / 'preflight.json',
                      {'status': 'PASS_CAPACITY_PREFLIGHT',
                       'rss_cap_bytes': cold.RSS_CAP_BYTES})
    executed = output / 'executed_bins'
    executed.mkdir()
    ic_copy = executed / 'koblitz_orbit_dlp_fast'
    rho_copy = executed / 'n83_rho_optimization_eval'
    shutil.copy2(binary, ic_copy)
    shutil.copy2(rho_binary, rho_copy)
    ic_copy.chmod(0o755)
    rho_copy.chmod(0o755)
    parent.write_json(output / 'source.json', {
        'schema': 'n83.target-cancel-tier1-source/v1',
        'hashes': hashes, 'source_sha256': SOURCE_SHA256,
        'binary_sha256': parent.digest(ic_copy),
        'rho_binary_sha256': parent.digest(rho_copy),
        'runner_sha256': parent.digest(Path(__file__)),
        'targets_sha256': TARGETS_SHA256,
        'answers_sha256': ANSWERS_SHA256,
        'fixture_manifest_sha256': parent.digest(HERE / 'manifest.json'),
        'stratum': stratum, 'arms': arms, 'choice': choice, 'events': events,
        'rho_seed_rule': '2026110400+2*point_index+fixed_candidate_index',
        'calibration_seed': int(cold.RHO_SEED),
        'calibration_repetitions_per_side': 3,
        'calibration_comparison': 'median_of_three_before_vs_after_10_percent',
        'calibration_target_sha256': parent.digest(cold.TARGETS),
        'calibration_answer_sha256': parent.digest(cold.ANSWERS),
        'global_producer_cap_seconds': GLOBAL_PRODUCER_CAP_SECONDS,
    })
    by_name = {arm[0]: arm for arm in arms}
    started = time.monotonic()
    for event in events:
        index, name, candidate, kind = (event[key] for key in
                                         ('point_index', 'name', 'candidate', 'kind'))
        point_output = output / 'points' / f'p{index:02d}'
        point_output.mkdir(parents=True, exist_ok=True)
        target = point_output / 'target.jsonl'
        if target.exists():
            if target.read_text() != lines[index]:
                raise ValueError('per-point target file differs')
        else:
            target.write_text(lines[index])
        folder = point_output / 'runs' / name
        if kind in ('rho', 'calibration'):
            rho_target = cold.TARGETS if kind == 'calibration' else target
            seed = (int(cold.RHO_SEED) if kind == 'calibration' else
                    PAIR_RHO_SEED_BASE + len(arms) * index +
                    next(i for i, arm in enumerate(arms) if arm[0] == candidate))
            command = [str(rho_copy), str(rho_target), str(seed),
                       str(cold.RHO_CONFIG), str(folder / 'receipt.json')]
            environment = {'PATH': os.environ.get('PATH', '/usr/bin:/bin')}
            cap = RHO_CAP_SECONDS
        else:
            _, k, layout, _ = by_name[candidate]
            command = [str(ic_copy), f'construct:83:1:{k}', str(target), '29',
                       str(folder / 'targets_result.jsonl')]
            environment = ic_environment(folder, layout, candidate)
            cap = IC_CAP_SECONDS
        elapsed = time.monotonic() - started
        if elapsed + cap + 30 > GLOBAL_PRODUCER_CAP_SECONDS:
            parent.write_json(output / 'producer_stop.json', {
                'status': 'GLOBAL_CAP_BEFORE_NEXT_ARM',
                'elapsed_seconds': elapsed, 'next_event': event,
                'global_producer_cap_seconds': GLOBAL_PRODUCER_CAP_SECONDS})
            break
        print(f'starting p{index:02d}/{name} at {parent.now()}', flush=True)
        usage = holdout.run_with_telemetry(folder, command, environment,
                                          cap, group)
        print(f"p{index:02d}/{name}: {usage['status']} after "
              f"{usage['wall_seconds']:.3f}s", flush=True)
    parent.write_json(output / 'host_telemetry_after.json',
                      holdout.cpu_telemetry())
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--stratum', type=int, choices=range(4), required=True)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--rho-binary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes, _ = fixture()
    arms, choice = choose()
    if args.preflight:
        print(json.dumps({'hashes': hashes, 'choice': choice,
                          'arms': arms, 'events': schedule(args.stratum, arms)},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error('--binary, --rho-binary and --output are required')
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.stratum)


if __name__ == '__main__':
    raise SystemExit(main())
