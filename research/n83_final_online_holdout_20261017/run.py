#!/usr/bin/env python3
"""Run one six-point N83 final online Tier 1 stratum."""

import argparse
import importlib.util
import json
import os
import shutil
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
TARGETS_SHA256 = '4d167f8b2b85bcf6c599494daf3b3a61210a9dc8c0d65ad3ebcdf07cd6d5c4c7'
ANSWERS_SHA256 = 'f94a313ea6e7c414a4fc85c3ec55e384f6c6be14febbe3304d4b9dde1e54e484'
PAIR_RHO_SEED_BASE = 2026102900
IC_CAP_SMALL_SECONDS = 480
IC_CAP_LARGE_SECONDS = 900
RHO_CAP_SECONDS = 120
GLOBAL_PRODUCER_CAP_SECONDS = 190 * 60


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


panel = module('n83_final_online_panel',
               ROOT / 'research/n83_smallk_challenge_20261017/run.py')
generator = module('n83_final_online_generator', HERE / 'generate.py')
choice_module = module('n83_final_online_choice', HERE / 'arms.py')
finegrid = panel.finegrid
holdout = panel.holdout
dynamic = panel.dynamic
cold = panel.cold
parent = panel.parent


def fixture():
    hashes = finegrid.frozen_inputs()
    targets, answers = generator.expected_rows()
    manifest = json.loads((HERE / 'fixture_manifest.json').read_text())
    lines = (HERE / 'targets.jsonl').read_text().splitlines(keepends=True)
    if (parent.digest(HERE / 'targets.jsonl') != TARGETS_SHA256 or
            parent.digest(HERE / 'answers.jsonl') != ANSWERS_SHA256 or
            manifest['schema'] != 'n83.final-online-holdout-fixture/v1' or
            manifest['targets_sha256'] != TARGETS_SHA256 or
            manifest['answers_sha256'] != ANSWERS_SHA256 or
            manifest['target_count'] != 12 or
            generator.canonical_lines(targets) != ''.join(lines) or
            generator.canonical_lines(answers) !=
            (HERE / 'answers.jsonl').read_text() or
            len(lines) != 12 or any(not line.endswith('\n') for line in lines)):
        raise ValueError('sealed final online fixture differs')
    return hashes, lines


def schedule(half, arms):
    if (half not in (0, 1) or len(arms) not in (3, 4) or
            len({arm[3] for arm in arms}) != len(arms)):
        raise ValueError('three or four distinct frozen arms and half 0 or 1 required')
    events = []
    for local in range(6):
        index = half * 6 + local
        ordered = arms[index % len(arms):] + arms[:index % len(arms)]
        for arm in ordered:
            name = arm[0]
            fixed = arms.index(arm)
            events.append({'point_index': index, 'candidate': name,
                           'name': f'cal_before_{name}', 'kind': 'calibration'})
            pair = ('ic', 'rho') if (index + fixed) % 2 == 0 else ('rho', 'ic')
            for kind in pair:
                events.append({'point_index': index, 'candidate': name,
                               'name': name if kind == 'ic' else f'rho_{name}',
                               'kind': kind})
            events.append({'point_index': index, 'candidate': name,
                           'name': f'cal_after_{name}', 'kind': 'calibration'})
    return events


def run(binary, rho_binary, output, half, paths):
    hashes, lines = fixture()
    arms, choice = choice_module.choose(*paths)
    if output.exists():
        raise ValueError('output directory already exists')
    output.mkdir(parents=True)
    events = schedule(half, arms)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host['schema'] = 'n83.final-online-holdout-host/v1'
    host['parent_instance_id'] = host['instance_id']
    host['instance_id'] = f"{host['instance_id']}-final-online-half{half}"
    host['stratum_half'] = half
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
        'schema': 'n83.final-online-holdout-source/v1',
        'hashes': hashes, 'source_sha256': SOURCE_SHA256,
        'binary_sha256': parent.digest(ic_copy),
        'rho_binary_sha256': parent.digest(rho_copy),
        'runner_sha256': parent.digest(Path(__file__)),
        'targets_sha256': TARGETS_SHA256,
        'answers_sha256': ANSWERS_SHA256,
        'fixture_manifest_sha256': parent.digest(HERE / 'fixture_manifest.json'),
        'half': half, 'arms': arms, 'choice': choice, 'events': events,
        'rho_seed_rule': '2026102900+arm_count*point_index+fixed_candidate_index',
        'calibration_seed': int(cold.RHO_SEED),
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
            environment = dynamic.environment(folder, layout)
            cap = IC_CAP_LARGE_SECONDS if k >= 2000 else IC_CAP_SMALL_SECONDS
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
    parser.add_argument('--fixture-preflight', action='store_true')
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--half', type=int, choices=(0, 1), required=True)
    for name in ('five-validated', 'five-selection',
                 'k992-validated', 'k992-selection'):
        parser.add_argument(f'--{name}', type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--rho-binary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes, _ = fixture()
    if args.fixture_preflight:
        print(json.dumps({'hashes': hashes, 'sealed_target_count': 12},
                         sort_keys=True))
        return 0
    paths = (args.five_validated, args.five_selection,
             args.k992_validated, args.k992_selection)
    if any(item is None for item in paths):
        parser.error('all four panel evidence paths are required')
    paths = tuple(path.resolve() for path in paths)
    arms, choice = choice_module.choose(*paths)
    if args.preflight:
        print(json.dumps({'hashes': hashes, 'choice': choice,
                          'arms': arms, 'events': schedule(args.half, arms)},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error('--binary, --rho-binary and --output are required')
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.half, paths)


if __name__ == '__main__':
    raise SystemExit(main())
