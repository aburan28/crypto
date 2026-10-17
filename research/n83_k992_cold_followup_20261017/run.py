#!/usr/bin/env python3
"""Run one five-point, three-arm N83 K992 Tier 1 stratum."""

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
TARGETS_SHA256 = '89c7208dd4fc7a5621cf727b3fac9a6714f04f66a312e516b418c31764e7547e'
ANSWERS_SHA256 = '363adc7f9bded3f4a25a80c7b27da5a2305afb30ee6741932cbe76f455c2865b'
K992 = ('k992root11', 992, 'root11_dynamic',
        'IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094')
K1024 = ('k1024root11', 1024, 'root11_dynamic',
         'IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289')
K1107 = ('k1107root11', 1107, 'root11_dynamic',
         'IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa')
ARMS = (K992, K1024, K1107)
PAIR_RHO_SEED_BASE = 2026102700
IC_CAP_SECONDS = 480
RHO_CAP_SECONDS = 120
GLOBAL_PRODUCER_CAP_SECONDS = 120 * 60


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


finalist = module('n83_k992_finalist_runner',
                  ROOT / 'research/n83_e2e_final_holdout_20261016/run.py')
generator = module('n83_k992_fixture_generator', HERE / 'generate.py')
finegrid = finalist.finegrid
holdout = finalist.holdout
dynamic = finalist.dynamic
cold = finalist.cold
parent = finalist.parent


def fixture():
    hashes = finegrid.frozen_inputs()
    targets, answers = generator.expected_rows()
    manifest = json.loads((HERE / 'fixture_manifest.json').read_text())
    lines = (HERE / 'targets.jsonl').read_text().splitlines(keepends=True)
    if (parent.digest(HERE / 'targets.jsonl') != TARGETS_SHA256 or
            parent.digest(HERE / 'answers.jsonl') != ANSWERS_SHA256 or
            manifest['schema'] != 'n83.k992-cold-followup-fixture/v1' or
            manifest['targets_sha256'] != TARGETS_SHA256 or
            manifest['answers_sha256'] != ANSWERS_SHA256 or
            manifest['target_count'] != 10 or
            generator.canonical_lines(targets) != ''.join(lines) or
            generator.canonical_lines(answers) !=
            (HERE / 'answers.jsonl').read_text() or
            len(lines) != 10 or any(not line.endswith('\n') for line in lines)):
        raise ValueError('sealed K992 followup fixture differs')
    return hashes, lines


def schedule(half, arms=ARMS):
    if half not in (0, 1) or arms != ARMS:
        raise ValueError('three frozen candidates and half 0 or 1 required')
    events = []
    for local in range(5):
        index = half * 5 + local
        ordered = arms[index % 3:] + arms[:index % 3]
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


def run(binary, rho_binary, output, half):
    hashes, lines = fixture()
    arms = ARMS
    if output.exists():
        raise ValueError('output directory already exists')
    output.mkdir(parents=True)
    events = schedule(half, arms)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host['schema'] = 'n83.k992-cold-followup-host/v1'
    host['parent_instance_id'] = host['instance_id']
    host['instance_id'] = f"{host['instance_id']}-k992-cold-half{half}"
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
        'schema': 'n83.k992-cold-followup-source/v1',
        'hashes': hashes, 'source_sha256': SOURCE_SHA256,
        'binary_sha256': parent.digest(ic_copy),
        'rho_binary_sha256': parent.digest(rho_copy),
        'runner_sha256': parent.digest(Path(__file__)),
        'targets_sha256': TARGETS_SHA256,
        'answers_sha256': ANSWERS_SHA256,
        'fixture_manifest_sha256': parent.digest(HERE / 'fixture_manifest.json'),
        'half': half, 'arms': arms, 'events': events,
        'rho_seed_rule': '2026102700+3*point_index+fixed_candidate_index',
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
                    PAIR_RHO_SEED_BASE + 3 * index +
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
    parser.add_argument('--fixture-preflight', action='store_true')
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--half', type=int, choices=(0, 1), required=True)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--rho-binary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes, _ = fixture()
    if args.fixture_preflight:
        print(json.dumps({'hashes': hashes, 'arms': ARMS,
                          'events': schedule(args.half)}, sort_keys=True))
        return 0
    if args.preflight:
        print(json.dumps({'hashes': hashes, 'arms': ARMS,
                          'events': schedule(args.half)},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error('--binary, --rho-binary and --output are required')
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.half)


if __name__ == '__main__':
    raise SystemExit(main())
