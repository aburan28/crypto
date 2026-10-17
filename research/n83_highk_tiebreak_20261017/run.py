#!/usr/bin/env python3
"""Run one five-point N83 high-K Tier 1 tie-break stratum."""

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
TARGETS_SHA256 = 'c31ed8c9b8ceec71cf4570919bf01aa9bb8168077e37851a9ec144a55058a4b2'
ANSWERS_SHA256 = '5c104cfd26be3077b3b0081e46547d2b4c7a040809bd61b2e936198774422fd6'
WIDE3 = ('k3450wide3', 3450, 'root11_wide3',
         'IC1N83Ckb1a1fb572700PDP4rootRCguidedLAgaussTDdirectISO0hbde4a785a36e')
ROOT4 = ('k3500root4', 3500, 'root4_standard',
         'IC1N83Ckb1a1fb581000PDP4rootRCguidedLAgaussTDdirectISO0hb3a616c4e22f')
PAIR_RHO_SEED_BASE = 2026102100
RHO_CAP_SECONDS = 120
GLOBAL_PRODUCER_CAP_SECONDS = 240 * 60


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


prior = module('n83_highk_prior_runner',
               ROOT / 'research/n83_k3072_current_source_20261017/run.py')
prior_selector = module('n83_highk_prior_selector',
                        ROOT / 'research/n83_k3072_current_source_20261017/select.py')
generator = module('n83_highk_generator', HERE / 'generate.py')
finalist = prior.finalist
holdout = prior.holdout
cold = prior.cold
parent = prior.parent


def fixture():
    hashes, _ = prior.fixture()
    targets, answers = generator.expected_rows()
    manifest = json.loads((HERE / 'fixture_manifest.json').read_text())
    lines = (HERE / 'targets.jsonl').read_text().splitlines(keepends=True)
    if (parent.digest(HERE / 'targets.jsonl') != TARGETS_SHA256 or
            parent.digest(HERE / 'answers.jsonl') != ANSWERS_SHA256 or
            manifest['schema'] != 'n83.highk-tiebreak-fixture/v1' or
            manifest['targets_sha256'] != TARGETS_SHA256 or
            manifest['answers_sha256'] != ANSWERS_SHA256 or
            manifest['target_count'] != 10 or
            generator.canonical_lines(targets) != ''.join(lines) or
            generator.canonical_lines(answers) !=
            (HERE / 'answers.jsonl').read_text() or
            len(lines) != 10 or any(not line.endswith('\n') for line in lines)):
        raise ValueError('sealed high-K fixture differs')
    hashes[str((HERE / 'targets.jsonl').relative_to(ROOT))] = TARGETS_SHA256
    hashes[str((HERE / 'answers.jsonl').relative_to(ROOT))] = ANSWERS_SHA256
    return hashes, lines


def choose_anchor(validated_path, selection_path):
    validated = json.loads(validated_path.read_text())
    selection = json.loads(selection_path.read_text())
    if (selection != json.loads(json.dumps(prior_selector.select(validated))) or
            validated.get('schema') != 'n83.k3072-current-source-validated/v2' or
            validated.get('source_sha256') != SOURCE_SHA256 or
            validated.get('targets_sha256') != prior.TARGETS_SHA256 or
            validated.get('candidate_ids', [])[:2] !=
            [prior.K3072[3], prior.K3450[3]]):
        raise ValueError('predecessor validated evidence or selection differs')
    online = selection['online']
    selected_id = (online['selected_candidate_id'] if
                   selection.get('admission_status') == 'PASS' and
                   online.get('status') == 'SELECTED' else None)
    anchor_id = selected_id or prior.K3072[3]
    options = [prior.K3072, prior.K3450, *finalist.ARMS]
    matches = [arm for arm in options if arm[3] == anchor_id]
    if len(matches) != 1 or anchor_id not in validated['candidate_ids']:
        raise ValueError('predecessor anchor candidate is unknown')
    return (WIDE3, ROOT4, matches[0]), {
        'predecessor_admission_status': selection.get('admission_status'),
        'online_selection_status': online.get('status'),
        'selected_candidate_id': selected_id,
        'anchor_candidate_id': anchor_id,
        'fallback_used': selected_id is None,
        'validated_sha256': parent.digest(validated_path),
        'selection_sha256': parent.digest(selection_path),
    }


def schedule(half, arms):
    if half not in (0, 1) or len(arms) != 3 or arms[:2] != (WIDE3, ROOT4):
        raise ValueError('three frozen candidates and half 0 or 1 required')
    events = []
    for local in range(5):
        index = half * 5 + local
        ordered = arms[index % 3:] + arms[:index % 3]
        for arm in ordered:
            name, fixed = arm[0], arms.index(arm)
            events.append({'point_index': index, 'candidate': name,
                           'name': f'cal_before_{name}', 'kind': 'calibration'})
            for kind in (('ic', 'rho') if (index + fixed) % 2 == 0 else
                         ('rho', 'ic')):
                events.append({'point_index': index, 'candidate': name,
                               'name': name if kind == 'ic' else f'rho_{name}',
                               'kind': kind})
            events.append({'point_index': index, 'candidate': name,
                           'name': f'cal_after_{name}', 'kind': 'calibration'})
    return events


def environment(folder, arm):
    if arm in (WIDE3, ROOT4):
        receipt_name = ('receipt_K3450_wide3.json' if arm == WIDE3 else
                        'receipt_K3500_root4.json')
        receipt = json.loads((ROOT / 'research/n83_k3072_current_source_20261017' /
                              receipt_name).read_text())
        flags = dict(receipt['environment_flags'])
        flags['KIC_DUMP_BASE'] = str(folder / 'base_dump.jsonl')
        return {'PATH': os.environ.get('PATH', '/usr/bin:/bin'), **flags}
    return prior.dynamic.environment(folder, arm[2])


def cap(arm):
    if arm in (WIDE3, ROOT4, prior.K3450):
        return 1200
    if arm == prior.K3072:
        return 900
    return 480


def run(binary, rho_binary, output, half, validated_path, selection_path):
    hashes, lines = fixture()
    arms, choice = choose_anchor(validated_path, selection_path)
    if output.exists():
        raise ValueError('output directory already exists')
    output.mkdir(parents=True)
    events = schedule(half, arms)
    group = parent.cgroup_root()
    host = parent.host_receipt(group)
    host['schema'] = 'n83.highk-tiebreak-host/v1'
    host['parent_instance_id'] = host['instance_id']
    host['instance_id'] = f"{host['instance_id']}-highk-half{half}"
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
        'schema': 'n83.highk-tiebreak-source/v1',
        'hashes': hashes, 'source_sha256': SOURCE_SHA256,
        'binary_sha256': parent.digest(ic_copy),
        'rho_binary_sha256': parent.digest(rho_copy),
        'runner_sha256': parent.digest(Path(__file__)),
        'targets_sha256': TARGETS_SHA256,
        'answers_sha256': ANSWERS_SHA256,
        'fixture_manifest_sha256': parent.digest(HERE / 'fixture_manifest.json'),
        'half': half, 'arms': arms, 'choice': choice, 'events': events,
        'rho_seed_rule': '2026102100+3*point_index+fixed_candidate_index',
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
            environment_flags = {'PATH': os.environ.get('PATH', '/usr/bin:/bin')}
            limit = RHO_CAP_SECONDS
        else:
            arm = by_name[candidate]
            command = [str(ic_copy), f'construct:83:1:{arm[1]}', str(target),
                       '29', str(folder / 'targets_result.jsonl')]
            environment_flags = environment(folder, arm)
            limit = cap(arm)
        elapsed = time.monotonic() - started
        if elapsed + limit + 30 > GLOBAL_PRODUCER_CAP_SECONDS:
            parent.write_json(output / 'producer_stop.json', {
                'status': 'GLOBAL_CAP_BEFORE_NEXT_ARM',
                'elapsed_seconds': elapsed, 'next_event': event,
                'global_producer_cap_seconds': GLOBAL_PRODUCER_CAP_SECONDS})
            break
        print(f'starting p{index:02d}/{name} at {parent.now()}', flush=True)
        usage = holdout.run_with_telemetry(folder, command, environment_flags,
                                          limit, group)
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
    parser.add_argument('--validated', type=Path)
    parser.add_argument('--selection', type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--rho-binary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    hashes, _ = fixture()
    if args.fixture_preflight:
        print(json.dumps({'hashes': hashes, 'events': schedule(
            args.half, (WIDE3, ROOT4, prior.K3072))}, sort_keys=True))
        return 0
    if not args.validated or not args.selection:
        parser.error('--validated and --selection are required')
    arms, choice = choose_anchor(args.validated.resolve(),
                                 args.selection.resolve())
    if args.preflight:
        print(json.dumps({'hashes': hashes, 'choice': choice,
                          'arms': arms, 'events': schedule(args.half, arms)},
                         sort_keys=True))
        return 0
    if not args.binary or not args.rho_binary or not args.output:
        parser.error('--binary, --rho-binary and --output are required')
    return run(args.binary.resolve(), args.rho_binary.resolve(),
               args.output.resolve(), args.half, args.validated.resolve(),
               args.selection.resolve())


if __name__ == '__main__':
    raise SystemExit(main())
