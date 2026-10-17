#!/usr/bin/env python3
"""Run one frozen K3072 target-worker policy with an empty process cache."""

import argparse
import hashlib
import importlib.util
import json
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGETS = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl'
ANSWERS = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_answers.jsonl'
POLICIES = {'serial': (1, '0'), 'parallel_cancel': (2, '1')}

spec = importlib.util.spec_from_file_location(
    'n83_parallel_cancel_environment',
    ROOT / 'research/n83_root11_dynamic_20261016/run.py')
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'existing receipt differs: {path}')
    path.write_text(blob)


def frozen_files(binary):
    return {'source': ROOT / 'examples/koblitz_orbit_dlp_fast.rs',
             'binary': binary, 'targets': TARGETS, 'answers': ANSWERS,
             'planner': HERE / 'seal.py',
             'runner': HERE / 'run.py', 'verifier': HERE / 'verify.py',
             'group_verifier': ROOT / 'research/n83_indexed5_20261013/analyze.py',
             'candidate_materializer': HERE / 'materialize_candidates.py',
             'candidate_ids': HERE / 'candidate_ids.json',
             'candidate_serial': HERE / 'candidate_serial.json',
             'candidate_parallel_cancel': HERE / 'candidate_parallel_cancel.json',
             'environment': ROOT / 'research/n83_root11_dynamic_20261016/run.py',
             'protocol': HERE / 'PROTOCOL.md',
             'holdout_generator':
                 ROOT / 'research/n83_k3072_target_cancel_holdout_20261017/generate.py',
             'holdout_manifest':
                 ROOT / 'research/n83_k3072_target_cancel_holdout_20261017/manifest.json',
             'holdout_targets':
                 ROOT / 'research/n83_k3072_target_cancel_holdout_20261017/targets.jsonl',
             'holdout_answers':
                 ROOT / 'research/n83_k3072_target_cancel_holdout_20261017/answers.jsonl',
             'holdout_protocol':
                 ROOT / 'research/n83_k3072_target_cancel_holdout_20261017/PROTOCOL.md'}


def frozen(binary):
    plan = json.loads((HERE / 'plan.json').read_text())
    files = frozen_files(binary)
    if (plan['schema'] != 'n83.k3072-target-cancel-plan/v1' or
            plan['k'] != 3072 or plan['target_count'] != 128 or
            plan['binary_path'] != str(binary) or
            plan['candidate_ids'] !=
            json.loads((HERE / 'candidate_ids.json').read_text()) or
            plan['policies'] != {k: list(v) for k, v in POLICIES.items()} or
            plan['cap_seconds'] != 1200 or
            any(digest(path) != plan['sha256'][name]
                for name, path in files.items())):
        raise ValueError('source, binary, inputs or experiment plan differs')
    return plan


def run(policy, binary, output):
    plan = frozen(binary)
    if output.exists():
        raise ValueError('run directory already exists; inspect before retrying')
    output.mkdir(parents=True)
    threads, cancel = POLICIES[policy]
    candidate_id = json.loads((HERE / 'candidate_ids.json').read_text())[policy]['candidate_id']
    environment = dynamic.environment(output, 'root11_dynamic')
    environment['KIC_TARGET_THREADS'] = str(threads)
    environment['KIC_TARGET_CANCEL'] = cancel
    command = [str(binary), 'construct:83:1:3072', str(TARGETS), '29',
               str(output / 'targets_result.jsonl')]
    write_once(output / 'command.json', {
        'schema': 'n83.k3072-target-cancel-command/v1',
        'policy': policy, 'candidate_id': candidate_id,
        'argv': command, 'environment': environment,
        'plan_sha256': digest(HERE / 'plan.json'),
        'cap_seconds': plan['cap_seconds']})
    started = time.monotonic()
    try:
        process = subprocess.run(command, cwd=ROOT, env=environment,
                                 capture_output=True, text=True,
                                 timeout=plan['cap_seconds'], check=False)
        status = 'EXIT_0' if process.returncode == 0 else 'PROCESS_FAILED'
        exit_code, stdout, stderr = (process.returncode, process.stdout,
                                     process.stderr)
    except subprocess.TimeoutExpired as error:
        status, exit_code = 'TIMEOUT', None
        stdout = (error.stdout or b'').decode(errors='replace')
        stderr = (error.stderr or b'').decode(errors='replace')
    wall = time.monotonic() - started
    (output / 'stdout.jsonl').write_text(stdout)
    (output / 'stderr.txt').write_text(stderr)
    summary = None
    if status == 'EXIT_0':
        try:
            summary = json.loads(stdout)
        except json.JSONDecodeError:
            status = 'MALFORMED_SUMMARY'
    receipt = {'schema': 'n83.k3072-target-cancel-run/v1',
               'policy': policy, 'candidate_id': candidate_id,
               'producer_status': status,
               'exit_code': exit_code, 'outer_wall_seconds': wall,
               'source_sha256': plan['sha256']['source'],
               'binary_sha256': plan['sha256']['binary'],
               'stdout_sha256': digest(output / 'stdout.jsonl'),
               'stderr_sha256': digest(output / 'stderr.txt'),
               'target_result_sha256':
                   digest(output / 'targets_result.jsonl') if
                   (output / 'targets_result.jsonl').is_file() else None,
               'base_dump_sha256':
                   digest(output / 'base_dump.jsonl') if
                   (output / 'base_dump.jsonl').is_file() else None,
               'summary': summary}
    write_once(output / 'run_receipt.json', receipt)
    print(json.dumps({'policy': policy, 'producer_status': status,
                      'outer_wall_seconds': wall}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--policy', choices=POLICIES, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.policy, args.binary.resolve(), args.output.resolve())


if __name__ == '__main__':
    main()
