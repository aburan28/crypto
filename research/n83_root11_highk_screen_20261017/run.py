#!/usr/bin/env python3
"""Run an ordered, frozen high-K N83 root11 factor-base screen."""
import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FIXTURE = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_targets.jsonl'
ANSWERS = ROOT / 'research/n83_root11_finegrid_holdout_20261017/local_batch/batch_answers.jsonl'
SOURCE = ROOT / 'examples/koblitz_orbit_dlp_fast.rs'
REFERENCE_ENV = ROOT / 'research/n83_root11_dynamic_20261016/run.py'
GROUP = ROOT / 'research/n83_indexed5_20261013/analyze.py'
BASELINE = ROOT / 'research/n83_root11_finegrid_holdout_20261017/results/batch_k3072_result.json'
ARMS = (3584, 4096, 4608, 5120)
CAP = 1800

FLAGS = {
    'KIC_INDEXED4_FIXED31': '0', 'KIC_INDEXED5_FIXED30': '0',
    'KIC_INDEXED5_LOAD75': '1', 'KIC_INDEXED_ROOT4': '0',
    'KIC_INDEXED_ROOT5': '1', 'KIC_PACKED_ROOT_TABLE': '1',
    'KIC_PACKED_ROOT_TABLE_DENSE': '0',
    'KIC_PACKED_ROOT_TABLE_LOAD75': '0', 'KIC_RANK_THREADS': '4',
    'KIC_ROOT_BLOOM': '0', 'KIC_ROOT_ONLY11': '1',
    'KIC_ROOT_ONLY21': '0', 'KIC_ROOT_ONLY22': '0',
    'KIC_SPARSE_CANONICAL': '1', 'KIC_STATE26': '0',
    'KIC_STATE28': '0', 'KIC_SYMMETRIC_S3': '1',
    'KIC_TARGET_THREADS': '1',
}


def environment(folder):
    return {'PATH': os.environ.get('PATH', '/usr/bin:/bin'),
            **FLAGS, 'KIC_DUMP_BASE': str(folder / 'base_dump.jsonl')}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'existing receipt differs: {path}')
    path.write_text(blob)


def frozen(binary):
    plan = json.loads((HERE / 'plan.json').read_text())
    paths = {'source': SOURCE, 'binary': binary, 'targets': FIXTURE,
             'answers': ANSWERS,
             'reference_environment': REFERENCE_ENV,
             'group': GROUP, 'baseline_k3072': BASELINE,
             'protocol': HERE / 'PROTOCOL.md', 'runner': HERE / 'run.py',
             'verifier': HERE / 'verify.py'}
    if (plan['schema'] != 'n83.root11-highk-screen-plan/v1' or
            plan['binary_path'] != str(binary) or
            plan['k_values'] != list(ARMS) or plan['cap_seconds_per_arm'] != CAP or
            any(digest(path) != plan['sha256'][name]
                for name, path in paths.items())):
        raise ValueError('frozen source, inputs, code, binary or plan differs')
    return plan


def run(binary, root):
    plan = frozen(binary)
    if root.exists():
        raise ValueError('raw root already exists; inspect before retrying')
    root.mkdir(parents=True)
    write_once(root / 'plan_receipt.json', {
        'schema': 'n83.root11-highk-screen-plan-receipt/v1',
        'plan_sha256': digest(HERE / 'plan.json'),
        'plan': plan,
        'binary_path': str(binary),
        'raw_root': str(root)})
    for k in ARMS:
        folder = root / f'k{k}'
        folder.mkdir()
        env = environment(folder)
        command = [str(binary), f'construct:83:1:{k}', str(FIXTURE), '29',
                   str(folder / 'targets_result.jsonl')]
        write_once(folder / 'command.json', {
            'schema': 'n83.root11-highk-screen-command/v1', 'k': k,
            'argv': command, 'environment': env,
            'cap_seconds': CAP,
            'plan_sha256': digest(HERE / 'plan.json')})
        print(f'K{k}: starting at {time.time():.3f}', flush=True)
        start = time.monotonic()
        rss_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        try:
            process = subprocess.run(command, cwd=ROOT, env=env,
                                     capture_output=True, text=True,
                                     timeout=CAP, check=False)
            status = 'EXIT_0' if process.returncode == 0 else 'PROCESS_FAILED'
            code, stdout, stderr = process.returncode, process.stdout, process.stderr
        except subprocess.TimeoutExpired as error:
            status, code = 'TIMEOUT', None
            stdout = (error.stdout or b'').decode(errors='replace')
            stderr = (error.stderr or b'').decode(errors='replace')
        wall = time.monotonic() - start
        rss_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        (folder / 'stdout.jsonl').write_text(stdout)
        (folder / 'stderr.txt').write_text(stderr)
        try:
            summary = json.loads(stdout) if status == 'EXIT_0' else None
        except json.JSONDecodeError:
            status, summary = 'MALFORMED_SUMMARY', None
        write_once(folder / 'run_receipt.json', {
            'schema': 'n83.root11-highk-screen-run/v1', 'k': k,
            'producer_status': status, 'exit_code': code,
            'outer_wall_seconds': wall,
            'children_maxrss_before_raw': rss_before,
            'children_maxrss_after_raw': rss_after,
            'children_maxrss_unit': 'bytes' if sys.platform == 'darwin' else 'kib',
            'stdout_sha256': digest(folder / 'stdout.jsonl'),
            'stderr_sha256': digest(folder / 'stderr.txt'),
            'base_dump_sha256': digest(folder / 'base_dump.jsonl') if
                (folder / 'base_dump.jsonl').is_file() else None,
            'target_results_sha256': digest(folder / 'targets_result.jsonl') if
                (folder / 'targets_result.jsonl').is_file() else None,
            'summary': summary})
        print(f'K{k}: {status} after {wall:.3f}s', flush=True)
        if status != 'EXIT_0':
            break
        verifier_spec = importlib.util.spec_from_file_location(
            'n83_highk_verifier', HERE / 'verify.py')
        verifier = importlib.util.module_from_spec(verifier_spec)
        verifier_spec.loader.exec_module(verifier)
        try:
            result = verifier.verify_one(folder, k, plan)
        except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
            result = {'schema': 'n83.root11-highk-screen-analysis/v1',
                      'status': 'VERIFICATION_FAILED_OR_CENSORED',
                      'k': k, 'error': str(error)}
        write_once(folder / 'analysis.json', result)
        print(f'K{k}: {result["status"]}', flush=True)
        if result['status'] != 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY':
            break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.binary.resolve(), args.output.resolve())


if __name__ == '__main__':
    main()
