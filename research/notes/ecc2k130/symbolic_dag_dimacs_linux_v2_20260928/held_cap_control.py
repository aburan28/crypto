#!/usr/bin/env python3
"""Harmless held Ubuntu control of release prerequisites under toy RLIMIT_AS.

This script never imports or launches the measured toy/panel/n131 producer.
Each command runs in its own child with the same hard 512-MiB address-space
limit that the first measured child would inherit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CAP = 536870912


def limits() -> None:
    resource.setrlimit(resource.RLIMIT_AS, (CAP, CAP))


def classify_capped_git(row: dict) -> str:
    if row['exit_code'] == 0:
        return 'PASS'
    stderr = row.get('stderr_tail', '')
    if (row['exit_code'] == 128 and 'packfile ' in stderr and
            'cannot be mapped' in stderr and 'Cannot allocate memory' in stderr):
        return 'PACK_MMAP_REFUSAL'
    raise RuntimeError('unexpected capped Git failure: ' + row['label'])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--expected-head', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    receipt = {
        'schema': 'k0-dag-linux-v2-held-toy-cap-control-v1',
        'decision': 'FAIL_HARMLESS_CONTROL',
        'release_status': frozen['status'],
        'measured_children_started': 0,
        'cap_bytes': CAP,
        'checkout_head': args.expected_head,
        'runner_image_os': os.environ.get('ImageOS'),
        'runner_image_version': os.environ.get('ImageVersion'),
        'commands': [],
    }

    def run(label: str, argv: list[str], timeout: int = 120,
            *, allow_failure: bool = False) -> bytes:
        start = time.monotonic()
        row: dict = {'label': label, 'argv': argv, 'cap_bytes': CAP}
        receipt['commands'].append(row)
        try:
            completed = subprocess.run(
                argv, cwd=ROOT, preexec_fn=limits, capture_output=True,
                timeout=timeout)
            row.update({
                'exit_code': completed.returncode,
                'stdout_bytes': len(completed.stdout),
                'stderr_bytes': len(completed.stderr),
                'stdout_sha256': hashlib.sha256(completed.stdout).hexdigest(),
                'stderr_sha256': hashlib.sha256(completed.stderr).hexdigest(),
            })
            if completed.returncode != 0:
                token = os.environ.get('GH_TOKEN', '')
                for stream in ('stdout', 'stderr'):
                    tail = getattr(completed, stream).decode('utf-8', 'replace')[-3000:]
                    row[stream + '_tail'] = tail.replace(token, '[REDACTED]') if token else tail
                if not allow_failure:
                    raise RuntimeError(f'{label} exited {completed.returncode}')
            return completed.stdout
        finally:
            row['wall_seconds'] = time.monotonic() - start

    try:
        if (HERE / frozen['archive_relative'] / 'receipt.json').is_file():
            # Outcome CI already replays the archive in its uncapped step.
            receipt['decision'] = 'SKIPPED_AFTER_ARCHIVED_OUTCOME'
        elif (sys.platform != 'linux' or platform.machine() != 'x86_64' or
                os.environ.get('ImageOS') != 'ubuntu24' or
                not os.environ.get('ImageVersion') or
                os.environ.get('GITHUB_ACTIONS') != 'true'):
            raise RuntimeError('requires hosted Ubuntu 24.04 x86-64 Actions image')
        else:
            if CAP != frozen['caps']['toy']['rss_bytes']:
                raise RuntimeError('toy hard cap changed')
            if subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                       text=True).strip() != args.expected_head:
                raise RuntimeError('checkout differs from event head')
            observed = json.loads(run('hard_limit_observation', [
                sys.executable, '-c',
                'import json,resource; print(json.dumps(resource.getrlimit(resource.RLIMIT_AS)))',
            ]))
            if observed != [CAP, CAP]:
                raise RuntimeError('hard address-space cap was not inherited')
            child = json.loads(run('capped_child_local_byte_gate', [
                sys.executable, str(HERE / 'v2_child.py'), '--gate-only',
            ]))
            if child['decision'] != 'HASH_ONLY_NO_MEASURED_CHILD':
                raise RuntimeError('actual child byte gate did not pass under toy cap')
            if run('git_cap_rev_parse', ['git', 'rev-parse', 'HEAD']).decode().strip() != args.expected_head:
                raise RuntimeError('capped Git head differs from event head')
            base_type = run('git_cap_base_object', [
                'git', 'cat-file', '-t', frozen['base_main_head'],
            ], allow_failure=True)
            base_row = receipt['commands'][-1]
            base_class = classify_capped_git(base_row)
            if base_class == 'PASS' and base_type.strip() != b'commit':
                raise RuntimeError('capped Git base object has wrong type')
            receipt['capped_git_base_object_classification'] = base_class
            run('git_cap_checkout_ancestry', [
                'git', 'merge-base', '--is-ancestor',
                frozen['base_main_head'], args.expected_head,
            ], allow_failure=True)
            ancestry_class = classify_capped_git(receipt['commands'][-1])
            receipt['capped_git_ancestry_classification'] = ancestry_class
            capped_refusal = 'PACK_MMAP_REFUSAL' in (base_class, ancestry_class)
            receipt['capped_git_classification'] = (
                'PACK_MMAP_REFUSAL' if capped_refusal else 'PASS')
            check_raw = run('preparation_and_v2_hash_gate', [
                sys.executable, str(HERE / 'ci_replay.py'),
            ], allow_failure=capped_refusal)
            check_row = receipt['commands'][-1]
            if check_row['exit_code'] == 0:
                check = json.loads(check_raw)
                if (check['decision'] != 'PASS_HASH_AND_ARCHIVE_REPLAY' or
                        check['release_status'] != frozen['status'] or
                        check['evidence'] is not None):
                    raise RuntimeError('preparation/static gate did not pass')
            elif (not capped_refusal or check_row['exit_code'] != 1 or
                  'AssertionError: v2 branch does not descend' not in
                  check_row.get('stderr_tail', '')):
                raise RuntimeError('capped static gate failed for another reason')
            prep_raw = run('preparation_archive_gate', [
                sys.executable,
                str(HERE.parent / 'symbolic_dag_dimacs_linux_attempt2_20260926' /
                    'verify_preparation.py'),
            ], allow_failure=capped_refusal)
            prep_row = receipt['commands'][-1]
            if prep_row['exit_code'] == 0:
                if b'PASS_HARMLESS_PREPARATION_ONLY' not in prep_raw:
                    raise RuntimeError('preparation archive gate did not pass')
            elif (not capped_refusal or
                  'NOT_ADMITTED: first-failure record is not an ancestor' not in
                  prep_row.get('stderr_tail', '')):
                raise RuntimeError('capped preparation gate failed for another reason')
            event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())
            pr_number = event['pull_request']['number']
            if event['pull_request']['head']['sha'] != args.expected_head:
                raise RuntimeError('event PR head differs from checkout')
            run_id = int(os.environ['GITHUB_RUN_ID'])
            first_raw = run('gh_parent_804', [
                'gh', 'pr', 'view', '804', '--repo', 'aburan28/crypto',
                '--json', 'state,headRefOid,mergeCommit',
            ], allow_failure=True)
            first_row = receipt['commands'][-1]
            if (first_row['exit_code'] == 2 and
                    'fatal error: failed to reserve page summary memory' in
                    first_row.get('stderr_tail', '')):
                receipt['capped_gh_classification'] = 'GO_PAGE_SUMMARY_REFUSAL'
                run('gh_actions_current_run', [
                    'gh', 'api', f'repos/aburan28/crypto/actions/runs/{run_id}',
                ], allow_failure=True)
                api_row = receipt['commands'][-1]
                if (api_row['exit_code'] != 2 or
                        'fatal error: failed to reserve page summary memory' not in
                        api_row.get('stderr_tail', '')):
                    raise RuntimeError('capped gh API failed for another reason')
            elif first_row['exit_code'] == 0:
                receipt['capped_gh_classification'] = 'PASS'
                first = json.loads(first_raw)
                if (first['state'] != 'MERGED' or
                        first['headRefOid'] != frozen['first']['pr_head'] or
                        first['mergeCommit']['oid'] != frozen['first']['merge_commit']):
                    raise RuntimeError('parent PR 804 differs from freeze')
                parent = json.loads(run('gh_parent_831', [
                    'gh', 'pr', 'view', '831', '--repo', 'aburan28/crypto',
                    '--json', 'state,headRefOid,mergeCommit',
                ]))
                if (parent['state'] != 'MERGED' or
                        parent['headRefOid'] != frozen['preparation']['pr_head'] or
                        parent['mergeCommit']['oid'] !=
                        frozen['preparation']['merge_commit']):
                    raise RuntimeError('parent PR 831 differs from freeze')
                live = json.loads(run('gh_live_pr', [
                    'gh', 'pr', 'view', str(pr_number),
                    '--repo', 'aburan28/crypto',
                    '--json', 'state,isDraft,baseRefName,headRefName,headRefOid',
                ]))
                if (live['state'] != 'OPEN' or live['baseRefName'] != 'main' or
                        live['headRefOid'] != args.expected_head or
                        live['headRefName'] != frozen['release_branch']):
                    raise RuntimeError('live PR differs from held exact head')
                current = json.loads(run('gh_actions_current_run', [
                    'gh', 'api', f'repos/aburan28/crypto/actions/runs/{run_id}',
                ]))
                if current['id'] != run_id or not isinstance(current['workflow_id'], int):
                    raise RuntimeError('Actions current-run API result drifted')
                listed = run('gh_actions_workflow_runs', [
                    'gh', 'api', '--paginate', '--jq', '.workflow_runs[] | @json',
                    'repos/aburan28/crypto/actions/workflows/' +
                    str(current['workflow_id']) + '/runs?event=pull_request&per_page=100',
                ])
                if run_id not in [json.loads(line)['id'] for line in listed.splitlines()]:
                    raise RuntimeError('Actions workflow run listing omitted current run')
            else:
                raise RuntimeError('unexpected capped gh prerequisite failure')
            if not capped_refusal:
                if run('git_clean_checkout', ['git', 'status', '--porcelain']).strip():
                    raise RuntimeError('held checkout is dirty')
                run('git_fetch_main', ['git', 'fetch', '--quiet', 'origin', 'main'])
                run('git_main_ancestry', [
                    'git', 'merge-base', '--is-ancestor',
                    frozen['base_main_head'], 'origin/main',
                ])
                run('git_guarded_diff', [
                    'git', 'diff', '--name-only',
                    frozen['base_main_head'] + '..origin/main', '--',
                    *frozen['upstream_guard_paths'],
                ])
                run('git_guarded_full_history', [
                    'git', 'log', '--full-history', '--format=%H',
                    frozen['base_main_head'] + '..origin/main', '--',
                    *frozen['upstream_guard_paths'],
                ])
            cap_limited = (capped_refusal or
                           receipt['capped_gh_classification'] == 'GO_PAGE_SUMMARY_REFUSAL')
            receipt['decision'] = (
                'PASS_CHILD_BYTE_GATE_WITH_CAPPED_PREREQUISITE_REFUSALS'
                if cap_limited else 'PASS_HARMLESS_TOY_CAP_CONTROL')
    except Exception as exc:
        receipt['error'] = f'{type(exc).__name__}: {exc}'
    out.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({
        'decision': receipt['decision'],
        'commands': len(receipt['commands']),
        'measured_children_started': 0,
        'error': receipt.get('error'),
    }, sort_keys=True))
    return 0 if receipt['decision'] in ('PASS_HARMLESS_TOY_CAP_CONTROL',
                                       'PASS_CHILD_BYTE_GATE_WITH_CAPPED_PREREQUISITE_REFUSALS',
                                       'SKIPPED_AFTER_ARCHIVED_OUTCOME') else 1


if __name__ == '__main__':
    raise SystemExit(main())
