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

    def run(label: str, argv: list[str], timeout: int = 120) -> bytes:
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
            run('git_cap_rev_parse', ['git', 'rev-parse', 'HEAD'])
            run('git_cap_base_object', [
                'git', 'cat-file', '-t', frozen['base_main_head'],
            ])
            run('git_cap_checkout_ancestry', [
                'git', 'merge-base', '--is-ancestor',
                frozen['base_main_head'], args.expected_head,
            ])
            check = json.loads(run('preparation_and_v2_hash_gate', [
                sys.executable, str(HERE / 'ci_replay.py'),
            ]))
            if (check['decision'] != 'PASS_HASH_AND_ARCHIVE_REPLAY' or
                    check['release_status'] != frozen['status'] or
                    check['evidence'] is not None):
                raise RuntimeError('held preparation/static gate did not pass')
            for number, key in ((804, 'first'), (831, 'preparation')):
                parent = json.loads(run(f'gh_parent_{number}', [
                    'gh', 'pr', 'view', str(number), '--repo', 'aburan28/crypto',
                    '--json', 'state,headRefOid,mergeCommit',
                ]))
                if (parent['state'] != 'MERGED' or
                        parent['headRefOid'] != frozen[key]['pr_head'] or
                        parent['mergeCommit']['oid'] != frozen[key]['merge_commit']):
                    raise RuntimeError(f'parent PR {number} differs from freeze')
            event = json.loads(Path(os.environ['GITHUB_EVENT_PATH']).read_text())
            pr_number = event['pull_request']['number']
            if event['pull_request']['head']['sha'] != args.expected_head:
                raise RuntimeError('event PR head differs from checkout')
            live = json.loads(run('gh_live_pr', [
                'gh', 'pr', 'view', str(pr_number),
                '--repo', 'aburan28/crypto',
                '--json', 'state,isDraft,baseRefName,headRefName,headRefOid',
            ]))
            if (live['state'] != 'OPEN' or live['baseRefName'] != 'main' or
                    live['headRefOid'] != args.expected_head or
                    live['headRefName'] != frozen['release_branch']):
                raise RuntimeError('live PR differs from held exact head')
            run_id = int(os.environ['GITHUB_RUN_ID'])
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
            receipt['decision'] = 'PASS_HARMLESS_TOY_CAP_CONTROL'
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
                                       'SKIPPED_AFTER_ARCHIVED_OUTCOME') else 1


if __name__ == '__main__':
    raise SystemExit(main())
