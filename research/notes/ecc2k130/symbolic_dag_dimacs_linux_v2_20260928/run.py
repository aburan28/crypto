#!/usr/bin/env python3
"""Held, one-shot Linux release gate and bounded three-phase supervisor."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from ci_replay import check_static

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FIRST = HERE.parent / 'symbolic_dag_dimacs_gate_20260925'
PREPARATION = HERE.parent / 'symbolic_dag_dimacs_linux_attempt2_20260926'
WORKFLOW = 'ecc2k130-symbolic-dag-dimacs-linux-v2-measure.yml'
LABEL = 'ecc2k130-dag-linux-v2-measure-once'
HEX40 = re.compile(r'[0-9a-f]{40}\Z')


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def git(*args: str) -> str:
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError('NOT_ADMITTED: ' + message)


def _gh_json(*args: str) -> dict:
    return json.loads(subprocess.check_output(['gh', *args], cwd=ROOT, text=True))


def _one_shot_run(frozen: dict, run_id: int) -> dict:
    """Refuse a second labeled workflow event even if its first run failed.

    The measured workflow is separate from hash-only CI and has only a
    pull_request.labeled trigger.  Any previous event for this PR branch is
    conservatively treated as an attempted dispatch, regardless of outcome.
    """
    current = _gh_json('api', f'repos/aburan28/crypto/actions/runs/{run_id}')
    require(current.get('id') == run_id and
            current.get('path') == '.github/workflows/' + WORKFLOW and
            current.get('head_branch') == frozen['release_branch'] and
            current.get('event') == 'pull_request' and
            isinstance(current.get('workflow_id'), int),
            'current labeled workflow identity drift')
    endpoint = ('repos/aburan28/crypto/actions/workflows/' +
                str(current['workflow_id']) + '/runs?event=pull_request&per_page=100')
    raw = subprocess.check_output(
        ['gh', 'api', '--paginate', '--jq', '.workflow_runs[] | @json', endpoint],
        cwd=ROOT, text=True)
    runs = [json.loads(line) for line in raw.splitlines() if line.strip()]
    pr_number = frozen['release_pr_number']
    branch = frozen['release_branch']
    matching = []
    for row in runs:
        numbers = {pr.get('number') for pr in row.get('pull_requests', [])}
        if pr_number in numbers or row.get('head_branch') == branch:
            matching.append(row)
    require(len(matching) == 1 and matching[0].get('id') == run_id,
            'measured workflow already has a labeled run for this PR')
    require(matching[0].get('event') == 'pull_request', 'wrong workflow event')
    return {'workflow': WORKFLOW, 'workflow_id': current['workflow_id'],
            'run_id': run_id,
            'matching_labeled_runs': [row['id'] for row in matching]}


def _event_context(frozen: dict, expected_head: str) -> dict:
    require(os.environ.get('GITHUB_ACTIONS') == 'true', 'not a hosted Actions run')
    require(os.environ.get('GITHUB_REPOSITORY') == 'aburan28/crypto', 'wrong repository')
    require(os.environ.get('GITHUB_EVENT_NAME') == 'pull_request', 'wrong event')
    require(os.environ.get('GITHUB_RUN_ATTEMPT') == '1', 'workflow rerun refused')
    workflow_ref = os.environ.get('GITHUB_WORKFLOW_REF', '')
    require('/.github/workflows/' + WORKFLOW + '@' in workflow_ref,
            'wrong measured workflow file')
    event_path = Path(os.environ.get('GITHUB_EVENT_PATH', ''))
    require(event_path.is_file(), 'missing immutable event payload')
    event = json.loads(event_path.read_text())
    pr = event.get('pull_request') or {}
    require(event.get('action') == 'labeled' and
            (event.get('label') or {}).get('name') == LABEL,
            'wrong labeled event')
    require((event.get('repository') or {}).get('full_name') == 'aburan28/crypto',
            'wrong event repository')
    require(pr.get('number') == frozen['release_pr_number'] and
            pr.get('base', {}).get('ref') == 'main' and
            pr.get('head', {}).get('sha') == expected_head and
            pr.get('head', {}).get('ref') == frozen['release_branch'] and
            pr.get('draft') is False,
            'event PR number/base/head/draft mismatch')
    run_id = int(os.environ.get('GITHUB_RUN_ID', '0'))
    require(run_id > 0, 'missing Actions run ID')
    return {'run_id': run_id, 'label': LABEL,
            'event_head': pr['head']['sha'], 'pr_number': pr['number']}


def _target_cap_probe(frozen: dict, expected_head: str) -> dict:
    value = os.environ.get('V2_CAP_PROBE_PATH', '')
    path = Path(value) if value else Path('/nonexistent-v2-cap-probe')
    require(path.is_file(), 'fresh target-image hard-cap receipt is missing')
    receipt = json.loads(path.read_text())
    require(receipt['schema'] == 'k0-dag-dimacs-linux-rlimit-probe-v1' and
            receipt['decision'] == 'PASS_HARMLESS_CAPS_ONLY' and
            receipt['checkout_head'] == expected_head and
            receipt['machine'] == 'x86_64' and
            receipt['runner_image_os'] == os.environ.get('ImageOS') == 'ubuntu24' and
            receipt['runner_image_version'] == os.environ.get('ImageVersion') and
            receipt['python'].startswith('3.12.'),
            'fresh target-image cap probe provenance')
    caps = [frozen['caps']['toy']['rss_bytes'],
            frozen['caps']['panel']['query_rss_bytes'],
            frozen['caps']['panel']['rss_bytes']]
    require(len(receipt['caps']) == len(caps), 'fresh cap observation count')
    for cap, row in zip(caps, receipt['caps']):
        require(row['cap_bytes'] == cap and
                row['true_exit_code'] == row['python_exit_code'] == 0 and
                row['observation']['limit'] == [cap, cap] and
                row['observation']['over_cap_allocation'] == 'REJECTED' and
                row['observation']['exception'] == 'OSError',
                'fresh hard RLIMIT_AS control failed')
    return {'receipt_sha256': sha(path),
            'stdout_sha256': sha(path.parent / 'target_rlimit_probe.stdout.txt'),
            'stderr_sha256': sha(path.parent / 'target_rlimit_probe.stderr.txt'),
            'runner_image_os': receipt['runner_image_os'],
            'runner_image_version': receipt['runner_image_version'],
            'caps_bytes': caps}


def release_gate(frozen: dict, expected_head: str) -> dict:
    """Validate exact source, live PRs, main lineage, host and one-shot event."""
    check_static()
    require(frozen['status'] == 'RELEASED' and
            isinstance(frozen['release_pr_number'], int) and
            frozen['release_pr_number'] > 0 and
            HEX40.fullmatch(frozen['release_main_head'] or '') is not None,
            'v2 remains held without a reviewed exact-main release')
    require(HEX40.fullmatch(expected_head or '') is not None,
            'explicit reviewed head is missing')
    checkout_head = git('rev-parse', 'HEAD')
    require(checkout_head == expected_head, 'checkout differs from reviewed head')
    require(not git('status', '--porcelain'), 'checkout is dirty')
    context = _event_context(frozen, expected_head)
    require(sys.platform.startswith('linux') and os.uname().machine == 'x86_64',
            'wrong measured host class')
    require(sys.version_info[:2] == (3, 12), 'wrong Python minor version')
    target_probe = _target_cap_probe(frozen, expected_head)
    live = _gh_json('pr', 'view', str(frozen['release_pr_number']), '--json',
                    'state,isDraft,baseRefName,headRefName,headRefOid')
    require(live['state'] == 'OPEN' and live['isDraft'] is False and
            live['baseRefName'] == 'main' and
            live['headRefName'] == frozen['release_branch'] and
            live['headRefOid'] == expected_head,
            'live PR differs from reviewed exact head')
    parent_states = {}
    for number, key in ((804, 'first'), (831, 'preparation')):
        pr = _gh_json('pr', 'view', str(number), '--json',
                      'state,headRefOid,mergeCommit')
        require(pr['state'] == 'MERGED' and
                pr['headRefOid'] == frozen[key]['pr_head'] and
                pr['mergeCommit']['oid'] == frozen[key]['merge_commit'],
                f'parent PR #{number} identity drift')
        parent_states[str(number)] = pr['mergeCommit']['oid']
    subprocess.run(['git', 'fetch', '--quiet', 'origin', 'main'], cwd=ROOT,
                   check=True)
    main_head = git('rev-parse', 'origin/main')
    for ancestor, descendant in ((frozen['release_main_head'], main_head),
                                 (frozen['release_main_head'], checkout_head),
                                 (parent_states['804'], main_head),
                                 (parent_states['831'], main_head)):
        require(subprocess.run(['git', 'merge-base', '--is-ancestor', ancestor,
                                descendant], cwd=ROOT, capture_output=True).returncode == 0,
                'release/parent is not an ancestor of live main and checkout')
    guarded = frozen['upstream_guard_paths']
    changed = git('diff', '--name-only',
                  frozen['release_main_head'] + '..' + main_head, '--', *guarded)
    history = git('log', '--full-history', '--format=%H',
                  frozen['release_main_head'] + '..' + main_head, '--', *guarded)
    require(not changed and not history,
            'guarded upstream path changed, including a reverted edit')
    binary = HERE.parent / frozen['linux_binary_path']
    require(binary.is_file() and os.access(binary, os.X_OK),
            'pinned Linux solver is missing or nonexecutable')
    once = _one_shot_run(frozen, context['run_id'])
    return {'reviewed_head': expected_head, 'checkout_head': checkout_head,
            'freeze_sha256': sha(HERE / 'FROZEN.json'),
            'pr_head': live['headRefOid'], 'pr_number': frozen['release_pr_number'],
            'checkout_root': str(ROOT), 'release_main_head': frozen['release_main_head'],
            'dispatch_main_head': main_head, 'parent_merge_commits': parent_states,
            'upstream_changed_paths': [], 'upstream_guarded_history': [],
            'preparation_freeze_sha256': frozen['preparation']['freeze_sha256'],
            'linux_binary_sha256': frozen['linux_binary_sha256'],
            'event': context, 'one_shot': once,
            'target_cap_probe': target_probe}


def gate_only(expected_head: str) -> dict:
    frozen = check_static()
    require(git('rev-parse', 'HEAD') == expected_head, 'hash-only checkout mismatch')
    return {'decision': 'HASH_ONLY_NO_MEASURED_CHILD',
            'release_status': frozen['status'],
            'checkout_head': expected_head,
            'freeze_sha256': sha(HERE / 'FROZEN.json')}


def _manifest_and_receipt(out: Path, frozen: dict, gate: dict | None,
                          attempts: list[dict], error: str | None,
                          total_wall_seconds: float) -> dict:
    # These two files are intentionally excluded so the ledger has no cycle.
    artifacts = sorted(p for p in out.rglob('*') if p.is_file() and
                       p not in (out / 'MANIFEST.json', out / 'receipt.json'))
    manifest = [{'path': p.relative_to(out).as_posix(), 'bytes': p.stat().st_size,
                 'sha256': sha(p)} for p in artifacts]
    (out / 'MANIFEST.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    caps = frozen['caps']
    complete = (gate is not None and error is None and len(attempts) == 3 and
                all(row['exit_code'] == 0 and row['stop_reason'] is None and
                    row['reported_decision'] == 'PASS' and
                    row['reported_peak_rss_bytes'] is not None and
                    row['reported_peak_rss_bytes'] <= caps[row['phase']]['rss_bytes'] and
                    row['sampled_peak_rss_bytes'] <= caps[row['phase']]['rss_bytes'] and
                    row['wall_seconds'] <= caps[row['phase']]['wall_seconds']
                    for row in attempts))
    receipt = {'domain': frozen['domain'], 'freeze_sha256': sha(HERE / 'FROZEN.json'),
               'release_gate': gate, 'attempts': attempts,
               'manifest_sha256': sha(out / 'MANIFEST.json'),
               'decision': 'PASS' if complete else 'FAIL_OR_CENSORED',
               'error': error, 'runner_total_wall_seconds': total_wall_seconds}
    (out / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--gate-only', action='store_true')
    parser.add_argument('--expected-head', required=True)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--fallback-refusal', type=Path)
    args = parser.parse_args()
    if args.gate_only:
        if args.out:
            parser.error('--gate-only cannot create a measured archive')
        print(json.dumps(gate_only(args.expected_head), sort_keys=True))
        return 0
    if args.out is None:
        parser.error('measured dispatch requires --out')
    runner_start = time.monotonic()
    out = args.out.resolve()
    try:
        out.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        refusal = {'decision': 'PRE_DISPATCH_REFUSAL', 'reason':
                   f'OUTPUT_COLLISION_OR_CREATE_ERROR: {type(exc).__name__}: {exc}',
                   'attempts': []}
        if args.fallback_refusal:
            args.fallback_refusal.write_text(json.dumps(refusal, sort_keys=True, indent=2) + '\n')
        print(json.dumps(refusal, sort_keys=True), file=sys.stderr)
        return 2
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    gate = None
    attempts: list[dict] = []
    error = None
    try:
        if frozen['status'] == 'RELEASED':
            probe_path = out / 'target_rlimit_probe.json'
            os.environ['V2_CAP_PROBE_PATH'] = str(probe_path)
            try:
                cap_probe = subprocess.run(
                    [sys.executable, str(PREPARATION / 'harmless_cap_probe.py'),
                     '--out', str(probe_path)], cwd=ROOT, capture_output=True,
                    timeout=90)
            except subprocess.TimeoutExpired as exc:
                (out / 'target_rlimit_probe.stdout.txt').write_bytes(exc.stdout or b'')
                (out / 'target_rlimit_probe.stderr.txt').write_bytes(exc.stderr or b'')
                raise RuntimeError('fresh target-image cap probe exceeded 90 s') from exc
            (out / 'target_rlimit_probe.stdout.txt').write_bytes(cap_probe.stdout)
            (out / 'target_rlimit_probe.stderr.txt').write_bytes(cap_probe.stderr)
            require(cap_probe.returncode == 0,
                    'fresh target-image hard-cap probe failed')
        gate = release_gate(frozen, args.expected_head)
        dispatch_path = out / 'DISPATCH.json'
        dispatch_path.write_text(json.dumps(gate, sort_keys=True, indent=2) + '\n')
        dispatch_sha256 = sha(dispatch_path)
        sys.path.insert(0, str(FIRST))
        from bounded import run_child
        for phase in ('toy', 'panel', 'n131'):
            cap = frozen['caps'][phase]
            command = [sys.executable, str(HERE / 'v2_child.py'), '--phase', phase,
                       '--out', str(out / phase), '--expected-head', args.expected_head,
                       '--dispatch', str(dispatch_path),
                       '--dispatch-sha256', dispatch_sha256]
            raw = run_child(command, cwd=ROOT,
                            stdout=out / f'{phase}.stdout.txt',
                            stderr=out / f'{phase}.stderr.txt',
                            wall_cap=cap['wall_seconds'], rss_cap=cap['rss_bytes'],
                            watched_file=out / phase / 'single_edge.cnf'
                            if phase == 'n131' else None,
                            watched_file_cap=frozen['caps']['n131']['cnf_bytes']
                            if phase == 'n131' else None,
                            file_cap=frozen['caps']['n131']['cnf_bytes']
                            if phase == 'n131' else None)
            result_path = out / phase / 'result.json'
            result = None
            if result_path.is_file():
                try:
                    result = json.loads(result_path.read_text())
                except ValueError as exc:
                    raw['result_parse_error'] = str(exc)
            raw.update({'phase': phase, 'result_sha256': sha(result_path),
                        'reported_decision': result.get('decision') if result else None,
                        'reported_cpu_seconds': result.get('cpu_seconds') if result else None,
                        'reported_peak_rss_bytes': result.get('peak_rss_bytes')
                        if result else None})
            attempts.append(raw)
            if (raw['exit_code'] != 0 or raw['stop_reason'] is not None or
                    raw['wall_seconds'] > cap['wall_seconds'] or
                    raw['sampled_peak_rss_bytes'] > cap['rss_bytes'] or
                    not result or result.get('decision') != 'PASS' or
                    result.get('peak_rss_bytes', cap['rss_bytes'] + 1) >
                    cap['rss_bytes']):
                break
    except Exception as exc:
        error = f'{type(exc).__name__}: {exc}'
        (out / ('PRE_DISPATCH_REFUSAL.json' if gate is None else 'RUNNER_ERROR.json')).write_text(
            json.dumps({'decision': 'PRE_DISPATCH_REFUSAL' if gate is None else 'RUNNER_ERROR',
                        'reason': error, 'attempts_started': len(attempts)},
                       sort_keys=True, indent=2) + '\n')
    receipt = _manifest_and_receipt(out, frozen, gate, attempts, error,
                                    time.monotonic() - runner_start)
    print(json.dumps({'decision': receipt['decision'],
                      'phases': [row['phase'] for row in attempts],
                      'error': error,
                      'runner_total_wall_seconds': receipt['runner_total_wall_seconds']},
                     sort_keys=True))
    return 0 if receipt['decision'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
