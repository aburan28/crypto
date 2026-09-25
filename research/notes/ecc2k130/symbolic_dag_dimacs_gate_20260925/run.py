#!/usr/bin/env python3
"""Fail-closed release gate and cold three-phase bounded experiment runner."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

from bounded import run_child, sha

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DOMAIN = 'k0-symbolic-dag-dimacs-gate-v1'
PARENT_REL = 'research/notes/ecc2k130/symbolic_dag_fullpoint_20260925'


def _git(*args: str) -> str:
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def _main_blob_sha(relative: str) -> str:
    data = subprocess.check_output(['git', 'show', f'origin/main:{relative}'], cwd=ROOT)
    return hashlib.sha256(data).hexdigest()


def release_gate(frozen: dict) -> dict:
    """Impossible to run a stacked outcome while parent #802 is unmerged."""
    if frozen['release_main_head'] is None:
        raise RuntimeError('NOT_ADMITTED: #802 must merge, then rebase and re-freeze')
    if sha(Path(frozen['cadical_path'])) != frozen['cadical_sha256']:
        raise RuntimeError('NOT_ADMITTED: pinned CaDiCaL binary drifted')
    pr = json.loads(subprocess.check_output(
        ['gh', 'pr', 'view', '802', '--json', 'state,mergeCommit'], cwd=ROOT, text=True))
    if pr['state'] != 'MERGED' or not pr['mergeCommit']:
        raise RuntimeError('NOT_ADMITTED: parent PR #802 is not merged')
    subprocess.run(['git', 'fetch', 'origin', 'main'], cwd=ROOT, check=True)
    main_head = _git('rev-parse', 'origin/main')
    if main_head != frozen['release_main_head']:
        raise RuntimeError('NOT_ADMITTED: main moved after the release freeze')
    subprocess.run(['git', 'merge-base', '--is-ancestor', main_head, 'HEAD'],
                   cwd=ROOT, check=True)
    if _main_blob_sha(f'{PARENT_REL}/FROZEN.json') != frozen['parent_freeze_sha256']:
        raise RuntimeError('NOT_ADMITTED: merged parent freeze bytes drifted')
    if _main_blob_sha(f'{PARENT_REL}/evidence/producer/rows.jsonl.gz') != frozen['parent_rows_sha256']:
        raise RuntimeError('NOT_ADMITTED: merged parent raw rows drifted')
    return {'parent_pr_state': pr['state'], 'parent_merge_commit': pr['mergeCommit']['oid'],
            'release_main_head': main_head}


def main() -> int:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(HERE / 'ci_replay.py')], cwd=ROOT, check=True)
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    gate = release_gate(frozen)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    attempts = []
    specs = (
        ('toy', frozen['toy_external_wall_cap_seconds'], frozen['toy_rss_cap_bytes']),
        ('panel', frozen['panel_external_wall_cap_seconds'], frozen['panel_rss_cap_bytes']),
        ('n131', frozen['n131_external_wall_cap_seconds'], frozen['n131_rss_cap_bytes']),
    )
    for phase, wall_cap, rss_cap in specs:
        command = [sys.executable, str(HERE / 'produce.py'), '--phase', phase,
                   '--out', str(out / phase)]
        raw = run_child(command, cwd=ROOT, stdout=out / f'{phase}.stdout.txt',
                        stderr=out / f'{phase}.stderr.txt', wall_cap=wall_cap,
                        rss_cap=rss_cap,
                        watched_file=out / phase / 'single_edge.cnf' if phase == 'n131' else None,
                        watched_file_cap=frozen['n131_cnf_byte_cap'] if phase == 'n131' else None,
                        file_cap=frozen['n131_cnf_byte_cap'] if phase == 'n131' else None)
        result_path = out / phase / 'result.json'
        raw['phase'] = phase
        raw['result_sha256'] = sha(result_path)
        result = None
        if result_path.is_file():
            try:
                result = json.loads(result_path.read_text())
            except ValueError as exc:
                raw['result_parse_error'] = str(exc)
        raw['reported_decision'] = result.get('decision') if result else None
        raw['reported_cpu_seconds'] = result.get('cpu_seconds') if result else None
        raw['reported_peak_rss_bytes'] = result.get('peak_rss_bytes') if result else None
        attempts.append(raw)
        if (raw['exit_code'] != 0 or raw['stop_reason'] is not None or
                raw['wall_seconds'] > wall_cap or
                raw['sampled_peak_rss_bytes'] > rss_cap or
                result is None or result.get('decision') != 'PASS' or
                result.get('peak_rss_bytes', rss_cap + 1) > rss_cap):
            break
    manifest = []
    for artifact in sorted(path for path in out.rglob('*') if path.is_file()):
        manifest.append({'path': artifact.relative_to(out).as_posix(),
                         'bytes': artifact.stat().st_size, 'sha256': sha(artifact)})
    (out / 'MANIFEST.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    caps = {phase: (wall, rss) for phase, wall, rss in specs}
    receipt = {'domain': DOMAIN, 'freeze_sha256': sha(HERE / 'FROZEN.json'),
               'release_gate': gate, 'attempts': attempts,
               'manifest_sha256': sha(out / 'MANIFEST.json'),
               'decision': 'PASS' if len(attempts) == 3 and all(
                   item['exit_code'] == 0 and item['stop_reason'] is None and
                   item['reported_decision'] == 'PASS' and
                   item['reported_peak_rss_bytes'] is not None and
                   item['reported_peak_rss_bytes'] <= caps[item['phase']][1] and
                   item['sampled_peak_rss_bytes'] <= caps[item['phase']][1] and
                   item['wall_seconds'] <= caps[item['phase']][0]
                   for item in attempts) else 'FAIL_OR_CENSORED'}
    (out / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'decision': receipt['decision'],
                      'phases': [item['phase'] for item in attempts]}, sort_keys=True))
    return 0 if receipt['decision'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
