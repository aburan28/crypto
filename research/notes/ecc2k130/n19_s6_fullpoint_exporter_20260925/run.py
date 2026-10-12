#!/usr/bin/env python3
"""Post-merge fail-closed S6 all-tuple semantic admission supervisor."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from bounded import run_child, sha  # noqa: E402


def git(*args: str) -> str:
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def release_gate(frozen: dict) -> dict:
    if frozen['release_main_head'] is None:
        raise RuntimeError('HELD: #802/#804/#807 must merge, then rebase and re-freeze')
    subprocess.run([sys.executable, str(HERE / 'ci_replay.py')], cwd=ROOT, check=True)
    subprocess.run(['git', 'fetch', 'origin', 'main'], cwd=ROOT, check=True)
    main_head = git('rev-parse', 'origin/main')
    if main_head != frozen['release_main_head']:
        raise RuntimeError('HELD: main moved after exact-head release freeze')
    subprocess.run(['git', 'merge-base', '--is-ancestor', main_head, 'HEAD'], cwd=ROOT, check=True)
    parents = {}
    for number in frozen['required_merged_prs']:
        pr = json.loads(subprocess.check_output(
            ['gh', 'pr', 'view', str(number), '--json', 'state,mergeCommit'],
            cwd=ROOT, text=True))
        if pr['state'] != 'MERGED' or not pr['mergeCommit']:
            raise RuntimeError(f'HELD: parent PR #{number} is unmerged')
        oid = pr['mergeCommit']['oid']
        subprocess.run(['git', 'merge-base', '--is-ancestor', oid, main_head],
                       cwd=ROOT, check=True)
        parents[str(number)] = oid
    for relative, expected in frozen['merged_parent_sha256'].items():
        actual = hashlib.sha256(subprocess.check_output(
            ['git', 'show', f'origin/main:{relative}'], cwd=ROOT)).hexdigest()
        if actual != expected:
            raise RuntimeError(f'HELD: merged source drifted: {relative}')
    return {'release_main_head': main_head, 'merged_parents': parents}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    gate = release_gate(frozen)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    attempts = []
    for beta in frozen['beta_order']:
        remaining = frozen['semantic_wall_seconds'] - (time.monotonic() - started)
        if remaining <= 0:
            attempts.append({'beta': beta, 'stop_reason': 'TOTAL_WALL_CAP', 'exit_code': None})
            break
        branch = out / f'beta-{beta}'
        cmd = [sys.executable, str(HERE / 'semantic.py'), '--internal-supervised',
               '--beta', str(beta), '--out', str(branch), '--max-seconds', str(remaining)]
        receipt = run_child(cmd, cwd=ROOT,
                            stdout=out / f'beta-{beta}.stdout.txt',
                            stderr=out / f'beta-{beta}.stderr.txt',
                            wall_cap=remaining,
                            rss_cap=frozen['semantic_rss_bytes'],
                            watched_file=branch / f'beta-{beta}.cnf',
                            watched_file_cap=frozen['cnf_byte_max'],
                            file_cap=frozen['cnf_byte_max'])
        result_file = branch / 'result.json'
        result = None
        if result_file.is_file():
            try:
                result = json.loads(result_file.read_text())
            except ValueError as exc:
                receipt['result_parse_error'] = str(exc)
        receipt.update({'beta': beta, 'result_sha256': sha(result_file),
                        'result_status': result.get('status') if result else None})
        attempts.append(receipt)
        if (receipt['exit_code'] != 0 or receipt['stop_reason'] is not None or
                result is None or result.get('status') != 'SEMANTIC_PASS'):
            break
    manifest = [{'path': path.relative_to(out).as_posix(), 'bytes': path.stat().st_size,
                 'sha256': sha(path)} for path in sorted(out.rglob('*')) if path.is_file()]
    (out / 'MANIFEST.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    decision = ('SEMANTIC_PASS' if len(attempts) == len(frozen['beta_order']) and
                all(a.get('result_status') == 'SEMANTIC_PASS' and a['exit_code'] == 0
                    and a['stop_reason'] is None for a in attempts)
                else 'CENSORED_OR_FAILED')
    receipt = {'decision': decision, 'release_gate': gate, 'attempts': attempts,
               'freeze_sha256': sha(HERE / 'FROZEN.json'),
               'manifest_sha256': sha(out / 'MANIFEST.json'),
               'total_wall_seconds': time.monotonic() - started}
    (out / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'decision': decision, 'attempted_betas': [a['beta'] for a in attempts]},
                     sort_keys=True))
    return 0 if decision == 'SEMANTIC_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
