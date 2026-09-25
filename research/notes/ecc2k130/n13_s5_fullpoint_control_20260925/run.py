#!/usr/bin/env python3
"""Fail-closed, bounded n13 S5 semantic prerequisite supervisor."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from bounded import run_child, sha  # noqa: E402


def gate(frozen: dict) -> dict:
    if frozen['release_main_head'] is None:
        raise RuntimeError('HELD: #802/#804/#807/#811 must merge, then rebase/re-freeze')
    subprocess.run([sys.executable, str(HERE / 'ci_replay.py')], cwd=ROOT, check=True)
    subprocess.run(['git', 'fetch', 'origin', 'main'], cwd=ROOT, check=True)
    main = subprocess.check_output(['git', 'rev-parse', 'origin/main'],
                                   cwd=ROOT, text=True).strip()
    if main != frozen['release_main_head']:
        raise RuntimeError('HELD: main moved after release freeze')
    subprocess.run(['git', 'merge-base', '--is-ancestor', main, 'HEAD'], cwd=ROOT, check=True)
    parents = {}
    for number in frozen['required_merged_prs']:
        pr = json.loads(subprocess.check_output(
            ['gh', 'pr', 'view', str(number), '--json', 'state,mergeCommit'],
            cwd=ROOT, text=True))
        if pr['state'] != 'MERGED' or not pr['mergeCommit']:
            raise RuntimeError(f'HELD: parent PR #{number} unmerged')
        oid = pr['mergeCommit']['oid']
        subprocess.run(['git', 'merge-base', '--is-ancestor', oid, main],
                       cwd=ROOT, check=True)
        parents[str(number)] = oid
    for relative, expected in frozen['merged_parent_sha256'].items():
        actual = hashlib.sha256(subprocess.check_output(
            ['git', 'show', f'origin/main:{relative}'], cwd=ROOT)).hexdigest()
        if actual != expected:
            raise RuntimeError(f'HELD: parent byte drift: {relative}')
    return {'release_main_head': main, 'parents': parents}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    release = gate(frozen)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    child = out / 'semantic'
    cmd = [sys.executable, str(HERE / 'semantic13.py'), '--internal-supervised',
           '--out', str(child)]
    receipt = run_child(cmd, cwd=ROOT, stdout=out / 'semantic.stdout.txt',
                        stderr=out / 'semantic.stderr.txt',
                        wall_cap=frozen['semantic_wall_seconds'],
                        rss_cap=frozen['semantic_rss_bytes'],
                        watched_file=child / 'base.cnf',
                        watched_file_cap=frozen['cnf_byte_max'],
                        file_cap=frozen['cnf_byte_max'])
    result_file = child / 'result.json'
    result = None
    if result_file.is_file():
        try:
            result = json.loads(result_file.read_text())
        except ValueError as exc:
            receipt['result_parse_error'] = str(exc)
    receipt['result_sha256'] = sha(result_file)
    receipt['result_status'] = result.get('status') if result else None
    files = [{'path': p.relative_to(out).as_posix(), 'bytes': p.stat().st_size,
              'sha256': sha(p)} for p in sorted(out.rglob('*')) if p.is_file()]
    (out / 'MANIFEST.json').write_text(json.dumps(files, sort_keys=True, indent=2) + '\n')
    decision = ('SEMANTIC_PASS' if receipt['exit_code'] == 0 and
                receipt['stop_reason'] is None and
                receipt['result_status'] == 'SEMANTIC_PASS'
                else 'CENSORED_OR_FAILED')
    final = {'decision': decision, 'release_gate': release, 'child': receipt,
             'freeze_sha256': sha(HERE / 'FROZEN.json'),
             'manifest_sha256': sha(out / 'MANIFEST.json')}
    (out / 'receipt.json').write_text(json.dumps(final, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'decision': decision}, sort_keys=True))
    return 0 if decision == 'SEMANTIC_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
