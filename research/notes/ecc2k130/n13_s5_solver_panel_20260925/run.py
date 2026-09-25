#!/usr/bin/env python3
"""Fail-closed release gate and externally bounded n13 32-label panel."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / 'symbolic_dag_dimacs_gate_20260925'))
from bounded import run_child  # noqa: E402


def sha(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def host() -> dict:
    return {'platform': platform.platform(), 'machine': platform.machine(),
            'python': platform.python_version(),
            'cc_path': shutil.which('cc'),
            'cc_version': subprocess.check_output(['cc', '--version'], text=True).splitlines()[:2],
            'cadical_path': shutil.which('cadical')}


def release_gate(frozen: dict) -> dict:
    if frozen['release_main_head'] is None:
        raise RuntimeError('HELD: #802/#804/#807/#811/#813 must merge, then rebase/re-freeze and review')
    subprocess.run([sys.executable, str(HERE / 'ci_replay.py')], cwd=ROOT, check=True)
    if host() != frozen['host_toolchain']:
        raise RuntimeError('HELD: exact frozen Mac/toolchain is unavailable')
    if sha(Path(frozen['cadical_path'])) != frozen['cadical_sha256']:
        raise RuntimeError('HELD: exact CaDiCaL binary SHA-256 drift')
    if sha(ROOT / frozen['checker_source_path']) != frozen['checker_source_sha256']:
        raise RuntimeError('HELD: checker source SHA-256 drift')
    # #804's watchdog needs working process-group RSS sampling in addition to
    # RLIMIT_AS. A host that cannot enumerate process RSS is not admitted.
    sample = subprocess.run(['ps', '-axo', 'pgid=,rss='], capture_output=True, text=True)
    if sample.returncode != 0 or not sample.stdout.strip():
        raise RuntimeError('HELD: process-tree RSS sampler unavailable')
    subprocess.run(['git', 'fetch', 'origin', 'main'], cwd=ROOT, check=True)
    main = subprocess.check_output(['git', 'rev-parse', 'origin/main'],
                                   cwd=ROOT, text=True).strip()
    if main != frozen['release_main_head']:
        raise RuntimeError('HELD: main moved after exact-head release freeze')
    subprocess.run(['git', 'merge-base', '--is-ancestor', main, 'HEAD'], cwd=ROOT, check=True)
    parents = {}
    for number in frozen['required_merged_prs']:
        pr = json.loads(subprocess.check_output(
            ['gh', 'pr', 'view', str(number), '--json', 'state,mergeCommit'],
            cwd=ROOT, text=True))
        if pr['state'] != 'MERGED' or not pr['mergeCommit']:
            raise RuntimeError(f'HELD: parent PR #{number} is unmerged')
        oid = pr['mergeCommit']['oid']
        subprocess.run(['git', 'merge-base', '--is-ancestor', oid, main],
                       cwd=ROOT, check=True)
        parents[str(number)] = oid
    for relative, expected in frozen['merged_parent_sha256'].items():
        actual = hashlib.sha256(subprocess.check_output(
            ['git', 'show', f'origin/main:{relative}'], cwd=ROOT)).hexdigest()
        if actual != expected:
            raise RuntimeError(f'HELD: merged parent source bytes drifted: {relative}')
    return {'release_main_head': main, 'merged_parents': parents,
            'host_toolchain': host(),
            'cadical_binary_sha256': frozen['cadical_sha256'],
            'checker_source_sha256': frozen['checker_source_sha256']}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    gate = release_gate(frozen)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    utc = datetime.now(timezone.utc).isoformat()
    child = out / 'panel'
    cmd = [sys.executable, str(HERE / 'panel.py'), '--internal-supervised',
           '--out', str(child)]
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    raw = run_child(cmd, cwd=ROOT, stdout=out / 'panel.stdout.txt',
                    stderr=out / 'panel.stderr.txt',
                    wall_cap=frozen['total_panel_wall_seconds'],
                    rss_cap=frozen['total_panel_rss_bytes'])
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    raw['children_plus_sampler_user_cpu_seconds'] = after.ru_utime - before.ru_utime
    raw['children_plus_sampler_system_cpu_seconds'] = after.ru_stime - before.ru_stime
    raw['utc_end'] = datetime.now(timezone.utc).isoformat()
    result_file = child / 'result.json'
    result = None
    if result_file.is_file():
        try:
            result = json.loads(result_file.read_text())
        except ValueError as exc:
            raw['result_parse_error'] = str(exc)
    raw['result_sha256'] = sha(result_file)
    raw['result_decision'] = result.get('decision') if result else None
    manifest = []
    for file in sorted(out.rglob('*')):
        if file.is_file():
            manifest.append({'path': file.relative_to(out).as_posix(),
                             'bytes': file.stat().st_size, 'sha256': sha(file)})
    (out / 'MANIFEST.json').write_text(json.dumps(manifest, sort_keys=True, indent=2) + '\n')
    total_wall = time.monotonic() - started
    decision = ('PANEL_PASS' if raw['exit_code'] == 0 and raw['stop_reason'] is None
                and result is not None and result.get('decision') == 'PANEL_PASS'
                and total_wall <= frozen['total_panel_wall_seconds']
                else 'CENSORED_OR_FAILED')
    receipt = {'decision': decision, 'domain': frozen['domain'],
               'release_gate': gate, 'frozen_sha256': sha(HERE / 'FROZEN.json'),
               'query_manifest_sha256': sha(HERE / 'QUERY_MANIFEST.json'),
               'child': raw, 'manifest_sha256': sha(out / 'MANIFEST.json'),
               'utc_start': utc, 'utc_end': datetime.now(timezone.utc).isoformat(),
               'total_wall_seconds_including_manifest': total_wall}
    (out / 'receipt.json').write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'decision': decision, 'reported_sat_lifted': result.get('sat_lifted') if result else None,
                      'reported_proved_unsat': result.get('proved_unsat') if result else None},
                     sort_keys=True))
    return 0 if decision == 'PANEL_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
