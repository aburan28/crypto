#!/usr/bin/env python3
"""Independent frozen-byte preflight and archive-only replay for Linux v2."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent / 'symbolic_dag_dimacs_linux_v2_20260928'
ROOT = HERE.parents[3]
FIRST = HERE.parent / 'symbolic_dag_dimacs_gate_20260925'
PREPARATION = HERE.parent / 'symbolic_dag_dimacs_linux_attempt2_20260926'
PARENT = HERE.parent / 'symbolic_dag_fullpoint_20260925'
DOMAIN = 'k0-symbolic-dag-dimacs-linux-v2'
HEX40 = re.compile(r'[0-9a-f]{40}\Z')
WORKFLOWS = ('ecc2k130-symbolic-dag-dimacs-linux-v2-check.yml',
             'ecc2k130-symbolic-dag-dimacs-linux-v2-measure.yml')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def load(path: Path):
    return json.loads(path.read_text())


def git(*args: str) -> str:
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def check_static() -> dict:
    frozen = load(HERE / 'FROZEN.json')
    require(frozen['domain'] == DOMAIN and frozen['schema'] == 'k0-dag-linux-v2-freeze-v1',
            'v2 freeze schema/domain')
    require(frozen['status'] in ('HELD', 'RELEASED'), 'unknown release state')
    require((frozen['release_main_head'] is None and
             frozen['release_pr_number'] is None) if frozen['status'] == 'HELD' else
            (HEX40.fullmatch(frozen['release_main_head'] or '') is not None and
             isinstance(frozen['release_pr_number'], int) and
             frozen['release_pr_number'] > 0), 'release fields disagree')
    require(HEX40.fullmatch(frozen['base_main_head']) is not None and
            subprocess.run(['git', 'merge-base', '--is-ancestor',
                            frozen['base_main_head'], 'HEAD'], cwd=ROOT,
                           capture_output=True).returncode == 0,
            'v2 branch does not descend from recorded merged main')
    require(frozen['release_branch'] ==
            'codex/ecc2k130-symbolic-dag-dimacs-linux-v2-held-20260928',
            'branch identity')
    require(frozen['archive_relative'] == 'evidence/run1', 'archive namespace')
    held = frozen['held_control']
    require(held['path'] == 'evidence/held_cap_control_36530704726.json' and
            held['run_id'] == 36530704726 and
            held['artifact_id'] == 11016158422 and
            held['artifact_digest'] ==
            'sha256:1753ec368784c8a8c5021e52f8d67b049874b8ba6638fdb4990b32c8b04f4e6e' and
            held['checkout_head'] == '981fdb39c6437b2cf70b6c744b0c172a8b5ba7ac' and
            held['decision'] ==
            'PASS_CHILD_BYTE_GATE_WITH_CAPPED_PREREQUISITE_REFUSALS',
            'held Ubuntu control identity')
    held_path = HERE / held['path']
    require(sha(held_path) == held['receipt_sha256'],
            'held Ubuntu control bytes changed')
    held_receipt = load(held_path)
    require(held_receipt['schema'] == 'k0-dag-linux-v2-held-toy-cap-control-v1' and
            held_receipt['decision'] == held['decision'] and
            held_receipt['checkout_head'] == held['checkout_head'] and
            held_receipt['cap_bytes'] == frozen['caps']['toy']['rss_bytes'] and
            held_receipt['runner_image_os'] == 'ubuntu24' and
            held_receipt['runner_image_version'] == '20260920.314.1' and
            held_receipt['measured_children_started'] == 0 and
            held_receipt['capped_git_classification'] == 'PACK_MMAP_REFUSAL' and
            held_receipt['capped_gh_classification'] == 'GO_PAGE_SUMMARY_REFUSAL',
            'held Ubuntu control classification')
    rows = {row['label']: row for row in held_receipt['commands']}
    require(rows['hard_limit_observation']['exit_code'] == 0 and
            rows['capped_child_local_byte_gate']['exit_code'] == 0 and
            rows['git_cap_rev_parse']['exit_code'] == 0 and
            rows['git_cap_checkout_ancestry']['exit_code'] == 128 and
            'cannot be mapped' in
            rows['git_cap_checkout_ancestry']['stderr_tail'] and
            rows['preparation_and_v2_hash_gate']['exit_code'] == 1 and
            rows['preparation_archive_gate']['exit_code'] == 1 and
            rows['gh_parent_804']['exit_code'] ==
            rows['gh_actions_current_run']['exit_code'] == 2 and
            'failed to reserve page summary memory' in
            rows['gh_parent_804']['stderr_tail'],
            'held Ubuntu control command outcomes')
    failed = frozen['release_cap_control_failure']
    require(failed['run_id'] == 36532472069 and
            failed['artifact_id'] == 11016841420 and
            failed['artifact_digest'] ==
            'sha256:7c3578ac495398016f7694450fac8e67d315a22154a01de3646bf3cda665fcf8' and
            failed['checkout_head'] == 'c4d4a65f7209bbc3c990113bc46e8b0d7239a043' and
            failed['receipt_path'] == 'evidence/release_cap_refusal_36532472069.json' and
            failed['runs_path'] == 'evidence/release_supervisor_runs_36532472069.jsonl',
            'release harmless-control failure identity')
    for path_key, sha_key in (('receipt_path', 'receipt_sha256'),
                              ('runs_path', 'runs_sha256')):
        require(sha(HERE / failed[path_key]) == failed[sha_key],
                'release harmless-control failure raw bytes changed')
    failure = load(HERE / failed['receipt_path'])
    require(failure['schema'] == 'k0-dag-linux-v2-held-toy-cap-control-v1' and
            failure['decision'] == 'FAIL_HARMLESS_CONTROL' and
            failure['release_status'] == 'RELEASED' and
            failure['measured_children_started'] == 0 and
            failure['checkout_head'] == failed['checkout_head'] and
            failure['cap_bytes'] == frozen['caps']['toy']['rss_bytes'] and
            failure['runner_image_os'] == 'ubuntu24' and
            failure['error'] == 'RuntimeError: git_cap_base_object exited 128',
            'release harmless-control failure semantics')
    failure_rows = failure['commands']
    require([row['label'] for row in failure_rows] ==
            ['hard_limit_observation', 'capped_child_local_byte_gate',
             'git_cap_rev_parse', 'git_cap_base_object'] and
            [row['exit_code'] for row in failure_rows] == [0, 0, 0, 128] and
            'packfile ' in failure_rows[-1]['stderr_tail'] and
            'cannot be mapped' in failure_rows[-1]['stderr_tail'] and
            'Cannot allocate memory' in failure_rows[-1]['stderr_tail'],
            'release capped Git base-object failure classification')
    run_rows = [json.loads(line) for line in
                (HERE / failed['runs_path']).read_text().splitlines() if line.strip()]
    require(any(row['id'] == failed['run_id'] and
                row['head_sha'] == failed['checkout_head'] for row in run_rows),
            'release control run missing from raw Actions listing')
    require(subprocess.run(['git', 'merge-base', '--is-ancestor',
                            failed['checkout_head'], 'HEAD'], cwd=ROOT,
                           capture_output=True).returncode == 0,
            'release control failure head is not ancestor of this PR')
    require(subprocess.run(['git', 'merge-base', '--is-ancestor',
                            held['checkout_head'], 'HEAD'], cwd=ROOT,
                           capture_output=True).returncode == 0,
            'held Ubuntu control head is not ancestor of this PR')
    require(frozen['first']['pr_head'] ==
            'ca5cdef8b24c6118b5c1c9bb2b8faa8d1774533b' and
            frozen['first']['merge_commit'] ==
            '8105788c1ddd4d4680f79f5d9041d5bc6a6a7150' and
            frozen['preparation']['pr_head'] ==
            '2f0489cf638168dd8a62b719f9c60fd1ea203bee' and
            frozen['preparation']['merge_commit'] ==
            'aa2c05c9ddbf15133878aa64e2cace8ee5cef1d8',
            'merged parent identity')
    if frozen['status'] == 'RELEASED':
        for ancestor, descendant in (
                (frozen['release_main_head'], 'HEAD'),
                (frozen['first']['merge_commit'], frozen['release_main_head']),
                (frozen['preparation']['merge_commit'], frozen['release_main_head'])):
            require(subprocess.run(['git', 'merge-base', '--is-ancestor',
                                    ancestor, descendant], cwd=ROOT,
                                   capture_output=True).returncode == 0,
                    'release main or parent is outside reviewed ancestry')
    for commit in (frozen['first']['merge_commit'],
                   frozen['preparation']['merge_commit']):
        require(subprocess.run(['git', 'merge-base', '--is-ancestor', commit,
                                'HEAD'], cwd=ROOT, capture_output=True).returncode == 0,
                'merged parent not in checkout ancestry')
    first = load(FIRST / 'FROZEN.json')
    prepare = load(PREPARATION / 'PREPARE.json')
    prep_freeze = load(PREPARATION / 'LINUX_PREPARATION_FREEZE.json')
    require(sha(FIRST / 'FROZEN.json') == frozen['first']['freeze_sha256'] ==
            prepare['attempt1_freeze_sha256'], 'first freeze changed')
    require(sha(FIRST / 'evidence/receipt.json') ==
            frozen['first']['receipt_sha256'] == prepare['attempt1_receipt_sha256'],
            'first failed receipt changed')
    require(sha(FIRST / 'evidence/MANIFEST.json') ==
            frozen['first']['manifest_sha256'] == prepare['attempt1_manifest_sha256'],
            'first failed manifest changed')
    require(frozen['first']['source_sha256'] == first['source_sha256'],
            'original source ledger changed')
    for name, digest in frozen['first']['source_sha256'].items():
        require(sha(FIRST / name) == digest, f'first source drift: {name}')
    require(frozen['first']['input_sha256'] == prepare['input_sha256'] ==
            sha(FIRST / 'INPUT.json'), 'input changed')
    require(frozen['first']['semantic_source_sha256'] ==
            prepare['semantic_source_sha256'], 'semantic source ledger changed')
    require(sha(ROOT / '.github/workflows/ecc2k130-symbolic-dag-dimacs-gate.yml') ==
            frozen['first']['workflow_sha256'] == first['workflow_sha256'],
            'first workflow drift')
    require(set(frozen['preparation']['source_sha256']) ==
            {'PROTOCOL.md', 'check_parent.py', 'harmless_cap_probe.py',
             'prepare_cadical.py', 'source_manifest.py',
             'verify_preparation.py'}, 'preparation source set')
    for name, digest in frozen['preparation']['source_sha256'].items():
        require(sha(PREPARATION / name) == digest,
                f'preparation source drift: {name}')
    require(sha(ROOT / '.github/workflows/ecc2k130-symbolic-dag-dimacs-linux-preparation.yml') ==
            frozen['preparation']['workflow_sha256'],
            'preparation workflow drift')
    require(sha(PREPARATION / 'PREPARE.json') ==
            frozen['preparation']['prepare_sha256'], 'PREPARE changed')
    require(sha(PREPARATION / 'LINUX_PREPARATION_FREEZE.json') ==
            frozen['preparation']['freeze_sha256'], 'preparation freeze changed')
    require(sha(PREPARATION / 'preparation_run_1/MANIFEST.json') ==
            frozen['preparation']['manifest_sha256'], 'preparation manifest changed')
    require(prep_freeze['status'] == 'PASS_HARMLESS_PREPARATION_ONLY' and
            prep_freeze['measured_attempt2_admitted'] is False and
            prep_freeze['measured_archive'] is None and
            prep_freeze['v2_release_head'] is None,
            'preparation falsely claims measured admission')
    for key in ('cap_probe_receipt_sha256', 'build_receipt_sha256',
                'linux_binary_sha256', 'source_commit', 'source_tree',
                'source_manifest_sha256'):
        require(prep_freeze[key] == frozen['preparation'][key],
                f'preparation identity drift: {key}')
    binary_rel = frozen['linux_binary_path']
    require(binary_rel == 'symbolic_dag_dimacs_linux_attempt2_20260926/' +
            prep_freeze['linux_binary_path'], 'Linux binary path changed')
    binary = HERE.parent / binary_rel
    require(binary.is_file() and sha(binary) == frozen['linux_binary_sha256'] ==
            prep_freeze['linux_binary_sha256'] and
            binary.stat().st_size == frozen['linux_binary_bytes'] ==
            prep_freeze['linux_binary_bytes'], 'Linux binary bytes changed')
    tree_path = binary.relative_to(ROOT).as_posix()
    mode, kind, blob, path = git('ls-tree', 'HEAD', tree_path).split(maxsplit=3)
    require((mode, kind, blob, path) ==
            ('100755', 'blob', frozen['linux_binary_git_blob'], tree_path),
            'Linux binary Git blob or executable mode changed')
    require(binary.read_bytes()[:6] == b'\x7fELF\x02\x01', 'Linux binary is not ELF64')
    require(frozen['linux_binary_git_blob'] ==
            '3f50d3a4218609874b01e02bdf61ad699e6b3f71', 'binary Git identity')
    require(set(frozen['source_sha256']) ==
            {'PROTOCOL.md', 'ci_replay.py', 'run.py', 'v2_child.py',
             'selftest.py', 'held_cap_control.py'}, 'v2 source set')
    for name, digest in frozen['source_sha256'].items():
        require(sha(HERE / name) == digest, f'v2 source drift: {name}')
    require(set(frozen['workflow_sha256']) == set(WORKFLOWS), 'workflow set')
    for name, digest in frozen['workflow_sha256'].items():
        require(sha(ROOT / '.github/workflows' / name) == digest,
                f'workflow drift: {name}')
    original_caps = {
        'toy': {'wall_seconds': first['toy_external_wall_cap_seconds'],
                'rss_bytes': first['toy_rss_cap_bytes'],
                'cnf_bytes': first['toy_cnf_byte_cap']},
        'panel': {'wall_seconds': first['panel_external_wall_cap_seconds'],
                  'rss_bytes': first['panel_rss_cap_bytes'],
                  'query_rss_bytes': first['n13_rss_cap_bytes'],
                  'cnf_bytes': first['n13_cnf_byte_cap'],
                  'solver_wall_seconds': first['solver_wall_cap_seconds'],
                  'checker_wall_seconds': first['checker_wall_cap_seconds'],
                  'checker_compile_wall_seconds':
                      first['checker_compile_wall_cap_seconds'],
                  'proof_bytes': first['proof_byte_cap']},
        'n131': {'wall_seconds': first['n131_external_wall_cap_seconds'],
                 'rss_bytes': first['n131_rss_cap_bytes'],
                 'cnf_bytes': first['n131_cnf_byte_cap'],
                 'node_cap': first['n131_node_cap']},
    }
    require(frozen['caps'] == original_caps, 'inherited phase caps changed')
    require(frozen['upstream_guard_paths'] == [
        '.github/workflows/ecc2k130-symbolic-dag-dimacs-gate.yml',
        '.github/workflows/ecc2k130-symbolic-dag-dimacs-linux-preparation.yml',
        'research/notes/ecc2k130/n13_oaware_sat_benchmark_20260925/evidence/panel/Q0T3-cadical.stdout',
        'research/notes/ecc2k130/symbolic_dag_dimacs_gate_20260925/',
        'research/notes/ecc2k130/symbolic_dag_dimacs_linux_attempt2_20260926/',
        'research/notes/ecc2k130/symbolic_dag_fullpoint_20260925/',
    ], 'upstream guard scope changed')
    # This existing verifier checks every raw preparation file, all hard-cap
    # observations, both independent builds and the unchanged first failure.
    check = subprocess.run([sys.executable, str(PREPARATION / 'verify_preparation.py')],
                           cwd=ROOT, capture_output=True, text=True, timeout=120)
    require(check.returncode == 0 and
            'PASS_HARMLESS_PREPARATION_ONLY' in check.stdout,
            'independent Linux preparation replay failed: ' + check.stderr[-500:])
    return frozen



def check_child_bytes() -> dict:
    """Git-free byte gate for a child already admitted by the supervisor.

    The 512-MiB toy limit prevents Git from mapping this checkout's packfile.
    The supervisor therefore performs the full ancestry, GitHub and
    preparation replay before dispatch; this child check pins the local code
    and solver bytes that can be used after that admission.
    """
    frozen = load(HERE / 'FROZEN.json')
    require(frozen['domain'] == DOMAIN and frozen['schema'] == 'k0-dag-linux-v2-freeze-v1',
            'child freeze domain/schema')
    require(frozen['status'] in ('HELD', 'RELEASED'), 'child release state')
    for name, digest in frozen['source_sha256'].items():
        require(sha(HERE / name) == digest, f'child v2 source drift: {name}')
    for name, digest in frozen['first']['source_sha256'].items():
        require(sha(FIRST / name) == digest, f'child predecessor source drift: {name}')
    require(sha(FIRST / 'INPUT.json') == frozen['first']['input_sha256'],
            'child predecessor input drift')
    binary = HERE.parent / frozen['linux_binary_path']
    require(binary.is_file() and sha(binary) == frozen['linux_binary_sha256'] and
            binary.stat().st_size == frozen['linux_binary_bytes'] and
            binary.read_bytes()[:6] == b'\x7fELF\x02\x01',
            'child Linux solver bytes drift')
    return frozen

def _expected_panel() -> list[tuple[str, str, dict]]:
    inputs = load(FIRST / 'INPUT.json')
    cases: list[tuple[str, str, dict]] = []
    for pair in inputs['pairs']:
        for suffix, expected, result_key in (('_sat', 'SAT', 'sat_r'),
                                             ('_wrong', 'UNSAT', 'unsat_r')):
            cases.append((pair['id'] + suffix, expected,
                          {'p': tuple(pair['p']), 'q': tuple(pair['q']),
                           'r': tuple(pair[result_key])}))
    for row in inputs['invalid']:
        cases.append((row['id'], 'UNSAT',
                      {key: tuple(row[key]) for key in ('p', 'q', 'r')}))
    require(len(cases) == 20, 'frozen panel case count')
    return cases


def _decompress_limited(path: Path, cap: int) -> bytes:
    with gzip.open(path, 'rb') as stream:
        raw = stream.read(cap + 1)
        require(len(raw) <= cap and stream.read(1) == b'', 'decompressed cap exceeded')
    return raw


def _check_process(row: dict, stdout: Path, stderr: Path, wall_cap: float,
                   rss_cap: int, *, exit_codes: tuple[int, ...]) -> None:
    require(row['exit_code'] in exit_codes and row['stop_reason'] is None and
            0 <= row['wall_seconds'] <= wall_cap and
            0 <= row['sampled_peak_rss_bytes'] <= rss_cap and
            row['rss_address_space_cap_bytes'] == rss_cap,
            'process exit/stop/wall/address-space cap')
    require(sha(stdout) == row['stdout_sha256'] and
            sha(stderr) == row['stderr_sha256'], 'process stream hash')


def _replay_toy(archive: Path, result: dict, frozen: dict) -> dict:
    sys.path.insert(0, str(FIRST))
    from verify import verify_parent_rows
    cnfs = {n: archive / 'toy' / f'n{n}.cnf' for n in (2, 3)}
    require(result['decision'] == 'PASS' and
            [row['n'] for row in result['fields']] == [2, 3], 'toy result')
    for row, n in zip(result['fields'], (2, 3)):
        require(cnfs[n].stat().st_size == row['bytes'] <=
                frozen['caps']['toy']['cnf_bytes'] and
                sha(cnfs[n]) == row['sha256'], 'toy CNF bytes')
    replayed = verify_parent_rows(cnfs)
    require([row['n'] for row in replayed['fields']] == [2, 3], 'toy row replay')
    return replayed


def _replay_panel(archive: Path, result: dict, frozen: dict,
                  phase_out: str, checkout_root: str) -> int:
    sys.path.insert(0, str(FIRST))
    from export import build_relation, fixed_point_units
    from verify import (ReferenceField, check_relation_cnf, evaluate_cnf,
                        lift_solver_model, parse_cnf, parse_solver_output)
    require(result['decision'] == 'PASS', 'panel stage did not pass')
    cap = frozen['caps']['panel']
    cases = _expected_panel()
    rows = result['queries']
    require(len(rows) == 20 and [row['id'] for row in rows] ==
            [item[0] for item in cases], 'panel query order or count')
    base = archive / 'panel/base.cnf'
    require(base.stat().st_size == result['base']['bytes'] <= cap['cnf_bytes'] and
            sha(base) == result['base']['sha256'], 'panel base CNF')
    relation = build_relation(13, 0x201B)
    field = ReferenceField(13, 0x201B)
    variables, clauses_count, clauses = parse_cnf(base)
    require((variables, clauses_count) ==
            (result['base']['variables'], result['base']['clauses']),
            'panel base dimensions')
    check_relation_cnf(relation, clauses)
    compile_row = result['checker_compile']
    _check_process(compile_row, archive / 'panel/checker_compile.stdout.txt',
                   archive / 'panel/checker_compile.stderr.txt',
                   cap['checker_compile_wall_seconds'], cap['query_rss_bytes'],
                   exit_codes=(0,))
    require(compile_row['argv'] ==
            ['cc', '-O2', '-o', phase_out + '/drat-trim',
             checkout_root + '/research/notes/ecc2k130/' +
             'symbolic_dag_dimacs_gate_20260925/third_party/drat-trim.c'],
            'checker compile command')
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        checker = temp / 'drat-trim'
        subprocess.run(['cc', '-O2', '-o', str(checker),
                        str(FIRST / 'third_party/drat-trim.c')],
                       check=True, capture_output=True, timeout=60)
        for row, (name, expected, points) in zip(rows, cases):
            require(row['id'] == name and row['expected'] == expected and
                    {key: tuple(row['points'][key]) for key in ('p', 'q', 'r')} == points,
                    'panel case identity')
            packed = archive / 'panel/queries' / f'{name}.cnf.gz'
            require(sha(packed) == row['cnf_gzip_sha256'], 'query gzip hash')
            query = temp / 'query.cnf'
            query.write_bytes(_decompress_limited(packed, cap['cnf_bytes']))
            require(sha(query) == row['cnf']['sha256'] and
                    query.stat().st_size == row['cnf']['bytes'], 'query CNF hash')
            vars_count, count, clauses = parse_cnf(query)
            require((vars_count, count) ==
                    (row['cnf']['variables'], row['cnf']['clauses']),
                    'query dimensions')
            units = fixed_point_units(relation, points['p'], points['q'], points['r'])
            require(row['units'] == units, 'query point units')
            check_relation_cnf(relation, clauses, units)
            solver = row['solver']
            out_log = archive / 'panel/logs' / f'{name}.solver.stdout.txt'
            err_log = archive / 'panel/logs' / f'{name}.solver.stderr.txt'
            _check_process(solver, out_log, err_log, cap['solver_wall_seconds'],
                           cap['query_rss_bytes'], exit_codes=(10, 20))
            expected_solver = [checkout_root + '/research/notes/ecc2k130/' +
                               frozen['linux_binary_path'], '--no-binary',
                               phase_out + '/query.cnf']
            if expected == 'UNSAT':
                expected_solver.append(phase_out + '/proof.drat')
            require(solver['argv'] == expected_solver, 'solver command drift')
            observed, values = parse_solver_output(out_log.read_text(),
                                                    vars_count, solver['exit_code'])
            require(observed == expected == row['observed'] and row['error'] is None,
                    'solver answer or producer error')
            if expected == 'SAT':
                require(row['status'] == 'SAT_LIFTED' and values is not None and
                        evaluate_cnf(clauses, values), 'SAT model not clause-valid')
                lifted = lift_solver_model(relation, values, points['p'],
                                           points['q'], points['r'], field)
                require(row['lifted']['lambda'] == lifted['lambda'] and
                        {key: tuple(row['lifted'][key]) for key in ('p', 'q', 'r')} ==
                        {key: lifted[key] for key in ('p', 'q', 'r')},
                        'SAT model did not lift to exact full point')
                require(row['checker'] is None and row['proof_gzip_sha256'] is None,
                        'unexpected SAT proof/checker')
            else:
                require(row['status'] == 'PROVED_UNSAT' and
                        row['proof_text_valid'] is True and row['lifted'] is None,
                        'UNSAT admission')
                packed_proof = archive / 'panel/proofs' / f'{name}.drat.gz'
                require(sha(packed_proof) == row['proof_gzip_sha256'],
                        'proof gzip hash')
                proof = temp / 'proof.drat'
                proof.write_bytes(_decompress_limited(packed_proof,
                                                      cap['proof_bytes']))
                proof_bytes = proof.read_bytes()
                require(bool(proof_bytes) and b'\x00' not in proof_bytes and
                        all(byte < 128 for byte in proof_bytes), 'nontext DRAT')
                checked = row['checker']
                require(checked is not None, 'missing archived checker receipt')
                _check_process(checked,
                               archive / 'panel/logs' / f'{name}.checker.stdout.txt',
                               archive / 'panel/logs' / f'{name}.checker.stderr.txt',
                               cap['checker_wall_seconds'], cap['query_rss_bytes'],
                               exit_codes=(0,))
                require(checked['argv'] == [phase_out + '/drat-trim',
                                            phase_out + '/query.cnf',
                                            phase_out + '/proof.drat'],
                        'checker command drift')
                fresh = subprocess.run([str(checker), str(query), str(proof)],
                                       capture_output=True,
                                       timeout=cap['checker_wall_seconds'])
                require(fresh.returncode == 0, 'fresh DRAT-trim did not verify')
    return len(rows)


def _replay_n131(archive: Path, result: dict, frozen: dict) -> str:
    sys.path.insert(0, str(FIRST))
    from export import build_relation, expected_clauses
    require(result['decision'] == 'PASS' and result['single_edge_only'] is True and
            result['n'] == 131 and result['modulus'] ==
            load(FIRST / 'INPUT.json')['n131_single_edge_modulus'],
            'n131 single-edge scope or field')
    cap = frozen['caps']['n131']
    packed = archive / 'n131/single_edge.cnf.gz'
    require(sha(packed) == result['cnf_gzip_sha256'] and
            packed.stat().st_size == result['cnf_gzip_bytes'], 'n131 gzip hash')
    with tempfile.TemporaryDirectory() as directory:
        cnf = Path(directory) / 'single_edge.cnf'
        cnf.write_bytes(_decompress_limited(packed, cap['cnf_bytes']))
        require(sha(cnf) == result['cnf']['sha256'] and
                cnf.stat().st_size == result['cnf']['bytes'] and
                result['cnf']['dag']['variables'] == 920 and
                result['cnf']['dag']['model_limbs'] == 15 and
                result['cnf']['dag']['variables'] + 2 <=
                result['cnf']['dag']['total_nodes'] ==
                result['cnf']['variables'] <= cap['node_cap'],
                'n131 CNF/width/node metadata')
        relation = build_relation(131, result['modulus'])
        counts = relation.dag.counts()
        require(counts == result['cnf']['dag'] and
                counts['total_nodes'] <= cap['node_cap'] and
                result['cnf']['variables'] == counts['total_nodes'] and
                result['cnf']['clauses'] == expected_clauses(relation.dag, []) and
                result['cnf']['relation_output_literal'] == relation.output + 1,
                'n131 reconstructed relation metadata/output')
        _check_relation_stream(relation, cnf, result['cnf']['clauses'])
    return result['cnf']['sha256']


def _check_relation_stream(relation, cnf: Path, expected_count: int) -> None:
    """Check canonical clauses against the frozen DAG without retaining a large CNF.

    Unlike a dimension-only parse, this binds every gate's inputs, operation,
    output ID and clause order to the reconstructed relation, then binds the
    asserted final unit to that relation's actual output node.
    """
    nodes = relation.dag.nodes
    with cnf.open('rt', encoding='ascii') as stream:
        require(stream.readline() == f'p cnf {len(nodes)} {expected_count}\n',
                'relation DIMACS header differs from reconstructed DAG')
        seen = 0

        def clause() -> tuple[int, ...]:
            nonlocal seen
            line = stream.readline()
            require(bool(line), 'relation DIMACS ended early')
            parts = line.split()
            require(len(parts) >= 2 and parts[-1] == '0',
                    'malformed relation clause')
            try:
                values = tuple(int(part) for part in parts[:-1])
            except ValueError as exc:
                raise AssertionError('noninteger relation literal') from exc
            require(all(lit != 0 and abs(lit) <= len(nodes) for lit in values),
                    'relation literal outside reconstructed DAG')
            require(line == ' '.join(map(str, values)) + ' 0\n',
                    'relation DIMACS clause is not canonical')
            seen += 1
            return values

        require(clause() == (-1,) and clause() == (2,),
                'relation constant clauses changed')
        for node_id, (op, a, b) in enumerate(nodes[2:], 2):
            if op == 'var':
                continue
            x, y, z = a + 1, b + 1, node_id + 1
            if op == 'xor':
                expected = ((-x, -y, -z), (x, y, -z),
                            (x, -y, z), (-x, y, z))
            elif op == 'and':
                expected = ((-x, -y, z), (x, -z), (y, -z))
            else:
                raise AssertionError('unknown reconstructed DAG gate')
            require(tuple(clause() for _ in expected) == expected,
                    f'relation {op} clauses changed at node {node_id}')
        require(clause() == (relation.output + 1,),
                'relation output assertion changed')
        require(seen == expected_count and stream.read(1) == '',
                'relation DIMACS has missing or extra clauses')


def check_archive(receipt_path: Path, frozen: dict) -> dict:
    receipt_path = receipt_path.resolve()
    archive = receipt_path.parent
    receipt = load(receipt_path)
    require(receipt['domain'] == DOMAIN and
            receipt['freeze_sha256'] == sha(HERE / 'FROZEN.json'),
            'outcome freeze identity')
    require(receipt['decision'] in ('PASS', 'FAIL_OR_CENSORED') and
            isinstance(receipt['runner_total_wall_seconds'], (int, float)) and
            receipt['runner_total_wall_seconds'] >= 0,
            'outcome decision/cost vocabulary')
    manifest_path = archive / 'MANIFEST.json'
    require(sha(manifest_path) == receipt['manifest_sha256'], 'manifest digest')
    manifest = load(manifest_path)
    actual = {p.relative_to(archive).as_posix(): p for p in archive.rglob('*')
              if p.is_file() and p not in (receipt_path, manifest_path)}
    require([row['path'] for row in manifest] == sorted(actual, key=lambda rel: Path(rel).parts),
            'manifest omits or invents an archive file')
    for row in manifest:
        rel = Path(row['path'])
        require(not rel.is_absolute() and '..' not in rel.parts and
                not actual[row['path']].is_symlink(), 'unsafe archive member')
        path = actual[row['path']]
        require(path.stat().st_size == row['bytes'] and
                sha(path) == row['sha256'], f'archive drift: {row["path"]}')
    gate = receipt['release_gate']
    attempts = receipt['attempts']
    require([row['phase'] for row in attempts] ==
            ['toy', 'panel', 'n131'][:len(attempts)] and len(attempts) <= 3,
            'phase order/count')
    require(receipt['runner_total_wall_seconds'] >=
            sum(row['wall_seconds'] for row in attempts),
            'runner wall does not include child walls')
    if gate is None:
        require(not attempts and receipt['decision'] == 'FAIL_OR_CENSORED' and
                receipt['error'] and (archive / 'PRE_DISPATCH_REFUSAL.json').is_file(),
                'pre-dispatch refusal was mislabeled')
        return {'decision': 'ARCHIVED_PRE_DISPATCH_REFUSAL', 'phases': 0}
    require(gate == load(archive / 'DISPATCH.json'), 'dispatch receipt differs')
    target_probe = load(archive / 'target_rlimit_probe.json')
    probe_gate = gate['target_cap_probe']
    require(sha(archive / 'target_rlimit_probe.json') ==
            probe_gate['receipt_sha256'] and
            sha(archive / 'target_rlimit_probe.stdout.txt') ==
            probe_gate['stdout_sha256'] and
            sha(archive / 'target_rlimit_probe.stderr.txt') ==
            probe_gate['stderr_sha256'], 'fresh cap probe archive hash')
    require(target_probe['schema'] == 'k0-dag-dimacs-linux-rlimit-probe-v1' and
            target_probe['decision'] == 'PASS_HARMLESS_CAPS_ONLY' and
            target_probe['checkout_head'] == gate['reviewed_head'] and
            target_probe['runner_image_os'] == probe_gate['runner_image_os'] ==
            'ubuntu24' and
            target_probe['runner_image_version'] ==
            probe_gate['runner_image_version'], 'fresh cap host provenance')
    require([row['cap_bytes'] for row in target_probe['caps']] ==
            probe_gate['caps_bytes'] == [frozen['caps']['toy']['rss_bytes'],
                                       frozen['caps']['panel']['query_rss_bytes'],
                                       frozen['caps']['panel']['rss_bytes']],
            'fresh cap schedule')
    for cap, row in zip(probe_gate['caps_bytes'], target_probe['caps']):
        require(row['true_exit_code'] == row['python_exit_code'] == 0 and
                row['observation']['limit'] == [cap, cap] and
                row['observation']['over_cap_allocation'] == 'REJECTED' and
                row['observation']['exception'] == 'OSError',
                'fresh hard address-space control')
    if receipt['error'] is not None:
        require((archive / 'RUNNER_ERROR.json').is_file(),
                'missing runner error receipt')
    require(gate['release_main_head'] == frozen['release_main_head'] and
            gate['freeze_sha256'] == sha(HERE / 'FROZEN.json') and
            gate['preparation_freeze_sha256'] ==
            frozen['preparation']['freeze_sha256'] and
            gate['linux_binary_sha256'] == frozen['linux_binary_sha256'] and
            gate['reviewed_head'] == gate['checkout_head'] == gate['pr_head'] and
            gate['upstream_changed_paths'] == gate['upstream_guarded_history'] == [],
            'release gate identity')
    for key in ('reviewed_head', 'dispatch_main_head', 'checkout_head'):
        require(HEX40.fullmatch(gate[key]) is not None, 'invalid gate SHA')
    require(gate['one_shot']['matching_labeled_runs'] ==
            [gate['event']['run_id']] and
            gate['one_shot']['run_id'] == gate['event']['run_id'] and
            gate['event']['pr_number'] == gate['pr_number'] ==
            frozen['release_pr_number'], 'one-shot archived gate')
    for ancestor, descendant in ((gate['release_main_head'],
                                  gate['dispatch_main_head']),
                                 (gate['checkout_head'], 'HEAD')):
        require(subprocess.run(['git', 'merge-base', '--is-ancestor', ancestor,
                                descendant], cwd=ROOT, capture_output=True).returncode == 0,
                'archive lineage no longer replays')
    changed = git('diff', '--name-only', gate['release_main_head'] + '..' +
                  gate['dispatch_main_head'], '--', *frozen['upstream_guard_paths'])
    history = git('log', '--full-history', '--format=%H',
                  gate['release_main_head'] + '..' + gate['dispatch_main_head'],
                  '--', *frozen['upstream_guard_paths'])
    require(not changed and not history, 'guarded main changed in archived interval')
    all_pass = bool(attempts) and len(attempts) == 3
    successful = {}
    for row in attempts:
        phase = row['phase']
        cap = frozen['caps'][phase]
        stdout = archive / f'{phase}.stdout.txt'
        stderr = archive / f'{phase}.stderr.txt'
        require(sha(stdout) == row['stdout_sha256'] and
                sha(stderr) == row['stderr_sha256'], 'phase stream hash')
        require(row['cwd'] == gate['checkout_root'] and
                row['rss_address_space_cap_bytes'] == cap['rss_bytes'] and
                row['wall_seconds'] >= 0 and
                row['sampled_peak_rss_bytes'] >= 0,
                'phase command context/cap')
        argv = row['argv']
        require(len(argv) == 12 and argv[1] == gate['checkout_root'] +
                '/research/notes/ecc2k130/symbolic_dag_dimacs_linux_v2_20260928/v2_child.py' and
                argv[2:4] == ['--phase', phase] and argv[4] == '--out' and
                argv[6:8] == ['--expected-head', gate['reviewed_head']] and
                argv[8:10] == ['--dispatch', str(Path(argv[5]).parent / 'DISPATCH.json')] and
                argv[10:] == ['--dispatch-sha256', sha(archive / 'DISPATCH.json')],
                'phase child command')
        result_path = archive / phase / 'result.json'
        result = None
        if result_path.is_file():
            require(sha(result_path) == row['result_sha256'], 'phase result hash')
            try:
                result = load(result_path)
            except ValueError:
                require('result_parse_error' in row,
                        'unreported malformed phase result')
        else:
            require(row['result_sha256'] is None, 'missing phase result bytes')
        require(row['reported_decision'] ==
                (result['decision'] if result else None) and
                row['reported_cpu_seconds'] ==
                (result['cpu_seconds'] if result else None) and
                row['reported_peak_rss_bytes'] ==
                (result['peak_rss_bytes'] if result else None),
                'phase result fields')
        if result:
            require(result['phase'] == phase and result['v2_wrapper'] == DOMAIN and
                    result['domain'] == 'k0-symbolic-dag-dimacs-gate-v1' and
                    result['cpu_seconds'] >= 0 and result['wall_seconds'] >= 0,
                    'phase result provenance')
        passed = (row['exit_code'] == 0 and row['stop_reason'] is None and
                  row['wall_seconds'] <= cap['wall_seconds'] and
                  row['sampled_peak_rss_bytes'] <= cap['rss_bytes'] and
                  result is not None and result['decision'] == 'PASS' and
                  result['peak_rss_bytes'] <= cap['rss_bytes'])
        all_pass = all_pass and passed
        if passed:
            if phase == 'toy':
                successful['toy'] = _replay_toy(archive, result, frozen)
            elif phase == 'panel':
                successful['panel_queries'] = _replay_panel(
                    archive, result, frozen, argv[5], gate['checkout_root'])
            else:
                successful['n131_cnf_sha256'] = _replay_n131(archive, result,
                                                              frozen)
        else:
            require(row is attempts[-1], 'runner continued after a failed phase')
    require((receipt['decision'] == 'PASS') == all_pass and
            (receipt['error'] is None or receipt['decision'] == 'FAIL_OR_CENSORED'),
            'outcome decision not recomputed from all phases')
    if receipt['decision'] == 'PASS':
        return {'decision': 'PASS_INDEPENDENT_REPLAY', **successful}
    return {'decision': 'ARCHIVED_FAILURE_ONLY',
            'phases': len(attempts), 'successful_prefix': list(successful)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--evidence', type=Path)
    args = parser.parse_args()
    frozen = check_static()
    if args.evidence:
        evidence = check_archive(args.evidence, frozen)
    else:
        committed = HERE / frozen['archive_relative'] / 'receipt.json'
        evidence = check_archive(committed, frozen) if committed.is_file() else None
        # RELEASED must pass hash-only CI before the unique label starts a
        # measured child, so an absent archive is valid before that event.
    print(json.dumps({'decision': 'PASS_HASH_AND_ARCHIVE_REPLAY',
                      'release_status': frozen['status'],
                      'freeze_sha256': sha(HERE / 'FROZEN.json'),
                      'evidence': evidence}, sort_keys=True))


if __name__ == '__main__':
    main()
