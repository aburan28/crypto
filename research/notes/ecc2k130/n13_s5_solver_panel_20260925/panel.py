#!/usr/bin/env python3
"""Held cold 32-label S5 SAT/DRAT panel. Only run through release supervisor."""
from __future__ import annotations

import argparse
import fcntl
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


def supervised(command: list[str], *, stdout: Path, stderr: Path,
               wall: float, rss: int, watched: Path | None = None,
               watched_cap: int | None = None, file_cap: int | None = None) -> dict:
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt = run_child(command, cwd=ROOT, stdout=stdout, stderr=stderr,
                        wall_cap=wall, rss_cap=rss,
                        watched_file=watched, watched_file_cap=watched_cap,
                        file_cap=file_cap)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    # Includes watchdog ps samplers. It is a raw platform diagnostic, not an
    # exclusive solver CPU count or a common operation unit against rho.
    receipt['children_plus_sampler_user_cpu_seconds'] = after.ru_utime - before.ru_utime
    receipt['children_plus_sampler_system_cpu_seconds'] = after.ru_stime - before.ru_stime
    receipt['utc_end'] = datetime.now(timezone.utc).isoformat()
    return receipt


def read_json(path: Path) -> dict | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except ValueError:
        return None


def progress(out: Path, frozen: dict, rows: list[dict], compile_row: dict,
             *, status: str) -> None:
    document = {'domain': frozen['domain'], 'status': status,
                'compile': compile_row, 'cases': rows,
                'completed': len(rows),
                'sat_lifted': sum(r['status'] == 'SAT_LIFTED' for r in rows),
                'proved_unsat': sum(r['status'] == 'PROVED_UNSAT' for r in rows)}
    temp = out / 'progress.json.tmp'
    temp.write_text(json.dumps(document, sort_keys=True, indent=2) + '\n')
    temp.replace(out / 'progress.json')


def run(out: Path, frozen: dict) -> dict:
    if frozen['release_main_head'] is None:
        raise RuntimeError('HELD: no merged-parent release')
    if host() != frozen['host_toolchain']:
        raise RuntimeError('HELD: host/toolchain differs from exact frozen Mac build')
    solver = Path(frozen['cadical_path'])
    source = ROOT / frozen['checker_source_path']
    fixed_checker = Path(frozen['checker_build_path'])
    if sha(solver) != frozen['cadical_sha256'] or sha(source) != frozen['checker_source_sha256']:
        raise RuntimeError('HELD: solver binary or checker source SHA-256 drift')
    manifest = json.loads((HERE / 'QUERY_MANIFEST.json').read_text())
    labels = manifest['labels']
    if len(labels) != 32 or sum(row['expected_status'] == 'SAT' for row in labels) != 5:
        raise AssertionError('exact 5/27 query manifest drift')
    # This exact fixed output path is part of the freeze. Differently named
    # Mach-O checker outputs hash differently on the observed Apple toolchain.
    fixed_checker.parent.mkdir(parents=True, exist_ok=True)
    fixed_checker.unlink(missing_ok=True)
    compile_command = [frozen['cc_path'], '-O2', '-o', str(fixed_checker), str(source)]
    compiled = supervised(compile_command,
                          stdout=out / 'checker_compile.stdout.txt',
                          stderr=out / 'checker_compile.stderr.txt',
                          wall=frozen['checker_compile_wall_seconds'],
                          rss=frozen['checker_compile_rss_bytes'])
    compiled['source_sha256'] = sha(source)
    compiled['binary_sha256'] = sha(fixed_checker)
    compiled['toolchain'] = host()
    if (compiled['exit_code'] != 0 or compiled['stop_reason'] is not None or
            compiled['binary_sha256'] != frozen['checker_binary_sha256']):
        progress(out, frozen, [], compiled, status='HELD_CHECKER_BUILD')
        raise RuntimeError('HELD: exact checker build failed or binary hash drifted')
    archived_checker = out / 'drat-trim.exact-binary'
    shutil.copyfile(fixed_checker, archived_checker)
    archived_checker.chmod(0o700)
    if sha(archived_checker) != frozen['checker_binary_sha256']:
        raise AssertionError('archived checker binary drift')
    compiled['archived_binary_sha256'] = sha(archived_checker)
    compiled['archived_binary_bytes'] = archived_checker.stat().st_size
    rows = []
    progress(out, frozen, rows, compiled, status='RUNNING')
    for expected in labels:
        label = expected['label']
        case_dir = out / f"{expected['index']:02d}-{label}"
        case_dir.mkdir(parents=True, exist_ok=False)
        case = {'label': label, 'index': expected['index'],
                'target': expected['target'],
                'archived_exact_tuple_count': expected['archived_exact_tuple_count'],
                'expected_status': expected['expected_status'],
                'query_sha256': expected['sha256'], 'status': 'CENSORED'}
        export_dir = case_dir / 'export'
        query = export_dir / 'query.cnf'
        export_command = [sys.executable, str(HERE / 'export_one.py'),
                          '--internal-supervised', '--label', label,
                          '--out', str(export_dir)]
        export = supervised(export_command,
                            stdout=case_dir / 'export.stdout.txt',
                            stderr=case_dir / 'export.stderr.txt',
                            wall=frozen['export_wall_seconds'],
                            rss=frozen['export_rss_bytes'],
                            watched=query, watched_cap=frozen['query_cnf_byte_max'],
                            file_cap=frozen['query_cnf_byte_max'])
        case['export'] = export
        case['query_actual_sha256'] = sha(query)
        export_meta = read_json(export_dir / 'query_meta.json')
        if (export['exit_code'] != 0 or export['stop_reason'] is not None or
                case['query_actual_sha256'] != expected['sha256'] or
                not export_meta or export_meta.get('status') != 'EXACT_QUERY'):
            case['status'] = 'EXPORT_FAILURE_OR_CAP'
            rows.append(case)
            progress(out, frozen, rows, compiled, status='CENSORED_OR_FAILED')
            break  # A wrong query invalidates the remaining frozen panel.
        proof = case_dir / 'proof.drat'
        solver_out = case_dir / 'solver.stdout.txt'
        solver_command = [str(solver), '--no-binary', str(query), str(proof)]
        solver_receipt = supervised(solver_command,
                                    stdout=solver_out,
                                    stderr=case_dir / 'solver.stderr.txt',
                                    wall=frozen['solver_wall_seconds'],
                                    rss=frozen['solver_rss_bytes'],
                                    watched=proof, watched_cap=frozen['proof_byte_max'],
                                    file_cap=frozen['proof_byte_max'])
        case['solver'] = solver_receipt
        case['solver_stdout_sha256'] = sha(solver_out)
        case['solver_stderr_sha256'] = sha(case_dir / 'solver.stderr.txt')
        case['proof_sha256'] = sha(proof)
        case['proof_bytes'] = proof.stat().st_size if proof.is_file() else None
        code = solver_receipt['exit_code']
        kind = 'SAT' if code == 10 else 'UNSAT' if code == 20 else None
        if solver_receipt['stop_reason'] is None and kind is not None:
            admit_dir = case_dir / 'admission'
            admit_command = [sys.executable, str(HERE / 'admit_one.py'),
                             '--internal-supervised', '--label', label,
                             '--kind', kind, '--query', str(query),
                             '--solver-stdout', str(solver_out),
                             '--solver-exit', str(code), '--out', str(admit_dir)]
            if kind == 'UNSAT':
                admit_command += ['--proof', str(proof), '--checker', str(archived_checker)]
            admit_wall = (frozen['sat_lift_wall_seconds'] if kind == 'SAT'
                          else frozen['unsat_admission_wall_seconds'])
            admit_rss = (frozen['sat_lift_rss_bytes'] if kind == 'SAT'
                         else frozen['checker_rss_bytes'])
            admitted = supervised(admit_command,
                                  stdout=case_dir / 'admission.stdout.txt',
                                  stderr=case_dir / 'admission.stderr.txt',
                                  wall=admit_wall, rss=admit_rss)
            case['admission'] = admitted
            case['admission_result_sha256'] = sha(admit_dir / 'admission.json')
            admission_row = read_json(admit_dir / 'admission.json')
            needed = 'SAT_LIFTED' if kind == 'SAT' else 'PROVED_UNSAT'
            if (admitted['exit_code'] == 0 and admitted['stop_reason'] is None and
                    admission_row and admission_row.get('status') == needed):
                case['status'] = needed if kind == expected['expected_status'] else 'EXACT_COUNT_MISMATCH'
                case['verified_witness'] = admission_row.get('verified') if kind == 'SAT' else None
            else:
                case['status'] = 'ADMISSION_FAILURE_OR_CAP'
        else:
            case['status'] = 'SOLVER_UNKNOWN_OR_CAP'
        rows.append(case)
        progress(out, frozen, rows, compiled, status='RUNNING')
    sat = sum(row['status'] == 'SAT_LIFTED' for row in rows)
    unsat = sum(row['status'] == 'PROVED_UNSAT' for row in rows)
    decision = ('PANEL_PASS' if len(rows) == 32 and sat == 5 and unsat == 27
                else 'CENSORED_OR_FAILED')
    result = {'decision': decision, 'domain': frozen['domain'],
              'branch_order_sha256': sha(HERE / 'QUERY_MANIFEST.json'),
              'solver_binary_sha256': sha(solver),
              'checker_source_sha256': sha(source),
              'checker_binary_sha256': sha(archived_checker),
              'host_toolchain': host(), 'compile': compiled,
              'sat_lifted': sat, 'proved_unsat': unsat,
              'cases': rows}
    (out / 'result.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    progress(out, frozen, rows, compiled, status=decision)
    return result


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--internal-supervised', action='store_true')
    args = p.parse_args()
    frozen = json.loads((HERE / 'FROZEN.json').read_text())
    if not args.internal_supervised or frozen['release_main_head'] is None:
        raise SystemExit('HELD: merged-parent exact-head release and bounded supervisor required')
    from run import release_gate
    release_gate(frozen)
    lock_path = Path(frozen['panel_lock_path'])
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    # Hold this through fixed-path checker build/copy and all 32 labels.
    # A concurrent attempt fails before creating an outcome directory.
    with lock_path.open('a+b') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise SystemExit('HELD: concurrent n13 S5 panel owns fixed-path checker') from exc
        out = args.out.resolve()
        out.mkdir(parents=True, exist_ok=False)
        result = run(out, frozen)
    print(json.dumps({'decision': result['decision'],
                      'sat_lifted': result['sat_lifted'],
                      'proved_unsat': result['proved_unsat']}, sort_keys=True))
    return 0 if result['decision'] == 'PANEL_PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
