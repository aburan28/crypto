#!/usr/bin/env python3
"""One-shot, disclosed-target F4/F5 encoder-dispatch pilot (stage costs only)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from generic_bases import verify_base
from generic_build import build, verify_binding, verify_build_record
from generic_query_law import verify_query_law
from generic_solver_feasibility import assess, check_source_checkout
from oracle import require
from tournament import read, write

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / 'goal_20260924/generic-backend-disclosed-pilot'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = '48db71af3c5f297b834f2476056201caec8378f96f6da6992d1aa8ecbe523685'
CELLS = ('n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0')
SOLVERS = ('f4', 'f5')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fixture_map(panel):
    path = HERE / panel['exposed_fixture_file']
    require(digest(path) == panel['exposed_fixture_sha256'],
            'disclosed fixture corpus changed')
    exposures = read(path)
    cases = [row for row in exposures['attempts']
             if row['stage'] == panel['fixture_stage'] and row['accepted']]
    require(len(cases) == 5, 'expected five disclosed A/A points')
    result = {f"n{case['fixture']['degree']}a{case['fixture']['curve_a']}": case['fixture']
              for case in cases}
    require(set(result) == set(CELLS), 'changed disclosed curve cells')
    return result


def static_preflight(panel):
    candidates = [dict(id=solver, config=dict(solver=solver,
                   summands=panel['config']['summands'], factor_base=panel['factor_base']))
                  for solver in SOLVERS]
    result = assess(dict(cells=list(CELLS), candidates=candidates))
    require(result['status'] == 'PASS_STATIC_LAYOUT_ONLY',
            'dimension-six pilot exceeds algebraic encoder cap')
    return result


def run_bounded_worker(worker, job, directory, panel):
    """Preserve raw streams and apply a diagnostic macOS-compatible RSS watch."""
    stdout_path, stderr_path = directory/'stdout.json', directory/'stderr.txt'
    # Exclusive-phase workers reject every undeclared KIC_ / F4_ / SOLVER_
    # override. Do not inherit algorithm controls from the interactive shell.
    environment = {key: os.environ[key] for key in
                   ('PATH', 'HOME', 'TMPDIR', 'LANG', 'LC_ALL') if key in os.environ}
    environment['RAYON_NUM_THREADS'] = '1'
    ps = shutil.which('ps')
    require(ps is not None, 'RSS monitor requires ps')
    started = time.monotonic_ns()
    disposition, peak_rss = None, 0
    with stdout_path.open('x') as stdout, stderr_path.open('x') as stderr:
        process = subprocess.Popen([str(worker)], stdin=subprocess.PIPE,
                                   stdout=stdout, stderr=stderr, text=True,
                                   env=environment)
        try:
            process.stdin.write(json.dumps(job))
        except BrokenPipeError:
            pass
        finally:
            process.stdin.close()
        while process.poll() is None:
            if (time.monotonic_ns() - started) >= panel['timeout_seconds'] * 10**9:
                disposition = 'TIMEOUT'
                process.kill()
                break
            sampled = subprocess.run([ps, '-o', 'rss=', '-p', str(process.pid)],
                                     text=True, capture_output=True, check=False)
            if sampled.returncode == 0 and sampled.stdout.strip():
                rss = int(sampled.stdout.strip()) * 1024
                peak_rss = max(peak_rss, rss)
                if rss > panel['memory_bytes']:
                    disposition = 'RSS_LIMIT'
                    process.kill()
                    break
            time.sleep(panel['memory_poll_ms'] / 1000)
        exit_code = process.wait()
    if disposition is None:
        disposition = 'EXITED' if exit_code == 0 else 'PROCESS_FAILURE'
    return (stdout_path.read_text(), stderr_path.read_text(), exit_code,
            disposition, time.monotonic_ns() - started, peak_rss)


def one_job(out, cell, solver, fixture, panel, worker, record, source):
    directory = out / 'jobs' / cell / solver
    directory.mkdir(parents=True, exist_ok=False)
    job = dict(mode='ic', degree=fixture['degree'], curve_a=fixture['curve_a'],
               public_targets=fixture['targets'], algorithm_seed=panel['algorithm_seed'],
               factor_base=panel['factor_base'],
               config=dict(panel['config'], solver=solver), exclusive_phases=True)
    write(directory/'job.json', job, exclusive=True)
    stdout, stderr, exit_code, disposition, wall_ns, peak_rss = run_bounded_worker(
        worker, job, directory, panel)
    row = dict(cell=cell, solver=solver, disposition=disposition, exit_code=exit_code,
               diagnostic_process_wall_ns=wall_ns, sampled_peak_rss_bytes=peak_rss,
               memory_policy=panel['memory_policy'], dispatched=False,
               competitive_total=None, competitive_speedup=None)
    write(directory/'process.json', row, exclusive=True)
    if disposition == 'EXITED':
        try:
            report = json.loads(stdout)
            build_receipt = verify_binding(report, record, source, executable=worker)
            base_receipt = verify_base(report, report['fixture'], job)
            query_receipt = verify_query_law(report, report['fixture'], job)
            write(directory/'build-audit.json', build_receipt, exclusive=True)
            write(directory/'base-audit.json', base_receipt, exclusive=True)
            write(directory/'query-audit.json', query_receipt, exclusive=True)
            attempts = [attempt for batch in report['collection_reports']
                        for attempt in batch['attempts']]
            require(len(attempts) == 1 and query_receipt['collection_queries'] == 1,
                    'pilot must contain exactly one natural ordinary query')
            pdp = attempts[0]['pdp']
            stats = pdp['stats']
            engine = 'MatrixF4' if solver == 'f4' else 'MatrixF5'
            require(stats['family'] == 'groebner'
                    and stats['engine'] == {engine: {'max_degree': 3}},
                    'wrong algebraic dispatch engine')
            native = stats['stats']
            row.update(pdp_outcome=pdp['outcome'], pdp_stats=native,
                       usable_base_points=base_receipt['inventory']['usable_point_count'],
                       folded_columns=base_receipt['inventory']['effective_columns'],
                       dispatched=(native['unsupported'] is False
                                   and pdp['outcome'] in ('witness', 'proved_unsat', 'incomplete')),
                       audit_status='PASS')
        except (ValueError, KeyError, TypeError) as exc:
            row.update(audit_status='FAIL', reason=str(exc))
    else:
        row.update(audit_status='NOT_AVAILABLE', reason='process did not exit with a report')
    write(directory/'result.json', row, exclusive=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    require(digest(PANEL) == PANEL_SHA256, 'pilot panel bytes changed after registration')
    panel = read(PANEL)
    require(panel['schema_version'] == 3
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['fixture_stage'] == 'smoke'
            and panel['algorithm_seed'] == 2026092918
            and panel['stage_a_cells'] == ['n17a1']
            and panel['stage_b_cells'] == list(CELLS[1:])
            and panel['solvers'] == list(SOLVERS)
            and panel['memory_policy'] == 'poll-child-rss-and-kill-above-limit'
            and panel['memory_poll_ms'] == 100, 'changed pilot schedule')
    check_source_checkout(args.source_root)
    fixtures = fixture_map(panel)
    static = static_preflight(panel)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copy2(PANEL, out/'registered-panel.json')
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    write(out/'static-preflight.json', static, exclusive=True)
    record = build(args.source_root, out/'build')
    source = read(out/'build/source-manifest.json')
    verify_build_record(record, source)
    worker = out/'build/worker'
    results = []
    for cell in panel['stage_a_cells']:
        for solver in SOLVERS:
            results.append(one_job(out, cell, solver, fixtures[cell], panel,
                                   worker, record, source))
    stage_b = all(row['dispatched'] and row['audit_status'] == 'PASS' for row in results)
    for cell in panel['stage_b_cells']:
        for solver in SOLVERS:
            if stage_b:
                results.append(one_job(out, cell, solver, fixtures[cell], panel,
                                       worker, record, source))
            else:
                results.append(dict(cell=cell, solver=solver,
                                    disposition='NOT_RUN_BY_PREDECLARED_GATE',
                                    dispatched=False, audit_status='NOT_RUN',
                                    competitive_total=None, competitive_speedup=None))
    summary = dict(schema_version=1, panel_sha256=PANEL_SHA256,
                   source_commit=panel['source_commit'],
                   source_manifest_sha256=record['source_manifest_sha256'],
                   worker_sha256=record['worker_sha256'],
                   worker_environment_policy='PATH/HOME/TMPDIR/LANG/LC_ALL plus RAYON_NUM_THREADS=1',
                   stage_b_opened=stage_b, rows=results,
                   status='ALL_DISPATCHED_AUDITED' if all(r['dispatched'] for r in results)
                          else 'FEASIBILITY_NOT_ESTABLISHED',
                   scope='disclosed-target one-query stage control; no competitive result')
    write(out/'summary.json', summary, exclusive=True)
    print(json.dumps(dict(status=summary['status'], stage_b_opened=stage_b,
                          dispatched=sum(r['dispatched'] for r in results),
                          scheduled=len(results)), sort_keys=True))


if __name__ == '__main__':
    main()
