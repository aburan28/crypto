#!/usr/bin/env python3
"""Continue only the preregistered stage B after auditing frozen stage-A raw data.

The original stage-A summary is preserved. A valid incomplete worker report
exits 2; its raw content, rather than that exit code alone, determines whether
the original encoder-dispatch gate was met. No new inputs or limits are chosen.
"""
import argparse
import json
from pathlib import Path
import shutil

from generic_build import verify_binding, verify_build_record
from generic_solver_feasibility import check_source_checkout
from oracle import require
from run_generic_backend_disclosed_pilot import (
    CELLS, PANEL, PANEL_SHA256, REGISTRATION, SOLVERS, audit_report,
    declared_job, digest, fixture_map, one_job, static_preflight,
)
from tournament import read, write

ORIGINAL_SUMMARY_SHA256 = 'ff58ad86bebd13000347c5971f2c29da42374cd4ef1700cee4035be28e3e4a23'
BUILD_RECORD_SHA256 = 'fe6378068da024ffcfdd4d5d23d7d2c1221284ffd6fbe810ca28510177188a93'
STAGE_A_STDOUT_SHA256 = {
    'f4': '2f24cd7914634665c03dbf53102276bbcc4a6790c491c193381b3a129509a550',
    'f5': 'bde8532abe74a1c073600413349db16d4c318fdefed181a1d3bf5b65dd34dc46',
}
EVIDENCE = REGISTRATION / 'RESULT-v3-stage-a.json'


def validate_original(original, panel, fixtures):
    """Reject changed or missing original measurements before creating output."""
    require(original.is_dir(), 'missing original stage-A bundle')
    require(digest(original/'registered-panel.json') == PANEL_SHA256,
            'original registration differs from frozen panel')
    require(digest(original/'summary.json') == ORIGINAL_SUMMARY_SHA256,
            'original summary changed')
    require(digest(original/'build/build-record.json') == BUILD_RECORD_SHA256,
            'original build record changed')
    require(not any((original/'jobs'/cell).exists() for cell in CELLS[1:]),
            'stage B was already attempted in original bundle')
    summary = read(original/'summary.json')
    require(summary['stage_b_opened'] is False
            and len(summary['rows']) == len(CELLS)*len(SOLVERS)
            and all(row['disposition'] == 'NOT_RUN_BY_PREDECLARED_GATE'
                    for row in summary['rows'][2:]),
            'original summary does not show unrun stage B')
    record = read(original/'build/build-record.json')
    source = read(original/'build/source-manifest.json')
    verify_build_record(record, source)
    require(record['source_manifest_sha256'] == summary['source_manifest_sha256']
            and record['worker_sha256'] == summary['worker_sha256'],
            'original summary refers to another build')
    worker = original/'build/worker'
    require(digest(worker) == record['worker_sha256'], 'original executable changed')

    evidence = read(EVIDENCE)
    require(evidence['panel_sha256'] == PANEL_SHA256
            and evidence['original_summary_sha256'] == ORIGINAL_SUMMARY_SHA256
            and evidence['build_record_sha256'] == BUILD_RECORD_SHA256
            and evidence['source_manifest'] == source
            and evidence['build_record'] == record
            and evidence['original_summary'] == summary,
            'committed stage-A evidence differs from original')
    for solver in SOLVERS:
        directory = original/'jobs/n17a1'/solver
        frozen = evidence['raw_stage_a'][solver]
        require(digest(directory/'stdout.json') == STAGE_A_STDOUT_SHA256[solver]
                and frozen['stdout_sha256'] == STAGE_A_STDOUT_SHA256[solver],
                'stage-A stdout changed')
        job = read(directory/'job.json')
        process = read(directory/'process.json')
        result = read(directory/'result.json')
        stdout = read(directory/'stdout.json')
        stderr = (directory/'stderr.txt').read_text()
        require(job == declared_job(solver, fixtures['n17a1'], panel)
                and frozen['job'] == job and frozen['process'] == process
                and frozen['result'] == result and frozen['stdout'] == stdout
                and frozen['stderr'] == stderr,
                'stage-A raw record differs from registered evidence')
        require(process['exit_code'] == 2 and process['disposition'] == 'PROCESS_FAILURE'
                and 0 < process['diagnostic_process_wall_ns'] <= panel['timeout_seconds']*10**9
                and process['sampled_peak_rss_bytes'] <= panel['memory_bytes'],
                'stage-A worker did not exit inside resource cap')
        require(stderr == '' and stdout['status'] == 'incomplete',
                'stage-A report is not a clean bounded incomplete result')
        verify_binding(stdout, record, source, executable=worker)
    return record, source


def continue_pilot(original, out, source_root):
    require(digest(PANEL) == PANEL_SHA256, 'registered panel changed')
    panel = read(PANEL)
    require(panel['schema_version'] == 3
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['stage_a_cells'] == ['n17a1']
            and panel['stage_b_cells'] == list(CELLS[1:])
            and panel['solvers'] == list(SOLVERS)
            and panel['algorithm_seed'] == 2026092918
            and panel['timeout_seconds'] == 60
            and panel['memory_bytes'] == 8*1024**3
            and panel['memory_poll_ms'] == 100,
            'registered schedule or resource envelope changed')
    check_source_checkout(source_root)
    fixtures = fixture_map(panel)
    static_preflight(panel)
    record, source = validate_original(original, panel, fixtures)
    require(not out.exists(), 'continuation output already exists; do not retry')
    out.mkdir(parents=True)
    shutil.copytree(original/'build', out/'build')
    shutil.copytree(original/'jobs/n17a1', out/'jobs/n17a1')
    for name in ('registered-panel.json', 'registered-runner.py',
                 'static-preflight.json', 'PROTOCOL.md'):
        shutil.copy2(original/name, out/name)
    shutil.copy2(original/'summary.json', out/'summary-original.json')
    shutil.copy2(Path(__file__), out/'continuation-runner.py')
    worker = out/'build/worker'
    rows = []
    for solver in SOLVERS:
        directory = out/'jobs/n17a1'/solver
        job = read(directory/'job.json')
        report = read(directory/'stdout.json')
        row = dict(read(directory/'process.json'))
        row.update(audit_report(report, job, worker, record, source, directory))
        require(row['audit_status'] == 'PASS' and row['dispatched'],
                'original stage-A dispatch gate was not met')
        row['disposition'] = 'BOUNDED_INCOMPLETE_REPORT'
        write(directory/'result-repaired.json', row, exclusive=True)
        rows.append(row)
    write(out/'stage-a-repair.json', dict(
        original_summary_sha256=ORIGINAL_SUMMARY_SHA256,
        stage_a_stdout_sha256=STAGE_A_STDOUT_SHA256,
        reason='Original runner discarded valid incomplete reports at exit code 2.',
        gate='PASS'), exclusive=True)
    print(json.dumps(dict(stage='A_REAUDITED', dispatched=2)), flush=True)
    for cell in panel['stage_b_cells']:
        for solver in SOLVERS:
            row = one_job(out, cell, solver, fixtures[cell], panel,
                          worker, record, source)
            rows.append(row)
            with (out/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps(row, sort_keys=True) + '\n')
            print(json.dumps(dict(stage='B', cell=cell, solver=solver,
                                  disposition=row['disposition'],
                                  dispatched=row['dispatched'],
                                  outcome=row.get('pdp_outcome'))), flush=True)
    summary = dict(schema_version=1, panel_sha256=PANEL_SHA256,
                   source_commit=panel['source_commit'],
                   source_manifest_sha256=record['source_manifest_sha256'],
                   worker_sha256=record['worker_sha256'],
                   original_summary_sha256=ORIGINAL_SUMMARY_SHA256,
                   worker_environment_policy='PATH/HOME/TMPDIR/LANG/LC_ALL plus RAYON_NUM_THREADS=1',
                   stage_b_opened=True, rows=rows,
                   status=('ALL_DISPATCHED_AUDITED' if all(row['dispatched'] for row in rows)
                           else 'FEASIBILITY_NOT_ESTABLISHED'),
                   scope='disclosed-target one-query stage control; no competitive result')
    write(out/'summary.json', summary, exclusive=True)
    print(json.dumps(dict(status=summary['status'], stage_b_opened=True,
                          dispatched=sum(row['dispatched'] for row in rows),
                          scheduled=len(rows)), sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--original', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    continue_pilot(args.original.resolve(), args.out.resolve(), args.source_root.resolve())


if __name__ == '__main__':
    main()
