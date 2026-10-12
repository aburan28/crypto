#!/usr/bin/env python3
"""Run the frozen disclosed-point F4/F5/SAT recovery panel once.

This is a local scientific-admission pilot, never a competitive tournament.
Every process and audit failure is retained and no job is retried.
"""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

from generic_admission import admit
from generic_bases import verify_base
from generic_build import build, verify_build_record
from generic_solver_feasibility import assess, check_source_checkout
from oracle import require
from run_generic_backend_disclosed_pilot import digest, run_bounded_worker
from tournament import read, write

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / 'goal_20260924/generic-backend-recovery-pilot'
PANEL = REGISTRATION / 'panel.json'
# Filled from the committed panel bytes before any worker is executed.
PANEL_SHA256 = 'bf586fa75f9b5a00c6de95793c73502a2fb064898f1d4eb09c84685082b7b03d'
CELLS = ('n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0')
SOLVERS = ('f4', 'f5', 'sat_xor', 'sat_cnf')


def validate_panel(panel):
    require(digest(PANEL) == PANEL_SHA256, 'registered panel bytes changed')
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_MEASUREMENT'
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['inventory_sha256'] == 'be3053b4b637311e8255f294807fde412510fbda2464e678a6a16fdea397d403'
            and panel['cells'] == list(CELLS)
            and panel['solvers'] == list(SOLVERS)
            and panel['factor_base'] == {'kind': 'standard_subspace', 'dimension': 6}
            and panel['algorithm_seed'] == 2026092921
            and panel['config'] == dict(linear_algebra='dense', summands=3,
                                        groebner_degree=3, node_budget=4096,
                                        conflict_budget=100000)
            and panel['worker_threads'] == 1
            and panel['memory_bytes'] == 8*1024**3
            and panel['memory_policy'] == 'poll-child-rss-and-kill-above-limit'
            and panel['memory_poll_ms'] == 100
            and panel['per_cell_limits'] == {
                'n17a1': dict(max_trials=128, batch_trials=8, timeout_seconds=300),
                'n19a0': dict(max_trials=24, batch_trials=8, timeout_seconds=300),
                'n23a0': dict(max_trials=4, batch_trials=4, timeout_seconds=180),
                'n23a1': dict(max_trials=4, batch_trials=4, timeout_seconds=180),
                'n31a0': dict(max_trials=1, batch_trials=1, timeout_seconds=180),
            }, 'pilot schedule or limits changed')


def fixtures_from_inventory(panel):
    inventory_path = HERE / panel['inventory_file']
    require(digest(inventory_path) == panel['inventory_sha256'],
            'disclosed inventory control bytes changed')
    inventory = read(inventory_path)
    require(inventory['schema_version'] == 1
            and inventory['status'] == 'BASE_CONSTRUCTION_CONTROL_ONLY'
            and inventory['reviewed_generic_commit'] == panel['source_commit']
            and [row['cell'] for row in inventory['rows']] == list(CELLS),
            'inventory source or five disclosed cells changed')
    fixtures = {}
    for row in inventory['rows']:
        report, job = row['worker_report'], row['job']
        receipt = verify_base(report, report['fixture'], job)
        require(receipt == row['independent_receipt']
                and job['mode'] == 'inventory'
                and job['public_targets'] == report['fixture']['targets']
                and len(job['public_targets']) == 1
                and report['fixture']['target_scalar_constructed'] is False,
                'disclosed point or independently checked base changed')
        fixtures[row['cell']] = report['fixture']
    return fixtures


def static_preflight(panel):
    candidates = [dict(id=solver, config=dict(solver=solver, summands=3,
                   factor_base=panel['factor_base'])) for solver in SOLVERS]
    result = assess(dict(cells=list(CELLS), candidates=candidates))
    require(result['status'] == 'PASS_STATIC_LAYOUT_ONLY'
            and result['impossible_algebraic_cells'] == 0,
            'F4/F5 layout exceeds reviewed source cap')
    return result


def declared_job(panel, cell, solver, fixture):
    require(cell in CELLS and solver in SOLVERS, 'unknown registered job')
    limits = panel['per_cell_limits'][cell]
    return dict(mode='ic', degree=fixture['degree'], curve_a=fixture['curve_a'],
                public_targets=fixture['targets'], algorithm_seed=panel['algorithm_seed'],
                factor_base=panel['factor_base'],
                config=dict(panel['config'], solver=solver,
                            max_trials=limits['max_trials'],
                            batch_trials=limits['batch_trials']),
                exclusive_phases=True)


def resources(panel, cell):
    return dict(host_class='macos-arm64-local-stage-diagnostic',
                process_timeout_seconds=panel['per_cell_limits'][cell]['timeout_seconds'],
                sampled_child_rss_limit_bytes=panel['memory_bytes'],
                memory_poll_ms=panel['memory_poll_ms'], rayon_threads=1)


def host_manifest():
    cpu = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'],
                         text=True, capture_output=True, check=False)
    return dict(system=platform.system(), release=platform.release(),
                version=platform.version(), machine=platform.machine(),
                processor=platform.processor(), logical_cpus=os.cpu_count(),
                cpu_brand=cpu.stdout.strip() if cpu.returncode == 0 else None,
                cpu_query_exit=cpu.returncode,
                scope='local diagnostic host, not calibrated Linux instruction platform')


def audit_report(report, fixture, job, build_record, source, worker,
                 process_wall_ns, panel, cell, number):
    require(report['status'] in {'complete', 'incomplete'},
            'worker returned an error rather than an IC report')
    receipt = admit(report, fixture, job, build_record, source,
                    executable=worker, process_wall_ns=process_wall_ns,
                    resources=resources(panel, cell), number=number)
    attempts = [attempt for batch in report['collection_reports']
                for attempt in batch['attempts']]
    outcomes = Counter(attempt['pdp']['outcome'] for attempt in attempts)
    solver = job['config']['solver']
    family = 'sat' if solver.startswith('sat_') else 'groebner'
    nonidentity = [a for a in attempts if a['pdp']['outcome'] != 'identity']
    dispatched = any(a['pdp']['stats']['family'] == family
                     and a['pdp']['stats']['stats'].get('unsupported') is False
                     for a in nonidentity)
    matrix = receipt['stages']['matrix']
    run = receipt['run']
    return receipt, dict(dispatched=dispatched, worker_status=report['status'],
                         attempts=len(attempts), pdp_outcomes=dict(outcomes),
                         verified_relations=matrix['accepted_rows'],
                         final_rank=matrix['rank'],
                         actual_base_points=receipt['stages']['base']['inventory']['usable_point_count'],
                         folded_columns=receipt['stages']['base']['inventory']['effective_columns'],
                         candidate_id=run['candidate_id'], workload_id=run['workload_id'],
                         run_id=run['run_id'],
                         verified_complete=report['status'] == 'complete',
                         online_wall_ns=run['online_wall_ns'] if report['status'] == 'complete' else None,
                         cold_wall_ns=run['cold_wall_ns'] if report['status'] == 'complete' else None,
                         paired_rho_speedup=None, promotion_eligible=False,
                         audit_status='PASS')


def one_job(out, panel, cell, solver, fixture, worker, build_record, source, number):
    directory = out / 'jobs' / cell / solver
    directory.mkdir(parents=True, exist_ok=False)
    job = declared_job(panel, cell, solver, fixture)
    write(directory/'job.json', job, exclusive=True)
    limits = dict(panel, timeout_seconds=panel['per_cell_limits'][cell]['timeout_seconds'])
    stdout, stderr, exit_code, disposition, wall_ns, peak_rss = run_bounded_worker(
        worker, job, directory, limits)
    process = dict(cell=cell, solver=solver, exit_code=exit_code,
                   disposition=disposition, diagnostic_process_wall_ns=wall_ns,
                   sampled_peak_rss_bytes=peak_rss, memory_policy=panel['memory_policy'])
    write(directory/'process.json', process, exclusive=True)
    row = dict(process, dispatched=False, verified_complete=False,
               online_wall_ns=None, cold_wall_ns=None,
               paired_rho_speedup=None, promotion_eligible=False)
    if disposition == 'EXITED' or (disposition == 'PROCESS_FAILURE' and exit_code == 2):
        try:
            report = json.loads(stdout)
            require((exit_code == 0 and report['status'] == 'complete')
                    or (exit_code == 2 and report['status'] == 'incomplete'),
                    'worker exit code and report status disagree')
            receipt, diagnostics = audit_report(report, fixture, job, build_record,
                source, worker, wall_ns, panel, cell, number)
            write(directory/'admission.json', receipt, exclusive=True)
            row.update(diagnostics)
            row['disposition'] = ('VERIFIED_COMPLETE' if report['status'] == 'complete'
                                  else 'BOUNDED_INCOMPLETE_REPORT')
        except Exception as exc:
            row.update(audit_status='FAIL', reason=f'{type(exc).__name__}: {exc}')
    else:
        row.update(audit_status='NOT_AVAILABLE', reason='process did not exit with a report')
    write(directory/'result.json', row, exclusive=True)
    return row


def run(panel, source_root, out):
    validate_panel(panel)
    check_source_checkout(source_root)
    fixtures = fixtures_from_inventory(panel)
    preflight = static_preflight(panel)
    require(not out.exists(), 'output already exists; this registration cannot be retried')
    out.mkdir(parents=True)
    shutil.copy2(PANEL, out/'registered-panel.json')
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    write(out/'host.json', host_manifest(), exclusive=True)
    write(out/'static-preflight.json', preflight, exclusive=True)
    record = build(source_root, out/'build')
    source = read(out/'build/source-manifest.json')
    verify_build_record(record, source)
    worker = out/'build/worker'
    results = []
    for cell in CELLS:
        for solver in SOLVERS:
            row = one_job(out, panel, cell, solver, fixtures[cell], worker,
                          record, source, len(results) + 1)
            results.append(row)
            with (out/'progress.jsonl').open('a') as stream:
                stream.write(json.dumps(row, sort_keys=True) + '\n')
            print(json.dumps(dict(cell=cell, solver=solver,
                                  disposition=row['disposition'],
                                  audit_status=row['audit_status'],
                                  complete=row['verified_complete'],
                                  attempts=row.get('attempts'),
                                  rank=row.get('final_rank'))), flush=True)
    n17 = [row for row in results if row['cell'] == 'n17a1']
    algebra_complete = any(row['solver'] in {'f4', 'f5'} and row['verified_complete']
                           and row['audit_status'] == 'PASS' for row in n17)
    sat_complete = any(row['solver'] in {'sat_xor', 'sat_cnf'} and row['verified_complete']
                       and row['audit_status'] == 'PASS' for row in n17)
    summary = dict(schema_version=1, panel_sha256=PANEL_SHA256,
                   source_commit=panel['source_commit'],
                   source_manifest_sha256=record['source_manifest_sha256'],
                   worker_sha256=record['worker_sha256'], rows=results,
                   n17_f4_f5_complete=algebra_complete, n17_sat_complete=sat_complete,
                   status='DISCLOSED_RECOVERY_PASS' if algebra_complete and sat_complete
                          else 'DISCLOSED_RECOVERY_NOT_ESTABLISHED',
                   scope='disclosed-point scientific admission; no matched rho or promotion')
    write(out/'summary.json', summary, exclusive=True)
    print(json.dumps(dict(status=summary['status'], complete=sum(r['verified_complete'] for r in results),
                          jobs=len(results)), sort_keys=True), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run(read(PANEL), args.source_root.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
