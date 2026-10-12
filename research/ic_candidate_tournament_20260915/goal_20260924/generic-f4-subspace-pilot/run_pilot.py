#!/usr/bin/env python3
"""One disclosed-point F4/F5 dispatch pilot on the dimension-6 subspace.

This runner does not generate targets and does not call ``tournament.py
prepare``. It replays the five public points already stored in the
standard-subspace inventory control. ``promotion_eligible`` is always
false. Native seconds are never written into an instruction-count field.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
TOURNAMENT = HERE.parents[1]
REPO = TOURNAMENT.parents[1]
INVENTORY = (TOURNAMENT / 'goal_20260924/generic-backend-qualification-v2'
             / 'standard-subspace-d6-inventory-control.json')
INVENTORY_SHA256 = 'be3053b4b637311e8255f294807fde412510fbda2464e678a6a16fdea397d403'
SOURCE_COMMIT = '765c3c5f19032bd852163805f257c56babef2040'
PANEL = HERE / 'panel.json'
SOLVERS = ('f4', 'f5')
MAX_TRIALS = 1
BATCH_TRIALS = 1
GROEBNER_DEGREE = 3
NODE_BUDGET = 4096
CONFLICT_BUDGET = 100_000
WALL_SECONDS = 180
ADDRESS_SPACE_BYTES = 8 * 1024 ** 3
sys.path.insert(0, str(TOURNAMENT))

from generic_build import build, verify_binding  # noqa: E402
from oracle import require  # noqa: E402


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inventory():
    raw = INVENTORY.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == INVENTORY_SHA256,
            'dimension-6 inventory control bytes changed')
    document = json.loads(raw)
    require(document['reviewed_generic_commit'] == SOURCE_COMMIT
            and document['dimension'] == 6
            and document['status'] == 'BASE_CONSTRUCTION_CONTROL_ONLY',
            'inventory control is not the sealed dimension-6 base check')
    require([row['cell'] for row in document['rows']]
            == ['n17a1', 'n19a0', 'n23a0', 'n23a1', 'n31a0'],
            'inventory control cells changed')
    return document


def static_preflight(source_root, out):
    receipt = out / 'static-feasibility.json'
    process = subprocess.run(
        [sys.executable, str(TOURNAMENT / 'generic_solver_feasibility.py'),
         '--panel', str(PANEL), '--source-repo', str(source_root),
         '--out', str(receipt), '--require-pass'],
        cwd=REPO, text=True, capture_output=True, check=False)
    (out / 'static-feasibility.log').write_text(process.stdout + process.stderr)
    require(process.returncode == 0, 'encoder preflight failed; see static-feasibility.log')
    parsed = json.loads(receipt.read_text())
    require(parsed['status'] == 'PASS_STATIC_LAYOUT_ONLY', 'static layout did not pass')
    return parsed


def ensure_build(source_root, out):
    record_path = out / 'build' / 'build-record.json'
    if record_path.is_file():
        record = json.loads(record_path.read_text())
        source = json.loads((out / 'build' / 'source-manifest.json').read_text())
        verify_binding({'generic_build': record['identity']}, record, source,
                       executable=out / 'build' / 'worker')
        return record
    return build(source_root, out / 'build')


def job_for(row, solver):
    stored = row['job']
    require(stored['factor_base'] == {'dimension': 6, 'kind': 'standard_subspace'}
            and stored['public_targets'] and len(stored['public_targets']) == 1
            and stored['algorithm_seed'] == 43,
            f'inventory job for {row["cell"]} is not the sealed subspace point')
    return {
        'mode': 'ic',
        'degree': stored['degree'],
        'curve_a': stored['curve_a'],
        'algorithm_seed': stored['algorithm_seed'],
        'public_targets': stored['public_targets'],
        'factor_base': {'kind': 'standard_subspace', 'dimension': 6},
        'exclusive_phases': True,
        'config': {
            'solver': solver,
            'linear_algebra': 'dense',
            'batch_trials': BATCH_TRIALS,
            'max_trials': MAX_TRIALS,
            'summands': 3,
            'groebner_degree': GROEBNER_DEGREE,
            'node_budget': NODE_BUDGET,
            'conflict_budget': CONFLICT_BUDGET,
            'rho_parallel_walks': 32,
            'factor_base_cube_root': False,
        },
    }


def child_env():
    env = {key: os.environ[key] for key in ('PATH', 'HOME', 'TMPDIR') if key in os.environ}
    env.update(RAYON_NUM_THREADS='1', LC_ALL='C')
    return env


def limit_address_space():
    resource.setrlimit(resource.RLIMIT_AS, (ADDRESS_SPACE_BYTES, ADDRESS_SPACE_BYTES))


def run_one(worker, row, solver, out):
    cell = row['cell']
    slot = out / 'runs' / f'{cell}-{solver}'
    slot.mkdir(parents=True)
    job = job_for(row, solver)
    (slot / 'job.json').write_text(json.dumps(job, sort_keys=True) + '\n')
    started = time.perf_counter_ns()
    status = 'completed'
    try:
        process = subprocess.run(
            [str(worker)], input=json.dumps(job), text=True, capture_output=True,
            env=child_env(), timeout=WALL_SECONDS, check=False,
            preexec_fn=limit_address_space)
    except subprocess.TimeoutExpired as error:
        status = 'timeout'
        process = error
    wall_ns = time.perf_counter_ns() - started
    stdout = process.stdout or ''
    stderr = process.stderr or ''
    (slot / 'stdout.txt').write_text(stdout)
    (slot / 'stderr.txt').write_text(stderr)
    report = None
    if status == 'completed':
        try:
            report = json.loads(stdout)
        except json.JSONDecodeError:
            status = 'unparseable_stdout'
        else:
            (slot / 'report.json').write_text(json.dumps(report, sort_keys=True) + '\n')
    receipt = {
        'schema_version': 1,
        'cell': cell,
        'solver': solver,
        'promotion_eligible': False,
        'runner_status': status,
        'exit_code': getattr(process, 'returncode', None),
        'wall_ns': wall_ns,
        'instruction_count': None,
        'valgrind': None,
        'address_space_cap_bytes': ADDRESS_SPACE_BYTES,
        'rayon_threads': 1,
        'comparison_kind': 'factor-base-policy',
    }
    if report is not None:
        receipt.update(summarize_report(report))
    (slot / 'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    return receipt


def summarize_report(report):
    attempts = []
    for collection in report.get('collection_reports') or []:
        for attempt in collection.get('attempts') or []:
            pdp = attempt.get('pdp') or {}
            stats = pdp.get('stats') or {}
            nested = stats.get('stats') or {}
            attempts.append({
                'trial': attempt.get('trial'),
                'outcome': pdp.get('outcome'),
                'family': stats.get('family'),
                'engine': stats.get('engine'),
                'unsupported': nested.get('unsupported'),
                'exhausted': nested.get('exhausted'),
            })
    table = report.get('log_table_report') or {}
    return {
        'report_status': report.get('status'),
        'report_reason': report.get('reason'),
        'collector_dispatch': report.get('collector_dispatch'),
        'descent_dispatch': report.get('descent_dispatch'),
        'columns': report.get('columns'),
        'geometric_points': None if report.get('factor_base') is None else len(report['factor_base']),
        'accepted_relations': table.get('relations', report.get('accepted_relations')),
        'solve_attempts': report.get('solve_attempts', table.get('solve_attempts')),
        'solutions': report.get('solutions'),
        'attempts': attempts,
        'online_wall_ns': report.get('online_wall_ns'),
        'generic_phase_policy': report.get('generic_phase_policy'),
        'instruction_count': None,
    }


def host_record():
    cpu = 'unknown'
    for line in Path('/proc/cpuinfo').read_text().splitlines():
        if line.startswith('model name'):
            cpu = line.split(':', 1)[1].strip()
            break
    return {
        'system': platform.system(),
        'machine': platform.machine(),
        'cpu': cpu,
        'logical_cpus': os.cpu_count(),
        'rayon_threads': 1,
        'valgrind': shutil.which('valgrind'),
        'instruction_unit': None,
        'note': 'Valgrind is absent, so no instruction count is recorded. '
                'wall_seconds is practicality only.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    # This budget already has a canonical series (RESULT.md, PR #961) and a
    # retained replay (run-20260929/). A new directory is not a new registration.
    # The body below is the historical runner and is not reached.
    raise SystemExit(
        'dispatch closed: the max_trials=1 dimension-6 budget already ran. '
        'Canonical receipts are RESULT.md (PR #961). '
        'run-20260929/ is a same-budget replay and is not a new registration. '
        'A different budget needs a new protocol amendment and a new runner.'
    )
    source_root = args.source_root.resolve()
    out = args.out.resolve()
    require(not out.exists(), 'pilot output directory already exists; use a new directory')
    out.mkdir(parents=True)
    inventory = load_inventory()
    preflight = static_preflight(source_root, out)
    record = ensure_build(source_root, out)
    worker = out / 'build' / 'worker'
    receipts = [run_one(worker, row, solver, out)
                for row in inventory['rows'] for solver in SOLVERS]
    summary = {
        'schema_version': 1,
        'promotion_eligible': False,
        'comparison_kind': 'factor-base-policy',
        'source_commit': SOURCE_COMMIT,
        'inventory_sha256': INVENTORY_SHA256,
        'panel_sha256': digest(PANEL),
        'static_feasibility_status': preflight['status'],
        'worker_sha256': record['worker_sha256'],
        'build_sha256': record['build_sha256'],
        'rustc': record['build']['rustc'].splitlines()[0],
        'budget': {
            'solvers': list(SOLVERS),
            'max_trials': MAX_TRIALS,
            'batch_trials': BATCH_TRIALS,
            'groebner_degree': GROEBNER_DEGREE,
            'node_budget': NODE_BUDGET,
            'conflict_budget': CONFLICT_BUDGET,
            'linear_algebra': 'dense',
            'summands': 3,
            'factor_base': {'kind': 'standard_subspace', 'dimension': 6},
            'wall_cap_seconds': WALL_SECONDS,
            'address_space_cap_bytes': ADDRESS_SPACE_BYTES,
        },
        'host': host_record(),
        's_operations': None,
        'rho_online_ratio': None,
        'receipts': receipts,
    }
    text = json.dumps(summary, indent=2, sort_keys=True) + '\n'
    (out / 'summary.json').write_text(text)
    print(json.dumps({
        'promotion_eligible': False,
        'jobs': len(receipts),
        'report_status': [row.get('report_status') or row['runner_status'] for row in receipts],
        'worker_sha256': record['worker_sha256'],
        'summary_sha256': hashlib.sha256(text.encode()).hexdigest(),
    }, sort_keys=True))


if __name__ == '__main__':
    main()
