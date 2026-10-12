#!/usr/bin/env python3
"""Classify every retained ordinary PDP query with an independent exact group oracle.

The input is the immutable disclosed recovery evidence; this does not run or
restart any measured worker and produces no competitive performance claim.
"""
from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile
import tempfile

from generic_build import verify_build_record
from oracle import Curve, require
from run_generic_backend_recovery_pilot import (
    CELLS, PANEL as PARENT_PANEL, SOLVERS, audit_report, declared_job,
    fixtures_from_inventory, validate_panel,
)
from tournament import read, write

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE / 'goal_20260924/generic-exact-yield-audit'
PANEL = REGISTRATION / 'panel.json'
PANEL_SHA256 = 'f795b7c896aca031eb37e630c14767b36e7f3922cef7f20fe9438748069cad67'


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def validate_registration(panel):
    require(digest_bytes(PANEL.read_bytes()) == PANEL_SHA256, 'registered panel changed')
    require(panel['schema_version'] == 1
            and panel['status'] == 'REGISTERED_BEFORE_ANALYSIS'
            and panel['parent_evidence_sha256']
                == 'e0ce19fc28c58e2dc1cae9649a16af74099016fff8183e5ce69f65b39c804f02'
            and panel['parent_panel_sha256']
                == 'bf586fa75f9b5a00c6de95793c73502a2fb064898f1d4eb09c84685082b7b03d'
            and panel['source_commit'] == '765c3c5f19032bd852163805f257c56babef2040'
            and panel['cells'] == list(CELLS)
            and panel['solvers'] == list(SOLVERS)
            and panel['summands'] == 3,
            'registered source, schedule or mathematical question changed')


def load_evidence(panel):
    archive_path = HERE / panel['parent_evidence_file']
    data = archive_path.read_bytes()
    require(digest_bytes(data) == panel['parent_evidence_sha256'],
            'source evidence archive changed')
    with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
        members = archive.getmembers()
        names = [m.name for m in members]
        require(len(names) == len(set(names))
                and all(m.isfile() and not m.name.startswith('/')
                        and '..' not in Path(m.name).parts for m in members),
                'unsafe or duplicate evidence archive entry')
        return {m.name: archive.extractfile(m).read() for m in members}


def record(files, name):
    require(name in files, f'missing retained evidence: {name}')
    return json.loads(files[name])


def pair_index(curve, base):
    """Every unordered pair is enough because the group law is commutative."""
    pairs = {}
    for i, first in enumerate(base):
        for j in range(i, len(base)):
            pairs.setdefault(curve.add(first, base[j]), (i, j))
    return pairs


def exact_three_sum(curve, base, pairs, target):
    """Return an independently checked triple, or None after exhaustive search."""
    for k, point in enumerate(base):
        pair = pairs.get(curve.add(target, curve.neg(point)))
        if pair is not None:
            indices = (*pair, k)
            reconstructed = None
            for index in indices:
                reconstructed = curve.add(reconstructed, base[index])
            require(reconstructed == target, 'exact-oracle witness does not re-add')
            return list(indices)
    return None


def audit(panel):
    validate_registration(panel)
    files = load_evidence(panel)
    parent = read(PARENT_PANEL)
    validate_panel(parent)
    require(files['registered-panel.json'] == PARENT_PANEL.read_bytes()
            and digest_bytes(files['registered-panel.json']) == panel['parent_panel_sha256'],
            'retained original panel changed')
    fixtures = fixtures_from_inventory(parent)
    summary = record(files, 'summary.json')
    require([(r['cell'], r['solver']) for r in summary['rows']]
            == [(c, s) for c in CELLS for s in SOLVERS],
            'missing or reordered original job')
    require(summary['source_commit'] == panel['source_commit'],
            'original worker source commit changed')
    source = record(files, 'build/source-manifest.json')
    build = record(files, 'build/build-record.json')
    verify_build_record(build, source)
    with tarfile.open(fileobj=io.BytesIO(files['build/root-source.tar.gz']),
                      mode='r:gz') as source_archive:
        require(set(source_archive.getnames()) == set(source['root_files']),
                'retained root source member list changed')
        for name, expected in source['root_files'].items():
            require(digest_bytes(source_archive.extractfile(name).read()) == expected,
                    f'retained root source changed: {name}')
    require(summary['source_manifest_sha256'] == build['source_manifest_sha256']
            and summary['worker_sha256'] == build['worker_sha256'],
            'original summary and build disagree')
    rows, per_query = [], {}
    geometries, query_prefixes = {}, {}
    with tempfile.TemporaryDirectory() as temp:
        worker = Path(temp) / 'worker'
        worker.write_bytes(files['build/worker'])
        require(digest_bytes(worker.read_bytes()) == build['worker_sha256'],
                'retained worker executable changed')
        for number, original in enumerate(summary['rows'], 1):
            cell, solver = original['cell'], original['solver']
            prefix = f'jobs/{cell}/{solver}/'
            require(record(files, prefix+'result.json') == original,
                    'original raw job and summary disagree')
            process = record(files, prefix+'process.json')
            require(process['diagnostic_process_wall_ns']
                    == original['diagnostic_process_wall_ns'],
                    'original process interval changed')
            job = record(files, prefix+'job.json')
            require(job == declared_job(parent, cell, solver, fixtures[cell])
                    and job['config']['summands'] == panel['summands'],
                    'original exact job changed')
            row = dict(cell=cell, solver=solver, original_run_id=original.get('run_id'),
                       original_status=original['disposition'],
                       attempts=None, observed_witnesses=None,
                       bounded_incomplete=None, exact_feasible=None,
                       feasible_but_incomplete=None, exact_infeasible=None,
                       final_rank=original.get('final_rank'),
                       analysis_status='UNKNOWN_NO_REPORT')
            if original['disposition'] == 'TIMEOUT':
                require(files[prefix+'stdout.json'] == b''
                        and prefix+'admission.json' not in files
                        and original['audit_status'] == 'NOT_AVAILABLE',
                        'timeout contains an unexamined report')
                rows.append(row)
                continue
            report = record(files, prefix+'stdout.json')
            receipt, diagnostic = audit_report(
                report, fixtures[cell], job, build, source, worker,
                process['diagnostic_process_wall_ns'], parent, cell, number)
            retained = record(files, prefix+'admission.json')
            receipt['run'].pop('independent_audit_wall_ns')
            retained['run'].pop('independent_audit_wall_ns')
            require(receipt == retained
                    and all(original[k] == v for k, v in diagnostic.items()),
                    'original independent scientific admission failed on replay')
            curve, base, pairs = geometries.get(cell, (None, None, None))
            if curve is None:
                curve = Curve(fixtures[cell])
                base = tuple(curve.decode(p) for p in report['factor_base'])
                require(base and all(point is not None for point in base),
                        'invalid geometric base')
                pairs = pair_index(curve, base)
                geometries[cell] = curve, base, pairs
            else:
                require(base == tuple(curve.decode(p) for p in report['factor_base']),
                        'different geometric base across solver arms')
            attempts = [a for batch in report['collection_reports']
                        for a in batch['attempts']]
            sequence = [(a['a'], a['b']) for a in attempts]
            previous = query_prefixes.get(cell, sequence)
            require(sequence[:min(len(sequence), len(previous))]
                    == previous[:min(len(sequence), len(previous))],
                    'solver arms used different ordinary query prefixes')
            query_prefixes[cell] = sequence if len(sequence) > len(previous) else previous
            labels = []
            for attempt in attempts:
                require(attempt['b'] == 0, 'target-dependent query in ordinary collection')
                scalar = attempt['a']
                key = cell, scalar
                if key not in per_query:
                    target = curve.mul(curve.g, scalar)
                    require(target is not None, 'identity ordinary query')
                    per_query[key] = exact_three_sum(curve, base, pairs, target)
                triple = per_query[key]
                outcome = attempt['pdp']['outcome']
                require(not (outcome == 'witness' and triple is None)
                        and not (outcome == 'proved_unsat' and triple is not None),
                        'measured PDP verdict contradicts exact group oracle')
                labels.append(dict(trial=attempt['trial'], a=scalar, b=0,
                                   measured_outcome=outcome,
                                   exact_feasible=triple is not None,
                                   exact_witness_indices=triple))
            outcomes = Counter(label['measured_outcome'] for label in labels)
            feasible = sum(label['exact_feasible'] for label in labels)
            missed = sum(label['exact_feasible'] and label['measured_outcome'] == 'incomplete'
                         for label in labels)
            row.update(attempts=len(labels), observed_witnesses=outcomes['witness'],
                       bounded_incomplete=outcomes['incomplete'], exact_feasible=feasible,
                       feasible_but_incomplete=missed,
                       exact_infeasible=len(labels)-feasible,
                       geometric_base_points=len(base),
                       usable_base_points=original['actual_base_points'],
                       folded_columns=original['folded_columns'],
                       analysis_status='PASS_EXACT_GROUP_ORACLE',
                       queries=labels)
            rows.append(row)
    require(len(rows) == 20 and sum(r['analysis_status'] == 'PASS_EXACT_GROUP_ORACLE'
                                   for r in rows) == 16,
            'original 16-report/four-timeout schedule changed')
    return dict(schema_version=1, status='EXACT_DISCLOSED_QUERY_AUDIT',
                scope='retained ordinary queries only; no new worker run or speedup',
                panel_sha256=PANEL_SHA256,
                parent_evidence_sha256=panel['parent_evidence_sha256'],
                source_manifest_sha256=build['source_manifest_sha256'],
                worker_sha256=build['worker_sha256'],
                rows=rows, exact_distinct_queries=len(per_query),
                timeout_rows=4, audited_report_rows=16)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    require(not out.exists(), 'audit output already exists; do not overwrite evidence')
    panel = read(PANEL)
    result = audit(panel)
    out.mkdir(parents=True)
    shutil.copy2(PANEL, out/'registered-panel.json')
    shutil.copy2(REGISTRATION/'PROTOCOL.md', out/'PROTOCOL.md')
    shutil.copy2(Path(__file__), out/'registered-runner.py')
    result['runner_sha256'] = digest_bytes((out/'registered-runner.py').read_bytes())
    write(out/'RESULT.json', result, exclusive=True)
    print(json.dumps(dict(status=result['status'], audited=result['audited_report_rows'],
                          timeouts=result['timeout_rows'],
                          distinct_queries=result['exact_distinct_queries'])), flush=True)


if __name__ == '__main__':
    main()
