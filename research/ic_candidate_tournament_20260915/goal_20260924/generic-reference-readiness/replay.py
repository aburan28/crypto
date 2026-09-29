"""Independently replay the registered readiness evidence without running workers."""
import argparse
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys

HARNESS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HARNESS))
from execution_ids import audit_runs
from generic_admission import admit, admit_rho
from generic_build import digest
from identity import sha256, write_immutable
from oracle import require


def replay(evidence):
    controls, build_dir = evidence / 'controls', evidence / 'build'
    build = json.loads((build_dir / 'build-record.json').read_text())
    source = json.loads((build_dir / 'source-manifest.json').read_text())
    environment = json.loads((controls / 'environment.json').read_text())
    resources = environment['resources']
    for name, expected_hash in environment['sources'].items():
        require(digest(HARNESS / name) == expected_hash, 'recorded checker source mismatch: ' + name)
    spec = importlib.util.spec_from_file_location('registered_controls',
        HARNESS / 'goal_20260924/generic-scientific-admission/run_controls.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    expected = module.reference_readiness_panel()
    rows = [json.loads(line) for line in (controls / 'worker-raw.jsonl').read_text().splitlines()]
    require(len(rows) == len(expected) == 15, 'readiness panel incomplete')
    records, runs = [], []
    for ordinal, (item, vector) in enumerate(zip(rows, expected, strict=True)):
        job, report = item['job'], item['report']
        require(item['name'] == vector['name'] and sha256({k: v for k, v in job.items()
            if k != 'public_targets'}) == sha256(vector['job']), 'registered job mismatch')
        require(report['status'] == 'complete', 'readiness solve incomplete')
        kwargs = dict(executable=build_dir / 'worker', process_wall_ns=item['process_wall_ns'])
        if job['mode'] == 'ic':
            fresh = admit(report, report['fixture'], job, build, source,
                          resources=resources, number=2000+ordinal, **kwargs)
            fresh['run']['independent_audit_wall_ns'] = item['receipt']['admission']['run']['independent_audit_wall_ns']
            runs.append(fresh['run'])
        else:
            fresh = admit_rho(report, report['fixture'], job, build, source, **kwargs)
        require(sha256(fresh) == sha256(item['receipt']['admission']), 'admission replay differs')
        record = dict(name=item['name'], cell=f'n{job["degree"]}a{job["curve_a"]}',
            mode=job['mode'], status='PASS', worker_status=report['status'],
            performance_qualified=False, online_speedup=None, normalized_S=None)
        if job['mode'] == 'ic':
            counts = Counter(a['pdp']['outcome'] for batch in report['collection_reports']
                             for a in batch['attempts'])
            target_counts = Counter(a['pdp']['outcome'] for s in report['solutions'] for a in s['attempts'])
            inv = fresh['stages']['base']['inventory']
            record.update(candidate_id=fresh['run']['candidate_id'], workload_id=fresh['run']['workload_id'],
                run_id=fresh['run']['run_id'], actual_usable_points=inv['usable_point_count'],
                effective_columns=inv['effective_columns'], final_rank=report['columns'],
                collection_outcomes=dict(counts), target_outcomes=dict(target_counts),
                relation_la=job['config']['linear_algebra'])
        records.append(record)
    return dict(schema_version=1, status='PASS', controls=15, complete_ic=10, complete_rho=5,
        run_key_audit=audit_runs(runs), workers_rerun=False, performance_qualified=False,
        promotion_eligible=False, records=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    receipt = replay(args.evidence.resolve())
    write_immutable(args.out, receipt)
    print(json.dumps({k: v for k, v in receipt.items() if k != 'records'}))
