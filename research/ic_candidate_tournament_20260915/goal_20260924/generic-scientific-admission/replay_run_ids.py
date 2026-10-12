"""Re-export retained executions with unique numbers; never launch a worker."""
import argparse
import copy
import gzip
import json
from pathlib import Path
import sys

HARNESS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(HARNESS))
from execution_ids import allocation, audit_runs
from generic_admission import admit, admit_rho
from generic_build import digest, verify_binding
from generic_stages import verify_stages
from identity import canonical, write_immutable
from oracle import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    build_dir = bundle / 'build-final'
    build = json.loads((build_dir / 'build-record.json').read_text())
    source = json.loads((build_dir / 'source-manifest.json').read_text())
    rows, inputs = [], {}
    for panel in ('final', 'failed-la', 'sparse-core'):
        directory = bundle / f'controls-{panel}'
        resources = json.loads((directory / 'environment.json').read_text())['resources']
        raw = directory / 'worker-raw.jsonl'
        inputs[str(raw.relative_to(bundle))] = digest(raw)
        rows.extend((panel, item, resources) for item in
                    [json.loads(line) for line in raw.read_text().splitlines()])
    plan = allocation(0, [f'{panel}/{item["name"]}' for panel, item, _ in rows])
    exports, results, audit_ns = [], [], 0
    for (panel, item, resources), execution in zip(rows, plan['executions'], strict=True):
        job, report = item['job'], item['report']
        fixture = report['fixture']
        if job['mode'] == 'inventory':
            verify_binding(report, build, source, executable=build_dir / 'worker')
            verify_stages(report, fixture, job)
        elif job['mode'] == 'rho':
            admit_rho(report, fixture, job, build, source, executable=build_dir / 'worker',
                      process_wall_ns=item['process_wall_ns'])
        else:
            fresh = admit(report, fixture, job, build, source, executable=build_dir / 'worker',
                process_wall_ns=item['process_wall_ns'], resources=resources, number=execution['number'])
            original = item['receipt']['admission']
            require(fresh['candidate'] == original['candidate']
                    and fresh['workload'] == original['workload'], 'replay changed method/workload')
            audit_ns += fresh['run']['independent_audit_wall_ns']
            corrected = copy.deepcopy(original['run'])
            corrected['run_id'] = fresh['run']['run_id']
            fresh['run']['independent_audit_wall_ns'] = original['run']['independent_audit_wall_ns']
            require(fresh['run'] == corrected, 'run changed beyond identity/audit replay time')
            exports.append(dict(panel=panel, name=item['name'], original_run_id=original['run']['run_id'],
                candidate=original['candidate'], workload=original['workload'], run=corrected))
        results.append(dict(panel=panel, name=item['name'], status='PASS'))
    audit = audit_runs([item['run'] for item in exports])
    args.out.mkdir(parents=True, exist_ok=False)
    payload = canonical(dict(schema_version=2, scope='corrected exports; no new worker executions',
                             allocation=plan, records=exports)) + b'\n'
    target = args.out / 'run-records-v2.json.gz'
    target.write_bytes(gzip.compress(payload, mtime=0))
    receipt = dict(schema_version=1, controls=len(results), results=results, run_key_audit=audit,
        workers_rerun=False, measurements_changed=False, candidate_or_workload_changed=False,
        archive_unchanged=True, input_sha256=inputs, corrected_export_sha256=digest(target),
        correction_audit_wall_ns=audit_ns, original_audit_times_preserved=True,
        supersedes='R0 run keys in original bundle only; original bundle retained unchanged',
        checker_sources={name: digest(HARNESS / name) for name in ('execution_ids.py',
            'generic_admission.py', 'generic_build.py', 'generic_stages.py', 'generic_bases.py',
            'generic_phases.py', 'generic_queries.py', 'generic_query_law.py', 'identity.py',
            'measurement.py', 'oracle.py')})
    receipt['reexport_source_sha256'] = digest(Path(__file__))
    write_immutable(args.out / 'run-key-correction.json', receipt)
    print(json.dumps(dict(controls=len(results), **audit)))


if __name__ == '__main__':
    main()
