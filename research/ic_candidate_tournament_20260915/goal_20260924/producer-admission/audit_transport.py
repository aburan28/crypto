"""Reconstruct every successful producer record from transported raw evidence."""
import argparse
import json
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument('artifact', type=Path)
args = parser.parse_args()
root = args.artifact.resolve()
prepared = root/'ic-producer'
evidence = root/'ic-producer-evidence'
sys.path.insert(0, str(evidence/'evaluator'))
from identity import candidate_manifest, sha256, workload_manifest
from measurement import measured_run, report_sha256
from oracle import require, verify
from producer.evidence import audit_stages, check_build_identity, method_record, scientific_ledger
from tournament import digest, parse_profiles, read

p = read(evidence/'protocol.json')
s = read(evidence/'summary.json')
preparation = read(prepared/'preparation.json')
manifest = read(prepared/'source-manifest.json')
source = sha256(manifest)
require(p['source_manifest_sha256'] == source == preparation['source_manifest_sha256'], 'source identity mismatch')
require(p['preparation'] == preparation, 'preparation mismatch')
require(digest(evidence/'worker') == p['worker_sha256'], 'binary changed in transport')
require(digest(evidence/'evaluator/producer/PROTOCOL.md') == p['protocol_sha256'], 'protocol changed')
for name, expected in p['evaluator_sha256'].items():
    require(digest(evidence/'evaluator'/name) == expected, 'evaluator changed: '+name)
require({str(f.relative_to(prepared/'source')) for f in (prepared/'source').rglob('*') if f.is_file()} == set(manifest), 'source inventory changed')
for name, expected in manifest.items():
    require(digest(prepared/'source'/name) == expected, 'source bytes changed: '+name)
reference = preparation['reference']
cells = [(13,0),(17,1),(23,1)] if reference == 'both' else [(13,0),(23,1),(37,0),(43,1),(61,1)]
require([(j['degree'],j['curve_a']) for j in p['jobs']] == cells, 'changed cell panel')
require(p['repetitions'] == 3 and p['resources']['worker_threads'] == 1, 'changed repetitions/resources')
require(s['verified'] == s['scheduled'] == len(cells)*3 and not s['promotion_eligible'], 'incomplete or promotional summary')
rows = []
for job in p['jobs']:
    cell = f"n{job['degree']}a{job['curve_a']}"
    directory = evidence/cell
    admission = read(directory/'admission/stdout.json')
    check_build_identity(admission,source)
    fixture = admission['fixture']
    method = method_record(job,fixture,manifest,source,reference,dict(p['build'],field_kernel=admission['field_kernel']))
    require(read(directory/'method.json') == method, 'method changed')
    candidate = candidate_manifest(fixture,admission,method)
    require(read(directory/'candidate.json') == candidate, 'candidate changed')
    workload = workload_manifest(fixture,input_law='ic-workflow-public-target-v1; hash-to-curve-cofactor-projection',algorithm_seed=job['algorithm_seed'],resource_envelope=p['resources'])
    require(read(directory/'workload.json') == workload, 'workload changed')
    for repetition in range(3):
        trial = directory/f'rep-{repetition}'
        require(read(trial/'job.json') == job, 'job changed')
        receipt = read(trial/'receipt.json')
        require(receipt['status'] == 'VERIFIED' and not receipt['promotion_eligible'], 'failed/promotional receipt')
        require({str(f.relative_to(trial)) for f in trial.rglob('*') if f.is_file()} == set(receipt['artifacts'])|{'receipt.json'}, 'trial inventory changed')
        for name,expected in receipt['artifacts'].items():
            require(digest(trial/name) == expected, 'trial bytes changed: '+name)
        native = read(trial/'native/process.json')
        profiled = read(trial/'profile/process.json')
        for process in (native,profiled):
            require(process['exit_code'] == 0 and process['process_status'] == 'EXITED','failed child')
            require(process['cpu'] == p['resources']['cpu'] and process['memory_cap_bytes'] == p['resources']['memory_bytes'], 'changed resource envelope')
        nr, pr = read(trial/'native/stdout.json'), read(trial/'profile/stdout.json')
        for report in (nr,pr):
            check_build_identity(report,source,admission['field_kernel'])
        require(verify(nr,fixture,summands=3) == verify(pr,fixture,summands=3), 'native/profile proofs differ')
        audit = audit_stages(pr,fixture,job['algorithm_seed'])
        require(audit_stages(nr,fixture,job['algorithm_seed']) == audit, 'native/profile diagnostics differ')
        ledger = scientific_ledger(pr,parse_profiles(trial/'profile',phase_schema=2))
        run = measured_run(candidate=candidate,workload=workload,report=pr,fixture=fixture,method=method,number=repetition,ledger=ledger,native_wall_ns=native['process_wall_ns'],status='complete',diagnostics=audit['queries'],admission_report=admission,provenance=dict(source_manifest_sha256=source,worker_sha256=p['worker_sha256'],host_id=sha256(p['host']),resource_envelope_id=sha256(p['resources']),calibration_id='valgrind-3.22-amd64-Ir',report_sha256=report_sha256(pr)))
        run['stage_audit'] = audit
        run['peak_rss_bytes'] = native['peak_rss_bytes']
        require(read(trial/'run.json') == run, 'canonical run record changed')
        require(receipt['run_id'] == run['run_id'] and receipt['candidate_id'] == candidate['candidate_id'], 'receipt identity changed')
        rows.append(dict(cell=cell,repetition=repetition,candidate_id=candidate['candidate_id'],run_id=run['run_id'],total_operations=run['total_operations'],native_wall_ns=run['native_wall_ns'],rank=audit['queries']['final_rank'],matrix_sha256=audit['matrix_sha256']))
expected_outcomes = [dict(cell=r['cell'], repetition=r['repetition'], status='VERIFIED',
    promotion_eligible=False, candidate_id=r['candidate_id'], run_id=r['run_id']) for r in rows]
require(s['outcomes'] == expected_outcomes, 'summary mismatch')
print(json.dumps(dict(status='VERIFIED',scope='transported producer integration evidence; no comparison',reference=reference,source_manifest_sha256=source,worker_sha256=p['worker_sha256'],verified_pairs=len(rows),rows=rows),sort_keys=True,indent=2))
