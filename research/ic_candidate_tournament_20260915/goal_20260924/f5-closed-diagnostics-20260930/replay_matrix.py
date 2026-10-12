"""Rebuild the closed matrix from witnesses and reject three artifact faults."""
import copy
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
sys.path.insert(0, str(PACKAGE))

from generic_stages import verify_stages
from identity import write_immutable
from oracle import InvalidEvidence, require
from replay_paired_n17_evidence import retained_files


def main(output):
    require(not output.exists(), 'matrix replay output exists; never overwrite evidence')
    result = json.loads((HERE/'RESULT-v2.json').read_text())
    for role, expected in result['analysis_sources'].items():
        require(hashlib.sha256((PACKAGE/role).read_bytes()).hexdigest() == expected,
                'matrix replay source differs from retained analysis: '+role)
    bundle = PACKAGE/'goal_20260924/paired-fresh-n17a1/results-20260929'
    receipt = json.loads((bundle/'receipt.json').read_text())
    require(receipt['archive_sha256'] == result['archive_sha256'], 'closed matrix archive changed')
    files = retained_files(bundle)
    report = json.loads(files['f5/stdout.json'])
    job = json.loads(files['f5/f5-job.json'])
    replay = verify_stages(report, report['fixture'], job)
    require(replay['matrix']['rank'] == 28 and replay['matrix']['certified_logs'] is False,
            'closed rank gap or incomplete log status changed')
    faults = {}
    for name in ('matrix_entry', 'matrix_modulus', 'query_witness'):
        changed = copy.deepcopy(report)
        if name == 'matrix_entry':
            entry = changed['relation_matrix']['rows'][0]['entries'][0]
            entry[1] = str((int(entry[1])+1) % int(changed['relation_matrix']['modulus']))
        elif name == 'matrix_modulus':
            changed['relation_matrix']['modulus'] = '2'
        else:
            attempt = next(row for batch in changed['collection_reports'] for row in batch['attempts']
                           if row['pdp']['points'] is not None)
            attempt['pdp']['points'][0] = (attempt['pdp']['points'][0]+1) % len(changed['factor_base'])
        try:
            verify_stages(changed, changed['fixture'], job)
        except InvalidEvidence as error:
            faults[name] = dict(status='REJECTED', reason=str(error))
        else:
            raise InvalidEvidence('altered closed matrix artifact was admitted: '+name)
    write_immutable(output, dict(schema_version=1, status='AUDITED_CLOSED_MATRIX_REPLAY',
                    archive_sha256=result['archive_sha256'], candidate_id=result['candidate_id'],
                    workload_id=result['workload_id'], run_id=result['run_id'],
                    analysis_receipt_sha256=hashlib.sha256((HERE/'RESULT-v2.json').read_bytes()).hexdigest(),
                    replay_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    input_sha256=result['input_sha256'], auditor_sources=result['analysis_sources'],
                    stages=replay, negative_controls=faults, native_solver_executed=False,
                    complete_solver_admitted=False, promotion_eligible=False, online_speedup=None))


if __name__ == '__main__':
    require(len(sys.argv) == 2, 'usage: replay_matrix.py /absolute/new-matrix-replay.json')
    main(Path(sys.argv[1]))
