#!/usr/bin/env python3
"""Give every verified K992 panel point canonical IC1 workload and run IDs."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
CURVE_ID = 'EC1N83Ckb1a1haacb69ee9f37'
DOMAIN = 'n83-k992-cold-followup-ten-20261017-v1'
SCALAR_RULE = '1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False).encode('utf-8')


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'sealed IC1 row file differs: {path}')
    path.write_text(blob)


def frozen_selection(validated):
    spec = importlib.util.spec_from_file_location('n83_k992_materialize_selector',
                                                  STUDY / 'select.py')
    selector = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(selector)
    return selector.select(validated)


def materialize(half_dirs):
    validated = json.loads((HERE / 'validated_rows.json').read_text())
    selection = json.loads((HERE / 'selection.json').read_text())
    if (validated['evidence_status'] != 'PASS_INDEPENDENT_REPLAY' or
            selection != frozen_selection(validated) or
            selection['admission_status'] != 'PASS' or
            validated['run_id'] != 39023971459 or
            len(validated['points']) != 10 or
            len(validated['host_strata']) != 2):
        raise ValueError('K992 validated panel is not admitted')
    public = [json.loads(line) for line in (STUDY / 'targets.jsonl').read_text().splitlines()]
    answers = [json.loads(line) for line in (STUDY / 'answers.jsonl').read_text().splitlines()]
    if len(public) != 10 or len(answers) != 10:
        raise ValueError('sealed fixture count differs')
    manifest_hashes = {}
    for candidate in validated['candidate_ids']:
        manifest = STUDY / f'candidate_{candidate}.json'
        if json.loads(manifest.read_text())['candidate_id'] != candidate:
            raise ValueError('exact candidate manifest differs')
        manifest_hashes[candidate] = sha256(manifest)
    workloads, rows = [], []
    for point in validated['points']:
        index = point['point_index']
        half = index // 5
        directory = half_dirs[half]
        stratum = validated['host_strata'][half]
        point_file = directory / 'points' / f'p{index:02d}' / 'target.jsonl'
        if (stratum['half'] != half or
                point_file.read_text() != json.dumps(public[index], separators=(',', ':')) + '\n' or
                answers[index]['point'] != public[index]):
            raise ValueError('raw point, host or sealed fixture differs')
        workload = {
            'schema_version': 1, 'curve_id': CURVE_ID,
            'target_count': 1, 'public_targets': [public[index]],
            'input_law': {'domain': DOMAIN, 'fixture_index': index,
                          'scalar_rule': SCALAR_RULE,
                          'seed_sha256': answers[index]['seed_sha256'],
                          'generator_sha256': sha256(STUDY / 'generate.py'),
                          'answer_sidecar_sha256': sha256(STUDY / 'answers.jsonl')},
            'target_sha256': sha256(point_file),
            'cold_or_warm': 'one_public_point_per_fresh_process_empty_base_index_and_logs',
        }
        workload_id = hashlib.sha256(canonical(workload)).hexdigest()[:12]
        workloads.append({'workload_id': workload_id, **workload})
        source = directory / 'source.json'
        host = directory / 'host.json'
        analysis = directory / 'analysis.json'
        source_record = json.loads(source.read_text())
        host_record = json.loads(host.read_text())
        if (sha256(analysis) != stratum['analysis_sha256'] or
                source_record['source_sha256'] != validated['source_sha256'] or
                host_record['instance_id'] != stratum['instance_id']):
            raise ValueError('raw source, host or replay receipt differs')
        for arm in point['arms']:
            candidate = arm['candidate_id']
            raw = arm['raw_replay']
            rho = point['rho_by_candidate'][candidate]
            if (candidate not in manifest_hashes or
                    arm['status'] != 'PASS_INDEPENDENT_REPLAY' or
                    rho['status'] != 'PASS_INDEPENDENT_REPLAY' or
                    raw['status'] != 'PASS_INDEPENDENT_REPLAY' or
                    arm['one_target_online_ms'] != raw['one_target_online_ms']):
                raise ValueError('unverified or unknown candidate row')
            rows.append({
                'candidate_id': candidate, 'workload_id': workload_id,
                'run_id': f'{candidate}W{workload_id}R1',
                'point_index': index, 'stratum_half': half,
                'status': arm['status'], 'hardware_tier': 'cloud_instance_tier1',
                'candidate_manifest_sha256': manifest_hashes[candidate],
                'source_snapshot_sha256': source_record['source_sha256'],
                'run_head_sha': host_record['git_sha'],
                'executed_binary_sha256': source_record['binary_sha256'],
                'rho_binary_sha256': source_record['rho_binary_sha256'],
                'host_receipt_sha256': sha256(host),
                'source_receipt_sha256': sha256(source),
                'correctness_receipt_sha256': sha256(analysis),
                'raw_artifact_run_id': validated['run_id'],
                'raw_artifact_sha256': stratum['raw_artifact_sha256'],
                'actual_usable_points': raw['actual_factor_base_points'],
                'effective_columns': raw['effective_columns'],
                'cold_process_ms': arm['cold_process_ms'],
                'one_target_online_ms': arm['one_target_online_ms'],
                'paired_rho_online_ms': rho['online_ms'],
                'rho_over_ic_online': rho['online_ms'] / arm['one_target_online_ms'],
                'online_phases_ms': raw['online_phases_ms'],
                'setup_ms': raw['setup_ms'],
                'rank_stage_ms': raw['rank_stage_ms'],
                'rank_attempts': raw['rank_attempts'],
                'rank_relations': raw['rank_relations'],
                'target_probes': raw['target_probes'],
                'peak_rss_bytes': raw['resources']['peak_rss_bytes'],
                'pair_calibration': point['calibration_by_candidate'][candidate],
            })
    if (len(workloads) != 10 or len({x['workload_id'] for x in workloads}) != 10 or
            len(rows) != 30 or
            len({x['run_id'] for x in rows}) != 30):
        raise ValueError('canonical IC1 workload or row count differs')
    output = {'schema': 'n83.k992-cold-followup-ic1/v1',
              'provider_run_id': validated['run_id'],
              'run_attempt': validated['run_attempt'],
              'run_head_sha': validated['head_sha'],
              'validated_rows_sha256': sha256(HERE / 'validated_rows.json'),
              'selection_sha256': sha256(HERE / 'selection.json'),
              'candidate_manifests': manifest_hashes,
              'host_strata': validated['host_strata'],
              'workloads': workloads, 'rows': rows,
              'online_selection': selection['online'],
              'cold_selection': selection['cold']}
    write_once(HERE / f'ic1_rows_{validated["run_id"]}.json', output)
    return {'status': 'PASS', 'workloads': len(workloads), 'rows': len(rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('half0_results', type=Path)
    parser.add_argument('half1_results', type=Path)
    args = parser.parse_args()
    print(json.dumps(materialize((args.half0_results.resolve(),
                                  args.half1_results.resolve())), sort_keys=True))


if __name__ == '__main__':
    main()
