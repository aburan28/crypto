#!/usr/bin/env python3
"""Materialize canonical IC1 workload and run records for 24 fresh points."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import check_dispatch
import run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CURVE_ID = 'EC1N83Ckb1a1haacb69ee9f37'
DOMAIN = 'n83-root11-finegrid-holdout24-20261017-v1'
SCALAR_RULE = '1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))'
spec = importlib.util.spec_from_file_location('n83_finegrid_materialize_selector',
                                              HERE / 'select.py')
selector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(selector)


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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(blob)


def materialize(validated_path, selection_path, stratum_dirs):
    dispatch, _ = check_dispatch.check()
    validated = json.loads(validated_path.read_text())
    selection = json.loads(selection_path.read_text())
    if (validated['schema'] != 'n83.root11-finegrid-holdout-validated/v1' or
            validated['source_sha256'] != run.SOURCE_SHA256 or
            validated['targets_sha256'] != run.TARGETS_SHA256 or
            validated['answers_sha256'] != run.ANSWERS_SHA256 or
            selection != selector.select(validated) or
            validated['choice'] != dispatch['choice'] or
            validated['candidate_ids'] != dispatch['choice']['ordered_candidate_ids'] or
            len(validated['points']) != 24 or len(validated['host_strata']) != 4):
        raise ValueError('validated result, dispatch or frozen selection differs')
    public = [json.loads(line) for line in (HERE / 'targets.jsonl').read_text().splitlines()]
    answers = [json.loads(line) for line in (HERE / 'answers.jsonl').read_text().splitlines()]
    if len(public) != 24 or len(answers) != 24:
        raise ValueError('sealed fixture count differs')
    workloads, rows = [], []
    for point in sorted(validated['points'], key=lambda x: x['point_index']):
        index = point['point_index']
        stratum_index = index // 6
        directory = stratum_dirs[stratum_index]
        stratum = validated['host_strata'][stratum_index]
        point_file = directory / 'points' / f'p{index:02d}' / 'target.jsonl'
        target_blob = json.dumps(public[index], separators=(',', ':')) + '\n'
        if (stratum['stratum'] != stratum_index or
                answers[index]['fixture_index'] != index or
                (point_file.exists() and point_file.read_text() != target_blob)):
            raise ValueError('raw point, host or sealed fixture differs')
        workload = {
            'schema_version': 1, 'curve_id': CURVE_ID,
            'target_count': 1, 'public_targets': [public[index]],
            'input_law': {'domain': DOMAIN, 'fixture_index': index,
                          'scalar_rule': SCALAR_RULE,
                          'seed_sha256': answers[index]['seed_sha256'],
                          'generator_sha256': sha256(HERE / 'generate.py'),
                          'answer_sidecar_sha256': sha256(HERE / 'answers.jsonl')},
            'target_sha256': hashlib.sha256(target_blob.encode()).hexdigest(),
            'cold_or_warm': 'one_public_point_per_fresh_process_empty_base_index_and_logs',
        }
        workload_id = hashlib.sha256(canonical(workload)).hexdigest()[:12]
        workloads.append({'workload_id': workload_id, **workload})
        source, host, analysis = (directory / name for name in
                                  ('source.json', 'host.json', 'analysis.json'))
        source_record = json.loads(source.read_text())
        host_record = json.loads(host.read_text())
        if (sha256(analysis) != stratum['analysis_sha256'] or
                source_record['source_sha256'] != validated['source_sha256'] or
                host_record['instance_id'] != stratum['instance_id']):
            raise ValueError('raw source, host or replay receipt differs')
        by_arm = {arm['candidate_id']: arm for arm in point['arms']}
        for candidate in validated['candidate_ids']:
            arm = by_arm.get(candidate, {'status': 'NOT_RUN'})
            raw = arm.get('raw_replay', {})
            rho = point['rho_by_candidate'].get(candidate, {'status': 'NOT_RUN'})
            if candidate not in dispatch['candidate_manifests']:
                raise ValueError('unexpected candidate ID')
            passed = (arm['status'] == 'PASS_INDEPENDENT_REPLAY' and
                      arm.get('resource_status') == 'PASS' and
                      rho['status'] == 'PASS_INDEPENDENT_REPLAY' and
                      rho.get('resource_status') == 'PASS')
            tier_qualified = (passed and
                              validated['evidence_status'] == 'PASS_INDEPENDENT_REPLAY' and
                              validated['different_cloud_instances'] and
                              stratum['tier_status'] == 'PASS_TIER1')
            rows.append({
                'candidate_id': candidate, 'workload_id': workload_id,
                'run_id': f'{candidate}W{workload_id}R1',
                'point_index': index, 'stratum': stratum_index,
                'producer_target_file_status': ('PRESENT_MATCHING' if
                                                point_file.exists() else 'NOT_WRITTEN'),
                'status': arm['status'],
                'rho_status': rho['status'],
                'pair_status': 'PASS_INDEPENDENT_REPLAY' if passed else
                               'CENSORED_OR_FAILED_PAIR',
                'hardware_tier': 'cloud_instance_tier1' if tier_qualified else
                                 'cloud_instance_exploratory',
                'candidate_manifest_sha256':
                    dispatch['candidate_manifests'][candidate]['sha256'],
                'source_snapshot_sha256': source_record['source_sha256'],
                'run_head_sha': host_record['git_sha'],
                'executed_binary_sha256': source_record['binary_sha256'],
                'rho_binary_sha256': source_record['rho_binary_sha256'],
                'host_receipt_sha256': sha256(host),
                'source_receipt_sha256': sha256(source),
                'correctness_receipt_sha256': sha256(analysis),
                'raw_artifact_run_id': validated['run_id'],
                'raw_artifact_sha256': stratum['raw_artifact_sha256'],
                'actual_usable_points': raw.get('actual_factor_base_points'),
                'effective_columns': raw.get('effective_columns'),
                'cold_process_ms': arm.get('cold_process_ms') if passed else None,
                'one_target_online_ms': arm.get('one_target_online_ms') if passed else None,
                'paired_rho_online_ms': rho.get('online_ms') if passed else None,
                'online_phases_ms': raw.get('online_phases_ms') if passed else None,
                'setup_ms': raw.get('setup_ms') if passed else None,
                'rank_stage_ms': raw.get('rank_stage_ms') if passed else None,
                'rank_attempts': raw.get('rank_attempts'),
                'rank_relations': raw.get('rank_relations'),
                'target_probes': raw.get('target_probes'),
                'peak_rss_bytes': raw.get('resources', {}).get('peak_rss_bytes'),
                'pair_calibration': point['calibration_by_candidate'].get(
                    candidate, {'status': 'NOT_RUN'}),
                'raw_evidence': raw,
                'rho_evidence': rho,
            })
    expected = 24 * len(validated['candidate_ids'])
    if (len(workloads) != 24 or len({x['workload_id'] for x in workloads}) != 24 or
            len(rows) != expected or len({x['run_id'] for x in rows}) != expected):
        raise ValueError('canonical IC1 workload or row count differs')
    output = {'schema': 'n83.root11-finegrid-holdout-ic1/v1',
              'provider_run_id': validated['run_id'],
              'run_attempt': validated['run_attempt'],
              'run_head_sha': validated['head_sha'],
              'validated_rows_sha256': sha256(validated_path),
              'selection_sha256': sha256(selection_path),
              'dispatch_sha256': sha256(HERE / 'dispatch.json'),
              'candidate_manifests': dispatch['candidate_manifests'],
              'host_strata': validated['host_strata'],
              'workloads': workloads, 'rows': rows,
              'online_selection': selection['online'],
              'cold_selection': selection['cold']}
    write_once(HERE / 'results' / f'ic1_rows_{validated["run_id"]}.json', output)
    return {'status': 'PASS', 'workloads': len(workloads), 'rows': len(rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('validated_rows', type=Path)
    parser.add_argument('selection', type=Path)
    for stratum in range(4):
        parser.add_argument(f'stratum{stratum}_results', type=Path)
    args = parser.parse_args()
    print(json.dumps(materialize(args.validated_rows.resolve(),
                                  args.selection.resolve(),
                                  tuple(getattr(args, f'stratum{number}_results').resolve()
                                        for number in range(4))), sort_keys=True))


if __name__ == '__main__':
    main()
