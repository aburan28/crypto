#!/usr/bin/env python3
"""Materialize bound high-K tie-break candidates, workloads, and run rows."""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DOMAIN = 'n83-highk-tiebreak-ten-20261017-v1'
SCALAR_RULE = '1+(little_endian(SHA256(domain||u32_le(index))) mod (r-1))'
SOURCE_COMMIT = '416489e09c3d8885b6c12f262935e6ba29373460'


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


runner = module('n83_highk_materialize_run', HERE / 'run.py')
selector = module('n83_highk_materialize_select', HERE / 'select.py')
finegrid = module('n83_highk_materialize_finegrid',
                  ROOT / 'research/n83_e2e_finegrid_20261016/materialize.py')
identity = module('n83_highk_materialize_identity',
                  ROOT / 'research/n83_k3072_current_source_20261017/verify_highk_archives.py')


def read(path):
    return json.loads(path.read_text())


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def write_once(path, value):
    blob = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists() and path.read_text() != blob:
        raise ValueError(f'sealed materialized file differs: {path}')
    path.write_text(blob)


def workload(index, encoded, answer, point_file):
    record = {
        'schema_version': 1, 'curve_id': finegrid.panel.CURVE_ID,
        'target_count': 1, 'public_targets': [encoded],
        'input_law': {
            'domain': DOMAIN, 'fixture_index': index,
            'scalar_rule': SCALAR_RULE,
            'seed_sha256': answer['seed_sha256'],
            'generator_sha256': sha256(HERE / 'generate.py'),
            'answer_sidecar_sha256': sha256(HERE / 'answers.jsonl'),
        },
        'target_sha256': sha256(point_file),
        'cold_or_warm': 'one_public_point_per_fresh_process_empty_base_index_and_logs',
    }
    return {'workload_id': finegrid.panel.digest(record)[:12], **record}


def row(item, candidate_id, workload_id, index, half, artifact_sha, run_id,
        source, host, analysis_path, tier_status, rho, calibration):
    passed = (item['status'] == 'PASS_INDEPENDENT_REPLAY' and
              item.get('resources', {}).get('status') == 'PASS')
    return {
        'candidate_id': candidate_id, 'workload_id': workload_id,
        'run_id': f'{candidate_id}W{workload_id}R1',
        'point_index': index, 'stratum_half': half,
        'status': item['status'],
        'hardware_tier': 'cloud_instance_tier1' if passed and
                         tier_status == 'PASS_TIER1' else
                         'cloud_instance_exploratory',
        'source_snapshot_sha256': source['source_sha256'],
        'run_head_sha': host['git_sha'],
        'executed_binary_sha256': source['binary_sha256'],
        'rho_binary_sha256': source['rho_binary_sha256'],
        'host_receipt_sha256': sha256(host['path']),
        'source_receipt_sha256': sha256(source['path']),
        'correctness_receipt_sha256': sha256(analysis_path),
        'raw_artifact_run_id': run_id,
        'raw_artifact_sha256': artifact_sha,
        'cold_process_ms': item.get('cold_process_ms_exploratory') if passed else None,
        'one_target_online_ms': item.get('one_target_online_ms') if passed else None,
        'paired_rho_online_ms': rho.get('online_ms'),
        'pair_calibration': calibration,
        'online_phases_ms': item.get('online_phases_ms') if passed else None,
        'setup_ms': item.get('setup_ms') if passed else None,
        'rank_stage_ms': item.get('rank_stage_ms') if passed else None,
        'targets_total_ms': item.get('targets_total_ms') if passed else None,
        'rank_attempts': item.get('rank_attempts'),
        'rank_relations': item.get('rank_relations'),
        'target_probes': item.get('target_probes'),
        'actual_usable_points': item.get('actual_factor_base_points'),
        'effective_columns': item.get('effective_columns'),
        'peak_rss_bytes': item.get('resources', {}).get('peak_rss_bytes'),
        'outer_wall_seconds': item.get('full_child_wall_seconds'),
        'raw_evidence': item,
    }


def materialize(validated_path, selection_path, predecessor_validated,
                predecessor_selection, half_dirs):
    validated, selected = read(validated_path), read(selection_path)
    if (validated['schema'] != 'n83.highk-k2543-tiebreak-validated/v1' or
            validated['source_sha256'] != runner.SOURCE_SHA256 or
            validated['targets_sha256'] != runner.TARGETS_SHA256 or
            validated['answers_sha256'] != runner.ANSWERS_SHA256 or
            selected != json.loads(json.dumps(selector.select(validated)))):
        raise ValueError('validated rows or selection differ from frozen rule')
    arms, choice = runner.choose_anchor(predecessor_validated, predecessor_selection)
    if (validated['choice'] != choice or
            validated['candidate_ids'] != [arm[3] for arm in arms] or
            len(half_dirs) != 2 or len(validated['host_strata']) != 2 or
            len(validated['points']) != 10):
        raise ValueError('predecessor choice, candidate IDs or coverage differs')
    target_lines = (HERE / 'targets.jsonl').read_text().splitlines(keepends=True)
    encoded = [json.loads(line) for line in target_lines]
    answers = [json.loads(line) for line in (HERE / 'answers.jsonl').read_text().splitlines()]
    if len(encoded) != 10 or len(answers) != 10:
        raise ValueError('sealed ten-point fixture count differs')
    field, curve = finegrid.panel.curve_record()
    by_id = {arm[3]: arm for arm in arms}
    manifests, rows, rho_rows, workloads = {}, [], [], []
    for half, directory in enumerate(half_dirs):
        analysis_path = directory / 'analysis.json'
        analysis, host, source = (read(directory / name) for name in
                                  ('analysis.json', 'host.json', 'source.json'))
        stratum = validated['host_strata'][half]
        if (analysis['schema'] != 'n83.highk-k2543-tiebreak-half-analysis/v1' or
                analysis['half'] != half or stratum['half'] != half or
                stratum['analysis_sha256'] != sha256(analysis_path) or
                host['instance_id'] != stratum['instance_id'] or
                source['source_sha256'] != runner.SOURCE_SHA256 or
                source['choice'] != choice):
            raise ValueError(f'half {half}: local raw receipts differ')
        host['path'], source['path'] = (directory / 'host.json',
                                       directory / 'source.json')
        for point in analysis['points']:
            index = point['point_index']
            if index not in range(half * 5, half * 5 + 5):
                raise ValueError('point assigned to wrong host stratum')
            raw_point = directory / 'points' / f'p{index:02d}'
            if (raw_point / 'target.jsonl').read_text() != target_lines[index]:
                raise ValueError('per-point public target differs')
            workload_record = workload(index, encoded[index], answers[index],
                                       raw_point / 'target.jsonl')
            workloads.append(workload_record)
            workload_id = workload_record['workload_id']
            for arm in point['arms']:
                candidate_id = arm['candidate_id']
                if candidate_id not in by_id:
                    raise ValueError('unexpected candidate ID')
                name, k, layout, _ = by_id[candidate_id]
                rho = point['rho_by_candidate'][candidate_id]
                calibration = point['calibration_by_candidate'][candidate_id]
                rho_rows.append({
                    'point_index': index, 'paired_candidate_id': candidate_id,
                    'workload_id': workload_id, 'status': rho['status'],
                    'online_ms': rho.get('online_ms'),
                    'calibration': calibration,
                    'host_instance_id': host['instance_id'],
                    'raw_artifact_sha256': stratum['raw_artifact_sha256'],
                    'raw_evidence': rho,
                })
                item = arm['raw_replay']
                if item['status'] == 'PASS_INDEPENDENT_REPLAY':
                    folder = raw_point / 'runs' / name
                    summary = read(folder / 'stdout.jsonl')
                    dump = folder / 'base_dump.jsonl'
                    if candidate_id in (runner.WIDE3[3], runner.ROOT4[3]):
                        kind = ('wide3' if candidate_id == runner.WIDE3[3]
                                else 'root4')
                        manifest = identity.candidate(kind, summary, dump)
                    else:
                        original_commit = finegrid.panel.SOURCE_COMMIT
                        if candidate_id in (runner.prior.K2543[3],
                                            runner.prior.K3072[3]):
                            finegrid.panel.SOURCE_COMMIT = SOURCE_COMMIT
                        try:
                            manifest = finegrid.dynamic.candidate(
                                k, layout, summary, dump, field, curve)
                        finally:
                            finegrid.panel.SOURCE_COMMIT = original_commit
                    if manifest['candidate_id'] != candidate_id:
                        raise ValueError('raw candidate manifest identity differs')
                    if candidate_id in manifests and manifests[candidate_id] != manifest:
                        raise ValueError('candidate manifest differs between points')
                    manifests[candidate_id] = manifest
                rows.append(row(item, candidate_id, workload_id, index, half,
                                stratum['raw_artifact_sha256'], validated['run_id'],
                                source, host, analysis_path, stratum['tier_status'],
                                rho, calibration))
    if (len(rows) != 30 or len(rho_rows) != 30 or len(workloads) != 10 or
            len({item['workload_id'] for item in workloads}) != 10):
        raise ValueError('comparison run or workload count differs')
    if (validated['evidence_status'] == 'PASS_INDEPENDENT_REPLAY' and
            set(manifests) != set(validated['candidate_ids'])):
        raise ValueError('admitted evidence lacks exact candidate manifests')
    output = {
        'schema': 'n83.highk-k2543-tiebreak-ic1/v1',
        'run_id': validated['run_id'],
        'run_attempt': validated['run_attempt'],
        'run_head_sha': validated['head_sha'],
        'validated_rows_sha256': sha256(validated_path),
        'selection_sha256': sha256(selection_path),
        'predecessor_choice': choice,
        'host_strata': validated['host_strata'],
        'workloads': workloads,
        'candidates': sorted(manifests),
        'rows': rows, 'rho_rows': rho_rows,
        'evidence_status': validated['evidence_status'],
        'online_selection': selected['online'],
        'cold_selection': selected['cold'],
    }
    for candidate_id, manifest in manifests.items():
        write_once(HERE / f'candidate_{candidate_id}.json', manifest)
    write_once(HERE / f"ic1_rows_{validated['run_id']}.json", output)
    return {'candidate_count': len(manifests), 'row_count': len(rows),
            'rho_count': len(rho_rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('validated_rows', type=Path)
    parser.add_argument('selection', type=Path)
    parser.add_argument('predecessor_validated', type=Path)
    parser.add_argument('predecessor_selection', type=Path)
    parser.add_argument('half0_results', type=Path)
    parser.add_argument('half1_results', type=Path)
    args = parser.parse_args()
    result = materialize(args.validated_rows.resolve(), args.selection.resolve(),
                         args.predecessor_validated.resolve(),
                         args.predecessor_selection.resolve(),
                         (args.half0_results.resolve(),
                          args.half1_results.resolve()))
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
