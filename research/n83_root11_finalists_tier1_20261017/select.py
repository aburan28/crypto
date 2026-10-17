#!/usr/bin/env python3
"""Frozen two-arm paired selector for 32 fresh N83 public points."""

import argparse
import json
import math
import statistics
from pathlib import Path

SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
TARGETS_SHA256 = '4eb139e73c99f0921041a911d30698997913d1e71155a8e09035ffa647e02a58'
N = 32
POINTS_PER_HOST = 8
T_975_DF31 = 2.0395134463964077


def positive(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value > 0)


def admitted_matrix(data):
    if (data.get('schema') != 'n83.root11-finalists-validated/v1' or
            data.get('source_sha256') != SOURCE_SHA256 or
            data.get('targets_sha256') != TARGETS_SHA256 or
            data.get('evidence_status') != 'PASS_INDEPENDENT_REPLAY' or
            data.get('different_cloud_instances') is not True or
            data.get('same_candidate_bases_across_strata') is not True):
        return None, 'SOURCE_TARGET_OR_REPLAY_GATE'
    ids, choice = data.get('candidate_ids'), data.get('choice')
    if (not isinstance(ids, list) or len(ids) != 2 or len(set(ids)) != 2 or
            not isinstance(choice, dict) or
            choice.get('schema') != 'n83.root11-finalists-choice/v1' or
            choice.get('ordered_candidate_ids') != ids):
        return None, 'CANDIDATE_ID_OR_CHOICE_GATE'
    hosts = data.get('host_strata')
    if (not isinstance(hosts, list) or len(hosts) != 4 or
            {tuple(host.get('point_indices', [])) for host in hosts} !=
            {tuple(range(i * POINTS_PER_HOST, (i + 1) * POINTS_PER_HOST))
             for i in range(4)} or
            any(host.get('tier_status') != 'PASS_TIER1' or
                host.get('aa_status') != 'PASS_10_PERCENT_PER_PAIR'
                for host in hosts)):
        return None, 'HOST_STRATA_GATE'
    points = data.get('points')
    if (not isinstance(points, list) or len(points) != N or
            {row.get('point_index') for row in points} != set(range(N))):
        return None, 'POINT_COVERAGE_GATE'
    matrix = []
    for point in sorted(points, key=lambda row: row['point_index']):
        arms, rho, cal = (point.get(key) for key in
                          ('arms', 'rho_by_candidate', 'calibration_by_candidate'))
        if (not isinstance(arms, list) or len(arms) != 2 or
                not isinstance(rho, dict) or set(rho) != set(ids) or
                not isinstance(cal, dict) or set(cal) != set(ids)):
            return None, 'PAIR_COVERAGE_GATE'
        by_id = {arm.get('candidate_id'): arm for arm in arms}
        if set(by_id) != set(ids) or len(by_id) != 2:
            return None, 'IC_CANDIDATE_ID_GATE'
        if any(arm.get('status') != 'PASS_INDEPENDENT_REPLAY' or
               arm.get('resource_status') != 'PASS' or
               not positive(arm.get('one_target_online_ms')) or
               not positive(arm.get('cold_process_ms'))
               for arm in arms):
            return None, 'IC_REPLAY_RESOURCE_OR_TIMING_GATE'
        if any(item.get('status') != 'PASS_INDEPENDENT_REPLAY' or
               item.get('resource_status') != 'PASS' or
               not positive(item.get('online_ms')) for item in rho.values()):
            return None, 'RHO_REPLAY_RESOURCE_OR_TIMING_GATE'
        if any(item.get('status') != 'PASS_10_PERCENT' for item in cal.values()):
            return None, 'CALIBRATION_GATE'
        matrix.append({'point_index': point['point_index'],
                       'arms': by_id, 'rho': rho})
    return matrix, 'PASS'


def metric_selection(matrix, ids, metric):
    first, second = ids
    if metric == 'one_target_online_ms':
        differences = [math.log(row['arms'][second][metric] /
                                row['arms'][first][metric]) for row in matrix]
    else:
        differences = [row['arms'][second][metric] -
                       row['arms'][first][metric] for row in matrix]
    positive_count = sum(value > 0 for value in differences)
    negative_count = sum(value < 0 for value in differences)
    ties = N - positive_count - negative_count
    dominant = max(positive_count, negative_count)
    p_value = (2 * sum(math.comb(N, k) for k in range(dominant, N + 1))
               / (1 << N)) if ties == 0 else None
    selected = (first if positive_count >= 23 else
                second if negative_count >= 23 else None)
    mean = statistics.mean(differences)
    halfwidth = T_975_DF31 * statistics.stdev(differences) / math.sqrt(N)
    return {
        'status': 'SELECTED' if selected is not None and p_value is not None
                  and p_value <= 0.05 else 'BOUNDED_UNKNOWN',
        'selected_candidate_id': selected if p_value is not None
                                 and p_value <= 0.05 else None,
        'rule': 'two_sided_exact_sign_test_32_non_tied_pairs_alpha_0.05',
        'first_candidate_id': first, 'second_candidate_id': second,
        'positive_second_over_first_count': positive_count,
        'negative_second_over_first_count': negative_count,
        'tie_count': ties, 'two_sided_exact_sign_p': p_value,
        'difference_definition': 'log_second_over_first' if
                                 metric == 'one_target_online_ms' else
                                 'second_minus_first',
        'mean_difference': mean, 'lower_95': mean - halfwidth,
        'upper_95': mean + halfwidth,
        'geometric_second_over_first': math.exp(mean) if
                                       metric == 'one_target_online_ms' else None,
    }


def select(data):
    matrix, admission = admitted_matrix(data)
    if matrix is None:
        return {'schema': 'n83.root11-finalists-selection/v1',
                'admission_status': admission,
                'online': {'status': 'BOUNDED_UNKNOWN',
                           'selected_candidate_id': None},
                'cold': {'status': 'BOUNDED_UNKNOWN',
                         'selected_candidate_id': None}}
    ids = data['candidate_ids']
    tier1 = {}
    for stratum in range(4):
        tier1[str(stratum)] = {}
        for candidate in ids:
            rows = matrix[stratum * POINTS_PER_HOST:(stratum + 1) * POINTS_PER_HOST]
            ratios = [row['rho'][candidate]['online_ms'] /
                      row['arms'][candidate]['one_target_online_ms'] for row in rows]
            tier1[str(stratum)][candidate] = {
                'status': 'PASS_ALL_EIGHT_SAME_POINT' if all(r > 1 for r in ratios)
                          else 'NO_ALL_POINT_SPEEDUP',
                'paired_rho_over_ic_ratios': ratios,
                'geometric_rho_over_ic': math.exp(statistics.mean(
                    math.log(ratio) for ratio in ratios)),
            }
    return {
        'schema': 'n83.root11-finalists-selection/v1',
        'admission_status': admission, 'candidate_ids': ids, 'point_count': N,
        'online': metric_selection(matrix, ids, 'one_target_online_ms'),
        'cold': metric_selection(matrix, ids, 'cold_process_ms'),
        'tier1_speedup_by_host': tier1,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('validated_rows', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = select(json.loads(args.validated_rows.read_text()))
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.output:
        if args.output.exists() and args.output.read_text() != blob:
            raise ValueError('sealed selection output differs')
        args.output.write_text(blob)
    else:
        print(blob, end='')


if __name__ == '__main__':
    main()
