#!/usr/bin/env python3
"""Preregistered exact paired selector for twenty N83 holdout points."""

import argparse
import json
import math
import statistics
from pathlib import Path

SOURCE_SHA256 = '8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12'
TARGETS_SHA256 = 'c548def473e11ef04e9350c89bad92540c3bfb215a1138847c7bc006078189de'
T_975_DF19 = 2.093024054408263
T_975_DF4 = 2.7764451051977987


def finite_positive(value):
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) and value > 0)


def sign_flip_p(differences):
    if len(differences) != 20 or not all(math.isfinite(x) for x in differences):
        raise ValueError('twenty finite paired differences required')
    observed = abs(math.fsum(differences))
    total = math.fsum(differences)
    tolerance = 1e-10 * max(1.0, math.fsum(abs(x) for x in differences))
    subset_sums = [0.0]
    for difference in differences:
        subset_sums += [value + difference for value in subset_sums]
    count = sum(abs(2 * value - total) + tolerance >= observed
                for value in subset_sums)
    return count / (1 << 20)


def summary(differences):
    mean = statistics.mean(differences)
    sd = statistics.stdev(differences)
    halfwidth = T_975_DF19 * sd / math.sqrt(20)
    return {'mean': mean, 'median': statistics.median(differences),
            'sample_sd': sd, 'lower_95': mean - halfwidth,
            'upper_95': mean + halfwidth,
            'two_sided_sign_flip_p': sign_flip_p(differences)}


def admitted_matrix(data):
    if (data.get('schema') != 'n83.intermediate-holdout-validated/v1' or
            data.get('source_sha256') != SOURCE_SHA256 or
            data.get('targets_sha256') != TARGETS_SHA256 or
            data.get('evidence_status') != 'PASS_INDEPENDENT_REPLAY'):
        return None, 'SOURCE_TARGET_OR_REPLAY_GATE'
    ids = data.get('candidate_ids')
    choice = data.get('choice')
    if (not isinstance(ids, list) or len(ids) not in (3, 4) or
            len(ids) != len(set(ids)) or
            not isinstance(choice, dict) or
            choice.get('schema') != 'n83.intermediate-holdout-choice/v1' or
            choice.get('ordered_candidate_ids') != ids):
        return None, 'CANDIDATE_ID_OR_CHOICE_GATE'
    hosts = data.get('host_strata')
    if (not isinstance(hosts, list) or len(hosts) != 4 or
            {tuple(host.get('point_indices', [])) for host in hosts} !=
            {tuple(range(5 * i, 5 * i + 5)) for i in range(4)} or
            any(host.get('tier_status') != 'PASS_TIER1' or
                host.get('aa_status') != 'PASS_10_PERCENT_PER_PAIR'
                for host in hosts)):
        return None, 'HOST_STRATA_GATE'
    points = data.get('points')
    if (not isinstance(points, list) or len(points) != 20 or
            {row.get('point_index') for row in points} != set(range(20))):
        return None, 'POINT_COVERAGE_GATE'
    matrix = []
    for point in sorted(points, key=lambda row: row['point_index']):
        arms, rho, cal = (point.get(key) for key in
                          ('arms', 'rho_by_candidate', 'calibration_by_candidate'))
        if (not isinstance(arms, list) or len(arms) != len(ids) or
                not isinstance(rho, dict) or set(rho) != set(ids) or
                not isinstance(cal, dict) or set(cal) != set(ids)):
            return None, 'PAIR_COVERAGE_GATE'
        by_id = {arm.get('candidate_id'): arm for arm in arms}
        if set(by_id) != set(ids) or len(by_id) != len(ids):
            return None, 'IC_CANDIDATE_ID_GATE'
        if any(arm.get('status') != 'PASS_INDEPENDENT_REPLAY' or
               arm.get('resource_status') != 'PASS' or
               not finite_positive(arm.get('one_target_online_ms')) or
               not finite_positive(arm.get('cold_process_ms'))
               for arm in arms):
            return None, 'IC_REPLAY_RESOURCE_OR_TIMING_GATE'
        if any(item.get('status') != 'PASS_INDEPENDENT_REPLAY' or
               item.get('resource_status') != 'PASS' or
               not finite_positive(item.get('online_ms')) for item in rho.values()):
            return None, 'RHO_REPLAY_RESOURCE_OR_TIMING_GATE'
        if any(item.get('status') != 'PASS_10_PERCENT' for item in cal.values()):
            return None, 'CALIBRATION_GATE'
        matrix.append({'point_index': point['point_index'],
                       'arms': by_id, 'rho': rho})
    return matrix, 'PASS'


def metric_selection(matrix, ids, metric):
    comparisons = []
    for i, first in enumerate(ids):
        for second in ids[i + 1:]:
            differences = ([math.log(row['arms'][second][metric] /
                                     row['arms'][first][metric]) for row in matrix]
                           if metric == 'one_target_online_ms' else
                           [row['arms'][second][metric] -
                            row['arms'][first][metric] for row in matrix])
            comparisons.append({'first': first, 'second': second,
                                'difference_definition':
                                'log_second_over_first' if
                                metric == 'one_target_online_ms' else
                                'second_minus_first', **summary(differences)})
    alpha = 0.05 / len(comparisons)
    winners = [candidate for candidate in ids if all(
        (item['mean'] > 0 if item['first'] == candidate else item['mean'] < 0)
        and item['two_sided_sign_flip_p'] <= alpha
        for item in comparisons
        if candidate in (item['first'], item['second']))]
    return {'status': 'SELECTED' if len(winners) == 1 else 'BOUNDED_UNKNOWN',
            'selected_candidate_id': winners[0] if len(winners) == 1 else None,
            'bonferroni_alpha_per_pair': alpha,
            'comparisons': comparisons}


def select(data):
    matrix, admission = admitted_matrix(data)
    if matrix is None:
        return {'schema': 'n83.intermediate-holdout-selection/v1',
                'admission_status': admission,
                'online': {'status': 'BOUNDED_UNKNOWN',
                           'selected_candidate_id': None},
                'cold': {'status': 'BOUNDED_UNKNOWN',
                         'selected_candidate_id': None}}
    ids = data['candidate_ids']
    ratios = {candidate: [row['rho'][candidate]['online_ms'] /
                          row['arms'][candidate]['one_target_online_ms']
                          for row in matrix] for candidate in ids}
    tier1 = {}
    for stratum in range(4):
        tier1[str(stratum)] = {}
        for candidate in ids:
            rows = matrix[stratum * 5:stratum * 5 + 5]
            ic = [row['arms'][candidate]['one_target_online_ms'] for row in rows]
            rho = [row['rho'][candidate]['online_ms'] for row in rows]
            differences = [reference - measured for reference, measured in
                           zip(rho, ic)]
            mean = statistics.mean(differences)
            halfwidth = T_975_DF4 * statistics.stdev(differences) / math.sqrt(5)
            tier1[str(stratum)][candidate] = {
                'status': 'PASS_ALL_FIVE_SAME_POINT' if all(x > 0 for x in differences)
                          else 'NO_ALL_POINT_SPEEDUP',
                'paired_ratios': ratios[candidate][stratum * 5:stratum * 5 + 5],
                'median_ratio': statistics.median(
                    ratios[candidate][stratum * 5:stratum * 5 + 5]),
                'mean_ic_online_ms': statistics.mean(ic),
                'mean_rho_online_ms': statistics.mean(rho),
                'mean_paired_rho_minus_ic_ms': mean,
                'paired_difference_lower_95_ms': mean - halfwidth,
                'paired_difference_upper_95_ms': mean + halfwidth,
            }
    return {'schema': 'n83.intermediate-holdout-selection/v1',
            'admission_status': admission,
            'candidate_ids': ids, 'point_count': len(matrix),
            'online': metric_selection(matrix, ids, 'one_target_online_ms'),
            'cold': metric_selection(matrix, ids, 'cold_process_ms'),
            'rho_over_ic_online_by_candidate': ratios,
            'tier1_speedup_by_host': tier1}


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
