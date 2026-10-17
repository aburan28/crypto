#!/usr/bin/env python3
"""Choose two complete high-K arms under the prospective selection rule."""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
K_VALUES = (3584, 4096, 4608, 5120)
STATUS = 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'
PLAN_SHA256 = '89812f798ad2625af7eb1069013abd21fe6d9b68af76bc0e4800dee753f11bb9'
BASELINE_SHA256 = 'b6c3cc7c092fd5a712c79acd199ee863899ea81f1167066cc3f37dc47ab57838'


def choose(summary):
    if (summary['schema'] != 'n83.root11-highk-screen-summary/v1' or
            summary['role'] != 'descriptive_candidate_discovery' or
            summary['hardware_tier'] != 'ordinary_local_exploratory' or
            summary['plan_sha256'] != PLAN_SHA256 or
            summary['baseline_k3072_sha256'] != BASELINE_SHA256 or
            [row['k'] for row in summary['rows']] != list(K_VALUES)):
        raise ValueError('high-K summary identity differs')
    eligible = []
    ended = False
    for row in summary['rows']:
        if row['status'] == STATUS:
            if (ended or row['actual_usable_points'] != 166*row['k'] or
                    row['effective_columns'] != row['k'] or
                    row['geometric_target_ms'] <= 0):
                raise ValueError('noncontiguous or invalid high-K arm')
            eligible.append(row)
        else:
            ended = True
    ordered = sorted(eligible, key=lambda row: (row['geometric_target_ms'], row['k']))
    return {'schema': 'n83.root11-highk-holdout-choice/v1',
            'status': 'READY_FOR_EXACT_MANIFESTS' if len(ordered) >= 2 else
                      'INSUFFICIENT_VERIFIED_HIGHK_ARMS',
            'fixed_k_values': [2816, 3072],
            'selected_highk_values': [row['k'] for row in ordered[:2]] if
                                     len(ordered) >= 2 else [],
            'eligible_highk_values': [row['k'] for row in eligible],
            'candidate_ids': None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('summary', type=Path)
    args = parser.parse_args()
    print(json.dumps(choose(json.loads(args.summary.read_text())), sort_keys=True))


if __name__ == '__main__':
    main()
