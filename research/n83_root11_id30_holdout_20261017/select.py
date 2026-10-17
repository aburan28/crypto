#!/usr/bin/env python3
"""Choose two complete high-K arms under the prospective selection rule."""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
K_VALUES = (4096, 4608, 5086)
STATUS = 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'
PLAN_SHA256 = 'ce79d7c6ebadb7a5d38a4d4b9609dabda3e3a224c4a7409862533159f0f07938'
BASELINE_SHA256 = 'b6c3cc7c092fd5a712c79acd199ee863899ea81f1167066cc3f37dc47ab57838'


def choose(summary):
    if (summary['schema'] != 'n83.root11-id30-screen-summary/v1' or
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
    return {'schema': 'n83.root11-id30-holdout-choice/v1',
            'status': 'READY_FOR_EXACT_MANIFESTS' if len(ordered) >= 2 else
                      'INSUFFICIENT_VERIFIED_ID30_ARMS',
            'fixed_k_values': [3072, 3584],
            'selected_id30_values': [row['k'] for row in ordered[:2]] if
                                     len(ordered) >= 2 else [],
            'eligible_id30_values': [row['k'] for row in eligible],
            'candidate_ids': None}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('summary', type=Path)
    args = parser.parse_args()
    print(json.dumps(choose(json.loads(args.summary.read_text())), sort_keys=True))


if __name__ == '__main__':
    main()
