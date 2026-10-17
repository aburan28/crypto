#!/usr/bin/env python3
"""Apply the prospective four-byte ID30 cloud-dispatch rule."""

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PASS = 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'
ARMS = ('indexed5_id30', 'indexed4_id30')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_plan():
    plan = json.loads((HERE / 'selection_plan.json').read_text())
    if (plan['schema'] != 'n83.indexed4-id30-holdout-selection-plan/v1' or
            plan['arm_names'] != list(ARMS) or plan['k'] != 4096 or
            plan['choice_rule'] != 'dispatch_only_if_both_complete_and_four_byte_geometric_online_lower' or
            any(digest(ROOT / path) != value
                for path, value in plan['sha256'].items())):
        raise ValueError('frozen holdout selection plan differs')
    return plan


def choose(summary, plan):
    if (summary.get('schema') != 'n83.indexed4-id30-screen-summary/v1' or
            summary.get('plan_sha256') != plan['screen_plan_sha256'] or
            summary.get('status') not in ('COMPLETE_BOTH_ARMS', 'PARTIAL_OR_CENSORED') or
            summary.get('hardware_tier') != 'ordinary_local_exploratory' or
            len(summary.get('rows', [])) != 2 or
            [row.get('arm') for row in summary['rows']] != list(ARMS)):
        raise ValueError('screen summary is outside the frozen rule')
    rows = summary['rows']
    complete = all(row.get('status') == PASS for row in rows)
    ratio = None
    if complete:
        for row in rows:
            if (row.get('actual_usable_points') != 679936 or
                    row.get('effective_columns') != 4096 or
                    not isinstance(row.get('analysis_sha256'), str) or
                    len(row['analysis_sha256']) != 64):
                raise ValueError('complete status has an incomplete replay receipt')
        comparison = summary.get('comparison')
        if (summary['status'] != 'COMPLETE_BOTH_ARMS' or
                not isinstance(comparison, dict) or
                comparison.get('paired_points') != 128 or
                comparison.get('ratio_definition') !=
                'five_byte_online_ms_over_four_byte_online_ms'):
            raise ValueError('complete status has an incomplete paired comparison')
        ratio = comparison.get('paired_geometric_ratio')
        if not isinstance(ratio, (int, float)) or not math.isfinite(ratio) or ratio <= 0:
            raise ValueError('invalid paired ratio')
    selected = complete and ratio > 1
    return {
        'schema': 'n83.indexed4-id30-holdout-choice/v1',
        'status': 'SELECTED_REQUIRES_IC1_MANIFESTS' if selected else 'BOUNDED_UNKNOWN',
        'reason': None if selected else 'INCOMPLETE_OR_NO_VERIFIED_LOCAL_GAIN',
        'selected_arm_names': list(ARMS) if selected else [],
        'candidate_ids': [None, None] if selected else [],
        'paired_geometric_ratio': ratio,
        'screen_plan_sha256': plan['screen_plan_sha256'],
        'selection_plan_sha256': digest(HERE / 'selection_plan.json'),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = choose(json.loads(args.summary.read_text()), frozen_plan())
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.output:
        if args.output.exists() and args.output.read_text() != blob:
            raise ValueError('existing choice differs')
        args.output.write_text(blob)
    else:
        print(blob, end='')


if __name__ == '__main__':
    main()
