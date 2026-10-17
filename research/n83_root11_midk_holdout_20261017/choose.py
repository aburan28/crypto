#!/usr/bin/env python3
"""Prospectively select two root11 capacity-edge challengers."""

import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KS = (3088, 3104, 3114, 3136)
CONTROL_K = 3072
CONTROL_ID = 'IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a'
PASS = 'PASS_INDEPENDENT_BASE_RELATION_SCALAR_REPLAY'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_plan():
    plan = json.loads((HERE / 'selection_plan.json').read_text())
    if (plan['schema'] != 'n83.root11-midk-holdout-selection-plan/v1' or
            plan['screen_k_values'] != list(KS) or
            plan['control_k'] != CONTROL_K or
            plan['control_candidate_id'] != CONTROL_ID or
            plan['choice_rule'] !=
            'two_lowest_complete_geometric_online_if_any_below_k3072_then_lower_k' or
            any(digest(ROOT / relative) != expected
                for relative, expected in plan['sha256'].items())):
        raise ValueError('frozen selection plan differs')
    return plan


def choose(summary):
    plan = frozen_plan()
    if (summary.get('schema') != 'n83.root11-midk-screen-summary/v1' or
            summary.get('plan_sha256') != plan['screen_plan_sha256'] or
            summary.get('baseline_k3072_sha256') != plan['baseline_sha256'] or
            summary.get('hardware_tier') != 'ordinary_local_exploratory' or
            summary.get('role') != 'descriptive_candidate_discovery' or
            len(summary.get('rows', [])) != len(KS) or
            [row.get('k') for row in summary['rows']] != list(KS)):
        raise ValueError('screen summary is outside the frozen rule')
    baseline = json.loads((ROOT / plan['baseline_path']).read_text())
    control = float(baseline['online_geometric_mean_ms'])
    if (digest(ROOT / plan['baseline_path']) != plan['baseline_sha256'] or
            summary['baseline_k3072_geometric_target_ms'] != control or
            not math.isfinite(control) or control <= 0):
        raise ValueError('frozen K3072 control differs')
    eligible = []
    for row in summary['rows']:
        if row['status'] != PASS:
            continue
        geometric = row.get('geometric_target_ms')
        if (not isinstance(geometric, (int, float)) or
                not math.isfinite(geometric) or geometric <= 0 or
                not isinstance(row.get('analysis_sha256'), str) or
                len(row['analysis_sha256']) != 64 or
                row.get('actual_usable_points') != 166 * row['k'] or
                row.get('effective_columns') != row['k']):
            raise ValueError(f"K{row['k']}: incomplete independent replay")
        eligible.append(row)
    eligible.sort(key=lambda row: (row['geometric_target_ms'], row['k']))
    selected = (eligible[:2] if len(eligible) >= 2 and
                eligible[0]['geometric_target_ms'] < control else [])
    return {
        'schema': 'n83.root11-midk-holdout-choice/v1',
        'status': 'SELECTED_K_VALUES_REQUIRES_IC1_MANIFESTS' if selected
                  else 'BOUNDED_UNKNOWN',
        'reason': None if selected else 'INSUFFICIENT_VERIFIED_IMPROVEMENT',
        'control_k': CONTROL_K,
        'control_candidate_id': CONTROL_ID,
        'selected_challenger_k': [row['k'] for row in selected],
        'challenger_candidate_ids': [None for _ in selected],
        'control_geometric_target_ms': control,
        'screen_plan_sha256': plan['screen_plan_sha256'],
        'selection_plan_sha256': digest(HERE / 'selection_plan.json'),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('summary', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = choose(json.loads(args.summary.read_text()))
    blob = json.dumps(result, sort_keys=True, indent=2) + '\n'
    if args.output:
        if args.output.exists() and args.output.read_text() != blob:
            raise ValueError('choice output differs')
        args.output.write_text(blob)
    else:
        print(blob, end='')


if __name__ == '__main__':
    main()
