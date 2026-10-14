#!/usr/bin/env python3
"""Describe useful-rank yield beside the independent occupancy controls."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

study = Path(__file__).parent
out = {'schema': 'rootfinding.coverage_efficiency/v1', 'panels': {},
       'inference_limit': 'Target-pinned walks stop at different probe counts; rank and cold-cost ratios are descriptive lower bounds, not fixed-budget yield or a qualified speed comparison.'}
for panel in ('train', 'holdout'):
    session_dir = study / 'sessions' / f'coverage-{panel}-local-20261013'
    session = json.loads((session_dir / 'session.json').read_text())
    null_file = study / f'coverage-{panel}-null.json'
    null = json.loads(null_file.read_text())
    records_file = session_dir / 'records.jsonl'
    rows = [json.loads(line) for line in records_file.read_text().splitlines()]
    by_arm = defaultdict(list)
    for row in rows:
        if not row['warmup'] and row['arm'].startswith('b48-stride-'):
            by_arm[row['arm']].append(row)
    arms = {}
    for arm, arm_rows in sorted(by_arm.items()):
        probes = null['arms'][arm]['lookups']
        rank = sum(row['counters']['matrix_rank'] for row in arm_rows)
        cost = sum(row['cost']['total_gae'] for row in arm_rows)
        rss = max(row['time']['max_rss_kib'] for row in arm_rows)
        arms[arm] = {
            'runs': len(arm_rows), 'verified': sum(row['outcome']['status'] == 'verified' for row in arm_rows),
            'total_independent_rank': rank, 'total_cold_gae_lower_bound': cost,
            'lookups': probes, 'peak_rss_kib': rss,
            'independent_rank_per_1000_gae_lower_bound': rank * 1000 / cost,
            'independent_rank_per_1000_lookups': rank * 1000 / probes,
            'mean_rank_per_peak_rss_mib': (rank / len(arm_rows)) / (rss / 1024),
            'relations_per_1000_lookups': null['arms'][arm]['relations'] * 1000 / probes,
            'unique_pair_keys': null['arms'][arm]['unique_pair_keys'],
            'observed_to_uniform_ratio': null['arms'][arm]['observed_to_uniform_ratio'],
            'target_cluster_ratio_interval_95': null['arms'][arm]['target_cluster_ratio_interval_95'],
        }
    out['panels'][panel] = {'session_id': session['session_id'],
                            'binary_sha256': session['binary_sha256'],
                            'records_sha256': hashlib.sha256(records_file.read_bytes()).hexdigest(),
                            'occupancy_sha256': hashlib.sha256(null_file.read_bytes()).hexdigest(),
                            'arms': arms}

(study / 'coverage-efficiency.json').write_text(json.dumps(out, indent=2) + '\n')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), layout='constrained')
colors = {'train': '#4c78a8', 'holdout': '#f58518'}
for i, panel in enumerate(('train', 'holdout')):
    for j, stride in enumerate((1, 7, 23)):
        arm = out['panels'][panel]['arms'][f'b48-stride-{stride}']
        position = j + (i - 0.5) * 0.25
        ratio = arm['observed_to_uniform_ratio']
        lower, upper = arm['target_cluster_ratio_interval_95']
        axes[0].errorbar(ratio, position, xerr=[[ratio - lower], [upper - ratio]],
                         fmt='o', color=colors[panel], capsize=3, label=panel if j == 0 else None)
        axes[1].scatter(arm['independent_rank_per_1000_gae_lower_bound'], position,
                        color=colors[panel], label=panel if j == 0 else None)
axes[0].axvline(1, color='#444444', linestyle='--', linewidth=1)
axes[0].set_xlabel('Observed / uniform table-hit probability')
axes[1].set_xlabel('Useful rank / 1,000 cold counted GAE')
for axis in axes:
    axis.set_yticks(range(3), ['stride 1', 'stride 7', 'stride 23'])
    axis.invert_yaxis()
    axis.grid(axis='x', alpha=0.2)
axes[0].legend(loc='lower right')
fig.suptitle('B48 coverage train and public holdout | L0 counted-cost scope')
fig.savefig(study / 'coverage-efficiency.svg')
fig.savefig(study / 'coverage-efficiency.png', dpi=160)
plt.close(fig)
svg = study / 'coverage-efficiency.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
print(study / 'coverage-efficiency.json')
