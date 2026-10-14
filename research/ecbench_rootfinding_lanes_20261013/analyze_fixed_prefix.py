#!/usr/bin/env python3
"""Aggregate one frozen equal-prefix control with paired target resampling."""
import hashlib
import json
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source = Path(sys.argv[1])
prefix = Path(sys.argv[2])
receipt = json.loads(source.read_text())
data = receipt['result']
targets = data['targets']
assert len(targets) == 8
levels = (512, 1024, 2048, 4096, 8192)
strides = (1, 7, 23)


def arm(target, stride):
    return next(value for value in target['arms'] if value['stride'] == stride)


def milestone(target, stride, level):
    return next(value for value in arm(target, stride)['milestones'] if value['probes'] == level)


def interval(values):
    ordered = sorted(values)
    return [ordered[499], ordered[19499]]


summary = {'schema': 'rootfinding.fixed_prefix_summary/v1',
           'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
           'targets': len(targets), 'strides': {}, 'paired': {},
           'inference_limit': 'The public targets were reused from the prior coverage holdout; this is an exploratory, correlated-walk component test. Bootstrap intervals are over eight targets and no multiple-look correction is applied.'}
for stride in strides:
    entries = arm(targets[0], stride)['unique_pair_keys']
    probability = arm(targets[0], stride)['uniform_hit_probability']
    assert all(arm(target, stride)['unique_pair_keys'] == entries for target in targets)
    levels_out = {}
    for level in levels:
        observations = [milestone(target, stride, level) for target in targets]
        hits = sum(value['hits'] for value in observations)
        rank = sum(value['independent_rank'] for value in observations)
        cost = sum(value['cold_gae_lower_bound'] for value in observations)
        levels_out[str(level)] = {
            'probes_per_target': level, 'total_probes': level * len(targets),
            'hits': hits, 'uniform_expected_hits': probability * level * len(targets),
            'observed_to_uniform': hits / (probability * level * len(targets)),
            'independent_rank': rank, 'mean_rank_per_target': rank / len(targets),
            'independent_rank_per_1000_probes': rank * 1000 / (level * len(targets)),
            'cold_gae_lower_bound': cost, 'independent_rank_per_1000_gae': rank * 1000 / cost,
            'row_ops': sum(value['row_ops'] for value in observations),
            'group_field_invs_only': sum(value['group_field_invs_only'] for value in observations),
        }
    summary['strides'][str(stride)] = {'unique_pair_keys': entries,
                                       'uniform_hit_probability': probability,
                                       'structural_bytes_floor': arm(targets[0], stride)['structural_bytes_floor'],
                                       'levels': levels_out}

rng = random.Random(20261024)
for stride in (7, 23):
    for level in (4096, 8192):
        differences = []
        efficiency_ratios = []
        for _ in range(20_000):
            draw = rng.choices(targets, k=len(targets))
            base = [milestone(target, 1, level) for target in draw]
            candidate = [milestone(target, stride, level) for target in draw]
            base_rank = sum(value['independent_rank'] for value in base)
            candidate_rank = sum(value['independent_rank'] for value in candidate)
            differences.append((candidate_rank - base_rank) / len(draw))
            base_cost = sum(value['cold_gae_lower_bound'] for value in base)
            candidate_cost = sum(value['cold_gae_lower_bound'] for value in candidate)
            efficiency_ratios.append((candidate_rank / candidate_cost) / (base_rank / base_cost))
        baseline = summary['strides']['1']['levels'][str(level)]
        candidate = summary['strides'][str(stride)]['levels'][str(level)]
        summary['paired'][f'{stride}_vs_1_at_{level}'] = {
            'mean_rank_difference_per_target': candidate['mean_rank_per_target'] - baseline['mean_rank_per_target'],
            'mean_rank_difference_interval_95': interval(differences),
            'rank_per_gae_ratio': candidate['independent_rank_per_1000_gae'] / baseline['independent_rank_per_1000_gae'],
            'rank_per_gae_ratio_interval_95': interval(efficiency_ratios),
            'method': 'paired target-cluster bootstrap, 20000 resamples',
        }

prefix.with_suffix('.json').write_text(json.dumps(summary, indent=2) + '\n')
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), layout='constrained')
colors = {1: '#4c78a8', 7: '#f58518', 23: '#54a24b'}
for stride in strides:
    values = summary['strides'][str(stride)]['levels']
    axes[0].plot(levels, [values[str(level)]['mean_rank_per_target'] for level in levels],
                 marker='o', color=colors[stride], label=f'stride {stride}')
    axes[1].plot(levels, [values[str(level)]['observed_to_uniform'] for level in levels],
                 marker='o', color=colors[stride], label=f'stride {stride}')
axes[0].set_ylabel('Mean independent rank per target')
axes[1].set_ylabel('Observed / exact uniform hit expectation')
axes[1].axhline(1, color='#444444', linestyle='--', linewidth=1)
for axis in axes:
    axis.set_xscale('log', base=2)
    axis.set_xticks(levels, [str(level) for level in levels], rotation=30)
    axis.set_xlabel('Same surviving probe prefix per target')
    axis.grid(alpha=0.2)
axes[0].legend()
fig.suptitle('Eight public targets | same batch-8 walk stream for every B48 base')
fig.savefig(prefix.with_suffix('.svg'))
fig.savefig(prefix.with_suffix('.png'), dpi=160)
plt.close(fig)
svg = prefix.with_suffix('.svg')
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
print(prefix.with_suffix('.json'))
