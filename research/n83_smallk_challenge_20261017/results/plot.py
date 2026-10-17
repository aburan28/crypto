#!/usr/bin/env python3
"""Render the admitted five-arm N83 comparison and its point-level table."""

import csv
import json
import statistics
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
validated = json.loads((HERE / 'validated_rows.json').read_text())
selection = json.loads((HERE / 'selection.json').read_text())
assert validated['evidence_status'] == 'PASS_INDEPENDENT_REPLAY'
assert selection['admission_status'] == 'PASS'
assert len(validated['points']) == 10 and len(validated['candidate_ids']) == 5

ids = validated['candidate_ids']
labels = ['K1248', 'K1280', 'K2543', 'K3072', 'K1120']
colors = ['#3d68a4', '#6192bc', '#b34f3c', '#d18b50', '#32a07d']
rows = []
for point in validated['points']:
    for candidate, label in zip(ids, labels):
        arm = next(item for item in point['arms']
                   if item['candidate_id'] == candidate)
        rho = point['rho_by_candidate'][candidate]
        rows.append({
            'point_index': point['point_index'], 'factor_base': label,
            'candidate_id': candidate,
            'actual_usable_points': arm['raw_replay']['actual_factor_base_points'],
            'effective_columns': arm['raw_replay']['effective_columns'],
            'ic_online_ms': arm['one_target_online_ms'],
            'rho_online_ms': rho['online_ms'],
            'paired_rho_over_ic': rho['online_ms'] / arm['one_target_online_ms'],
            'cold_process_ms': arm['cold_process_ms'],
            'tier_status': validated['host_strata'][point['point_index'] // 5]['tier_status'],
            'ic_replay_status': arm['status'], 'rho_replay_status': rho['status'],
        })

with (HERE / 'paired_rows.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                            lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout='constrained')
for position, (candidate, label, color) in enumerate(zip(ids, labels, colors), 1):
    arm_rows = [row for row in rows if row['candidate_id'] == candidate]
    measures = ([row['ic_online_ms'] for row in arm_rows],
                [row['paired_rho_over_ic'] for row in arm_rows],
                [row['cold_process_ms'] / 1000 for row in arm_rows])
    for axis, values in zip(axes, measures):
        axis.scatter([position] * len(values), values, color=color,
                     s=17, alpha=0.75, zorder=3)
        axis.hlines(statistics.median(values), position - 0.2,
                    position + 0.2, color='black', linewidth=1.8)

for axis in axes:
    axis.set_xticks(range(1, 6), labels, rotation=25)
    axis.grid(axis='y', alpha=0.2)
axes[0].set_ylabel('One-target IC online (ms)')
axes[0].set_yscale('log')
axes[0].set_title('Online size: bounded unknown')
axes[1].set_ylabel('Same-point rho / IC online')
axes[1].set_yscale('log')
axes[1].axhline(1, color='#666666', linewidth=0.8, linestyle='--')
axes[1].set_title('Fifty verified pairs')
axes[2].set_ylabel('Complete cold process (s)')
axes[2].set_title('Cold: K1120 selected')
fig.suptitle('N83 root11 factor bases on ten sealed points, two Tier 1 runners',
             fontsize=12)
for suffix in ('svg', 'pdf'):
    fig.savefig(HERE / f'panel.{suffix}')
svg = HERE / 'panel.svg'
svg.write_text(''.join(line.rstrip() + '\n'
                       for line in svg.read_text().splitlines()))
