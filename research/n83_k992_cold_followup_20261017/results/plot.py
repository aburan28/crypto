#!/usr/bin/env python3
"""Render the independently replayed ten-point K992 cloud panel."""

import json
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / 'validated_rows.json').read_text())
selection = json.loads((HERE / 'selection.json').read_text())
assert data['evidence_status'] == 'PASS_INDEPENDENT_REPLAY'
assert selection['admission_status'] == 'PASS'

ids = data['candidate_ids']
labels = ['K992', 'K1024', 'K1107']
colors = ['#2869a6', '#28a07a', '#cc6734']
fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.8), layout='constrained')

for idx, candidate in enumerate(ids):
    rows = [(next(arm for arm in point['arms'] if arm['candidate_id'] == candidate),
             point['rho_by_candidate'][candidate]) for point in data['points']]
    cold = [arm['cold_process_ms'] / 1000 for arm, _ in rows]
    online = [arm['one_target_online_ms'] for arm, _ in rows]
    speedup = [rho['online_ms'] / arm['one_target_online_ms'] for arm, rho in rows]
    for axis, values in zip(axes, (cold, online, speedup)):
        axis.scatter([idx + 1] * len(values), values, color=colors[idx],
                     s=17, alpha=0.72, zorder=3)
        axis.hlines(sorted(values)[4:6][0] / 2 + sorted(values)[4:6][1] / 2,
                    idx + 0.79, idx + 1.21, color='black', linewidth=1.8)

axes[0].set_ylabel('Complete cold process (s)')
axes[1].set_ylabel('One-target IC online (ms)')
axes[2].set_ylabel('Same-point rho / IC online')
axes[1].set_yscale('log')
axes[2].set_yscale('log')
axes[2].axhline(1, color='#666666', linewidth=0.8, linestyle='--')
for axis in axes:
    axis.set_xticks((1, 2, 3), labels)
    axis.grid(axis='y', alpha=0.2)
axes[0].set_title('Cold: K1107 selected')
axes[1].set_title('Online size: bounded unknown')
axes[2].set_title('Thirty verified same-point pairs')
fig.suptitle('N83 root11 factor bases on ten sealed points, two Tier 1 runners',
             fontsize=12)
for suffix in ('svg', 'png', 'pdf'):
    fig.savefig(HERE / f'panel.{suffix}', dpi=180)
svg = HERE / 'panel.svg'
svg.write_text(''.join(line.rstrip() + '\n' for line in svg.read_text().splitlines()))
