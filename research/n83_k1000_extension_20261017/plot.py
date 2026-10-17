#!/usr/bin/env python3
"""Plot independently replayed N83 K984/K1000/K1016 controls."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
data = json.loads((HERE / 'replay.json').read_text())
assert data['status'] == 'PASS_INDEPENDENT_ARITHMETIC_REPLAY'
rows = sorted(data['rows'], key=lambda row: row['k'])
positions = list(range(len(rows)))
ks = [row['k'] for row in rows]
index = [row['index_build_ms'] / 1000 for row in rows]
rank = [row['rank_stage_ms'] / 1000 for row in rows]
cold = [row['cold_process_ms'] / 1000 for row in rows]
online = [row['target_online_ms'] for row in rows]

plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'n83-k1000-extension-20261017'})
fig, (left, right) = plt.subplots(1, 2, figsize=(9.5, 4.0), dpi=160)
left.bar(positions, index, width=0.6, color='#4c78a8', label='Root index')
left.bar(positions, rank, bottom=index, width=0.6,
         color='#f58518', label='Rank stage')
left.plot(positions, cold, color='#263238', marker='o', linewidth=1.5,
          label='Complete cold process')
left.axhline(75.841, color='#198754', linestyle='--', linewidth=1.3,
             label='Prior K992 repeat (L0)')
left.set(xlabel='Folded columns K', ylabel='Seconds',
         title='Adjacent sizes spend more on rank')
left.set_xticks(positions, ks)
left.set_ylim(0, max(cold) * 1.16)
left.legend(loc='upper right', fontsize=8, frameon=False)
left.grid(axis='y', alpha=0.2)
for position, value in zip(positions, cold):
    left.annotate(f'{value:.1f}', (position, value), xytext=(0, 5),
                  textcoords='offset points', ha='center', fontsize=8)

right.plot(positions, online, color='#b279a2', marker='o', linewidth=1.5)
right.set(xlabel='Folded columns K', ylabel='Target online, ms',
          title='Same exposed point; none beats K992 lookup')
right.set_xticks(positions, ks)
right.set_ylim(0, max(online) * 1.19)
right.grid(axis='y', alpha=0.2)
for position, value in zip(positions, online):
    right.annotate(f'{value:.1f}', (position, value), xytext=(0, 6),
                   textcoords='offset points', ha='center', fontsize=8)

fig.text(0.5, 0.005,
         'L0 shared-host diagnostics; one previously exposed N83 point. '
         'Every base, relation and scalar independently replayed.',
         ha='center', fontsize=8, color='#424242')
fig.tight_layout(rect=(0, 0.055, 1, 1))
for suffix in ('svg', 'png', 'pdf'):
    fig.savefig(HERE / f'screen.{suffix}', metadata={'Date': None})
plt.close(fig)
