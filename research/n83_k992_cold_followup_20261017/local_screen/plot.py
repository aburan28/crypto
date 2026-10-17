#!/usr/bin/env python3
"""Plot the exposed-point local N83 lower-K screen from replayed rows."""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = Path(__file__).resolve().parent
data = json.loads((here / 'replay.json').read_text())
assert data['status'] == 'PASS_INDEPENDENT_ARITHMETIC_REPLAY'
rows = data['rows']
ks = [row['k'] for row in rows]
positions = list(range(len(rows)))
index = [row['index_build_ms'] / 1000 for row in rows]
rank = [row['rank_stage_ms'] / 1000 for row in rows]
cold = [row['cold_process_ms'] / 1000 for row in rows]
online = [row['target_online_ms'] for row in rows]

plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'n83-lower-k-20261017'})
fig, (left, right) = plt.subplots(1, 2, figsize=(10.2, 4.1), dpi=160)
left.bar(positions, index, width=0.6, color='#4c78a8', label='Root index')
left.bar(positions, rank, bottom=index, width=0.6, color='#f58518', label='Rank stage')
left.plot(positions, cold, color='#263238', marker='o', linewidth=1.7,
          label='Complete cold process')
left.set(xlabel='Folded columns K', ylabel='Seconds',
         title='Cold process: index and rank trade off')
left.set_xticks(positions, ks)
left.tick_params(axis='x', rotation=35)
left.legend(loc='upper right', fontsize=8, frameon=False)
left.grid(axis='y', alpha=0.2)

right.plot(positions, online, color='#b279a2', marker='o', linewidth=1.7)
right.set_yscale('log')
right.set(xlabel='Folded columns K', ylabel='Target online, ms (log scale)',
          title='One exposed point: target lookup varies')
right.set_xticks(positions, ks)
right.tick_params(axis='x', rotation=35)
right.grid(axis='y', alpha=0.2)
for position, value in zip(positions, online):
    right.annotate(f'{value:.1f}', (position, value), xytext=(0, 7),
                   textcoords='offset points', ha='center', fontsize=8)

fig.text(0.5, 0.005,
         'L0 shared-host diagnostics; one previously exposed N83 point, '
         'same source and root11 flags. Every scalar replayed.',
         ha='center', fontsize=8, color='#424242')
fig.tight_layout(rect=(0, 0.05, 1, 1))
for suffix in ('svg', 'png', 'pdf'):
    fig.savefig(here / f'screen.{suffix}', metadata={'Date': None})
plt.close(fig)
