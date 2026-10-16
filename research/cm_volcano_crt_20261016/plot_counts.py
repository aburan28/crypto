"""Plot only structural counts already recorded by validate.rs; no computation of results."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
rows = [json.loads(line) for line in (root / 'records/acceptance.jsonl').read_text().splitlines()]
labels = [str(row['ell']) for row in rows]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'path'})
fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout='constrained')
for ax, field, title, color in zip(axes, ['floor_classes', 'crt_primes'],
        ['Enumerated floor classes', 'CRT primes consumed'], ['#416c8d', '#398473']):
    values = [row[field] for row in rows]
    bars = ax.bar(labels, values, color=color, width=0.75)
    ax.bar_label(bars, fontsize=8, padding=3)
    ax.set(title=title, xlabel='Prime degree ell', ylabel='Count')
    ax.margins(y=0.18)
    ax.spines[['top', 'right']].set_visible(False)
fig.suptitle('Exact structural counts from the completed acceptance run', fontsize=12, weight='bold')
fig.savefig(root / 'counts.svg')
svg = root / 'counts.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
plt.close(fig)
