"""Lay out recorded D_s values (no computation beyond reading the records).

Usage: plot_cells.py CELLS.jsonl OUT.svg
Plotting only (AGENTS.md: Python never computes a result here).
"""
import json
import sys

import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt

cells = [json.loads(l) for l in open(sys.argv[1])]
arms = [('poly', 'poly: V = span{1, z, …, z^(ℓ−1)}', lambda c: True),
        ('frob', 'frob: V = ker g(σ), n = 21', lambda c: c['n'] == 21),
        ('rand', 'rand: random V, n = 21 (replication)', lambda c: True)]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=True)
colors = ['#2b6cb0', '#c05621', '#2f855a', '#6b46c1', '#b83280', '#4a5568']
for ax, (vt, title, keep) in zip(axes, arms):
    sel = [c for c in cells if c['vtype'] == vt and keep(c)]
    keys = sorted({(c['n'], c['a'], c['seed']) for c in sel})
    for i, k in enumerate(keys):
        s = sorted((c for c in sel if (c['n'], c['a'], c['seed']) == k), key=lambda c: c['ell'])
        ells = [c['ell'] for c in s]
        off = (i - (len(keys) - 1) / 2) * 0.06
        lab = f"n={k[0]} a={k[1]}" + (f" s={k[2] % 100}" if vt == 'frob' else '')
        ax.plot([e + off for e in ells], [c['curve']['d_s'] for c in s], '-o', color=colors[i % 6], ms=5, label=f"{lab}: curve")
        ax.plot([e + off for e in ells], [c['r']['d_s'] for c in s], ':x', color=colors[i % 6], ms=6, alpha=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('ℓ = dim V')
    ax.set_xticks(range(4, 9))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc='upper left')
axes[0].set_ylabel('sampled cut-out degree D_s  (K = 2000 samples)')
fig.suptitle('T2 cut-out degree: curve set (solid, o) vs random set of equal size (dotted, x). '
             'Run: ic_tree_split_structured_v_20261011 runs/registered', fontsize=10)
fig.tight_layout()
fig.savefig(sys.argv[2])
