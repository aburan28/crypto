"""Lay out recorded D_s values (no computation). Usage: plot_cells.py CELLS.jsonl OUT.svg"""
import json, sys
import matplotlib; matplotlib.use('svg')
import matplotlib.pyplot as plt
cells = [json.loads(l) for l in open(sys.argv[1])]
panels = [('T2 / frob2, n = 21 (subfield-free Frobenius kernels)', lambda c: c['object'] == 'T2' and c['vtype'] == 'frob2' and c['n'] == 21),
          ('T2 / nbw (normal-basis windows)', lambda c: c['object'] == 'T2' and c['vtype'] == 'nbw'),
          ('T2 / poly, n = 23, 25', lambda c: c['object'] == 'T2' and c['vtype'] == 'poly'),
          ('Z / poly, nbw, frob2 (n = 13, 15)', lambda c: c['object'] == 'Z')]
cols = ['#2b6cb0', '#c05621', '#2f855a', '#6b46c1', '#b83280', '#4a5568', '#d69e2e', '#319795']
fig, axes = plt.subplots(1, 4, figsize=(16, 4.4), sharey=True)
for ax, (title, keep) in zip(axes, panels):
    sel = [c for c in cells if keep(c)]
    keys = sorted({(c['vtype'], c['n'], c['a'], c['seed']) for c in sel})
    for i, k in enumerate(keys):
        s = sorted((c for c in sel if (c['vtype'], c['n'], c['a'], c['seed']) == k), key=lambda c: c['ell'])
        s = [c for c in s if c['curve']['d_s'] is not None]
        if not s: continue
        off = (i - (len(keys) - 1) / 2) * 0.05
        ctl = [min([c['r']['d_s']] + ([c['zr']['d_s']] if c.get('zr') else [])) for c in s]
        lab = (f"{k[0]} " if title.startswith('Z') else '') + f"n={k[1]} a={k[2]}" + (f" s={k[3] % 100}" if k[0] == 'frob2' and k[1] == 21 else '')
        ax.plot([c['ell'] + off for c in s], [c['curve']['d_s'] for c in s], '-o', ms=4, color=cols[i % 8], label=lab)
        ax.plot([c['ell'] + off for c in s], ctl, ':x', ms=5, color=cols[i % 8], alpha=0.8)
    ax.set_title(title, fontsize=9); ax.set_xlabel('ℓ = dim V'); ax.grid(alpha=0.3); ax.legend(fontsize=6.5, loc='upper left')
axes[0].set_ylabel('sampled cut-out degree D_s (K = 2000)')
fig.suptitle('Round 2: curve set (solid, o) vs smaller control (dotted, x); censored cells omitted. Run: ic_tree_split_structured_v2_20261011 runs/registered', fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[2])
