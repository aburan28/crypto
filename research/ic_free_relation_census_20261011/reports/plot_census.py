"""Plot recorded census values (no computation beyond reading records).
Usage: plot_census.py runs/ OUT.svg"""
import json, sys, glob
import matplotlib; matplotlib.use('svg')
import matplotlib.pyplot as plt
rows = [json.load(open(f)) for f in sorted(glob.glob(sys.argv[1] + '/census_n*.json'))]
fig, ax = plt.subplots(figsize=(8, 4.2))
x = range(len(rows)); w = 0.16
series = [('note: n⁴K⁴/r (predicted)', 'prediction_note_n4K4_over_r', '#a0aec0'),
          ('classes: K⁴n³/r (predicted)', 'prediction_classes_K4n3_over_r', '#90cdf4'),
          ('table collisions → nontrivial rows (measured)', 'nontrivial_relations', '#2b6cb0'),
          ('rank mod r (measured)', 'rank_mod_r', '#c05621')]
for s, (lab, key, col) in enumerate(series):
    vals = [max(r[key], 0.5) for r in rows]
    ax.bar([i + (s - 1.5) * w for i in x], vals, w, label=lab, color=col)
ax.axhline(599, color='k', ls='--', lw=1); ax.text(len(rows) - 0.55, 650, 'rank needed: K − 1 = 599', ha='right', fontsize=8)
ax.set_yscale('log'); ax.set_xticks(list(x)); ax.set_xticklabels([f"n = {r['n']}\nlog₂ r = {r['log2_r']:.1f}" for r in rows])
ax.set_ylabel('count (log scale; 0 drawn at 0.5)'); ax.legend(fontsize=7.5, loc='upper right')
ax.set_title('Free relations in the K = 600 pair-root tables (ic_free_relation_census_20261011 runs/)', fontsize=9)
fig.tight_layout(); fig.savefig(sys.argv[2])
