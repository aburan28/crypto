"""Lay out recorded d_F / d_last / D* values (no computation). Usage: plot_fall_degree.py STUDY_DIR OUT.svg"""
import json, sys, glob, os
import matplotlib; matplotlib.use('svg')
import matplotlib.pyplot as plt
D = sys.argv[1]
def rows(pat):
    out = []
    for f in sorted(glob.glob(os.path.join(D, pat))):
        out += [json.loads(l) for l in open(f) if l.strip().startswith('{')]
    return out
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6))
# panel 1: m = 2, d_F against n (Sage TABLE.md n <= 16 medians; native E3 draws)
sage = {8: 3, 10: 3, 12: 3, 14: 3, 16: 3}  # TABLE.md / RESULTS.md E2 Semaev medians (random base)
dreg = {8: 3, 10: 4, 12: 4, 14: 4, 16: 5, 18: 5, 20: 5, 22: 5}
a1.plot(list(sage), list(sage.values()), 's', color='#a0aec0', ms=7, label='Sage median d_F (committed, n ≤ 16)')
e3 = rows('e3/n*.jsonl')
import collections
by = collections.defaultdict(list)
for r in e3:
    by[r['n']].append(r['d_last'])
for n, vs in sorted(by.items()):
    cnt = collections.Counter(vs)
    for v, c in cnt.items():
        a1.scatter([n], [v], s=20 + 12 * c, color='#2b6cb0', zorder=3)
a1.scatter([], [], color='#2b6cb0', label='native E3 draws (marker size = count), rand + sat')
a1.plot(list(dreg), list(dreg.values()), ':', color='#c05621', label='semi-regular D_reg (RESULTS.md)')
a1.set_xlabel('n (m = 2, ℓ = n/2, N = n)'); a1.set_ylabel('degree'); a1.set_ylim(1.5, 7.5); a1.grid(alpha=0.3)
a1.set_title('E3: last fall degree d_F of descended S₃, m = 2', fontsize=9); a1.legend(fontsize=7)
# panel 2: m >= 3, against N: Semaev d_F, recorded plain D*, control d_last
lad = rows('e4/m3_*.jsonl') + rows('e4/x_*.jsonl') + rows('e6/A_*.jsonl')
cols = {3: '#2b6cb0', 4: '#2f855a', 5: '#6b46c1'}
for r in lad:
    a2.scatter(r['n_vars'], r['d_f'], color=cols[r['m']], s=28, zorder=3)
    d = r['d_star'].get('degree')
    if d is not None:
        a2.scatter(r['n_vars'], d, marker='^' if r['d_star'].get('kind') == 'at_least' else 'x', color=cols[r['m']], s=30, alpha=0.6)
for r in rows('e6/B_*.jsonl'):
    v = r['d_f'] if r['d_f'] is not None else r['d_f_at_least']
    a2.scatter(r['n_vars'], v, color=cols[r['m']], s=28, edgecolor='k', zorder=4)
for r in rows('e6/C_*.jsonl'):
    a2.scatter(r['n_vars'], r['d_last'], marker='D', facecolor='none', edgecolor='#2b6cb0', s=30)
ctl = {}
for r in rows('e6/Dp_*.jsonl'):
    ctl[(r['shape'], r['seed'])] = r
for r in ctl.values():
    v = r['d_last'] if r['d_last'] is not None else r['d_last_at_least']
    a2.scatter(r['n_vars'], v, marker='^' if r['d_last'] is None else 's', color='#c05621', s=34)
for m, c in cols.items():
    a2.scatter([], [], color=c, label=f'Semaev m = {m}: d_F (unsat; black edge = new Part B draws)')
a2.scatter([], [], marker='x', color='gray', label='ladder plain D* (recorded; ▲ = lower bound)')
a2.scatter([], [], marker='D', facecolor='none', edgecolor='#2b6cb0', label='Semaev m = 3 sat: d_last')
a2.scatter([], [], marker='s', color='#c05621', label='random control (same shape): d_last (▲ = lower bound)')
a2.set_xlabel('N (Boolean unknowns)'); a2.set_ylim(1.5, 8.5); a2.grid(alpha=0.3); a2.legend(fontsize=6.5, loc='upper left')
a2.set_title('E4/E6: d_F stays at 3–4 while plain D* and controls rise', fontsize=9)
fig.suptitle('Huang–Kosters–Yeo last fall degree of Weil-descended Semaev systems (research/fall_degree_bounds_20261008)', fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[2])
