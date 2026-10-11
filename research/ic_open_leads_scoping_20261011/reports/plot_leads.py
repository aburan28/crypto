"""Lay out values from data.json (no new computation). Usage: plot_leads.py data.json out.svg"""
import json, sys
import matplotlib; matplotlib.use('svg')
import matplotlib.pyplot as plt
d = json.load(open(sys.argv[1]))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={'width_ratios': [1.3, 1]})
cells = d['lead1_cells']
a1.barh([c['cell'] for c in cells], [c['columns'] for c in cells], color=['#a0aec0', '#2b6cb0', '#2b6cb0', '#c05621', '#c05621'])
a1.set_xscale('log'); a1.invert_yaxis()
a1.axvline(d['measured_d6_ceiling_columns'], color='k', ls='--', lw=1)
a1.text(d['measured_d6_ceiling_columns'] * 1.2, -0.45, '(9,3) D=7: already\n"large host"', fontsize=7)
for i, c in enumerate(cells):
    a1.text(c['columns'] * 1.3, i, c['decides'], va='center', fontsize=7)
a1.set_xlim(1e6, 1e14)
a1.set_xlabel('Macaulay columns at the deciding degree (derived)')
a1.set_title('Lead 1: price of the m = 4 solving-degree cells', fontsize=9)
pp = d['preprocessing']; x = range(len(pp)); w = 0.25
for k, (key, lab, col) in enumerate([('ic_target_probes', 'IC online target probes (recorded)', '#2b6cb0'),
                                    ('bl_online_steps', 'BL preprocessing rho online steps at equal advice (derived)', '#c05621'),
                                    ('plain_rho_ideal_steps', 'plain rho ideal steps (recorded)', '#a0aec0')]):
    a2.bar([i + (k - 1) * w for i in x], [p[key] for p in pp], w, label=lab, color=col)
a2.set_yscale('log'); a2.set_ylim(1e3, 1e10); a2.set_xticks(list(x)); a2.set_xticklabels([p['rung'] for p in pp])
a2.set_ylabel('operations per target (log)'); a2.legend(fontsize=7, loc='upper left')
a2.set_title('Lead 3A: online cost at equal advice S = root-table entries', fontsize=9)
fig.suptitle('Open IC leads, 2026-10-11 (ic_open_leads_scoping_20261011)', fontsize=10)
fig.tight_layout(); fig.savefig(sys.argv[2])
