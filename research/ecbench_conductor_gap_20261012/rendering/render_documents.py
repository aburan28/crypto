"""Document layout only: plots the native analyzer's already-recorded values."""
import json
import os
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import typst

root = Path(__file__).resolve().parent.parent
analysis = json.loads((root / 'analysis.json').read_text())
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(9, 3.4))
for i, session in enumerate(analysis['sessions']):
    for arm, delta, color in [('ic-mitm', 0.12, '#00695c'), ('rho-neg', -0.12, '#5e35b1')]:
        row = session['arms'][arm]['cross_curve']
        if 'ratio_floor_over_crater' not in row:
            continue
        lo, hi = row['ci99']
        ax.plot([lo, hi], [i + delta, i + delta], color=color, lw=2)
        ax.plot(row['ratio_floor_over_crater'], i + delta, 'o', color=color,
                label=arm if i == 0 else None)
ax.axvline(1, color='#777', ls='--', lw=1)
ax.set_yticks(range(len(analysis['sessions'])),
              [f"ell = {s['ell']}" for s in analysis['sessions']])
ax.set_xlabel('Floor / crater mean counted S (99% bootstrap interval)')
ax.set_title('Large prime conductor gaps: independent target samples on fixed curves', loc='left')
ax.legend(frameon=False)
ax.grid(axis='x', alpha=0.2)
fig.text(0.01, 0.01, '8 independent one-target workloads per curve; 20,000 resamples. Count diagnostics with unpriced native work.', fontsize=9)
fig.tight_layout(rect=(0, 0.065, 1, 1))
fig.savefig(root / 'conductor-gap.svg')
plt.close(fig)

def quote(s):
    return json.dumps(s, ensure_ascii=False)

def make_pdf(directory, title, text, graph):
    source = f'''#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= {title}

#text({quote(text)})

#image({quote(graph)}, width: 100%)
'''
    (directory / 'report.typ').write_text(source)
    typst.compile(str(directory / 'report.typ'), output=str(directory / 'report.pdf'), root=str(root))

make_pdf(root, 'Large-prime conductor-gap IC follow-up',
         (root / 'README.md').read_text(), 'conductor-gap.svg')
for session in analysis['sessions']:
    cell = Path(session['session']).name
    rd = root / 'reports' / cell
    rd.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 3.6))
    for arm, delta, color in [('ic-mitm', 0.12, '#00695c'), ('rho-neg', -0.12, '#5e35b1')]:
        for i, row in enumerate(session['arms'][arm]['curves']):
            if row['mean_s'] is None:
                continue
            ax.plot([row['min_s'], row['max_s']], [i + delta, i + delta], color=color, lw=2)
            ax.plot(row['mean_s'], i + delta, 'o', color=color, label=arm if i == 0 else None)
    ax.set_xscale('log')
    ax.set_yticks([0, 1], ['crater', 'floor'])
    ax.set_xlabel('Counted S = group-addition equivalents / sqrt(r); dot = mean, line = observed range')
    ax.set_title(f"{cell}: 8 measured workloads per curve and arm", loc='left')
    ax.legend(frameon=False)
    ax.grid(axis='x', alpha=0.2)
    fig.tight_layout()
    fig.savefig(rd / 'counts.svg')
    plt.close(fig)
    make_pdf(rd, f'Run report: {cell}', (rd / 'README.md').read_text(), 'counts.svg')
