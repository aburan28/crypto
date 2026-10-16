"""Document layout only: plots the native analyzer's already-recorded values."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import typst

root = Path(__file__).resolve().parent.parent
analysis = json.loads((root / 'analysis.json').read_text())
plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
fig, ax = plt.subplots(figsize=(9, 5.2))
labels = []
rows = []
for session in analysis['sessions']:
    rho = session['arms']['rho-neg']['cross_curve']
    if 'ratio_floor_over_crater' not in rho:
        continue
    labels.append(session['cell'])
    rows.append(rho)
for i, row in enumerate(rows):
    lo, hi = row['ci99']
    ax.plot([lo, hi], [i, i], color='#5e35b1', lw=2)
    ax.plot(row['ratio_floor_over_crater'], i, 'o', color='#5e35b1',
            label='rho-neg' if i == 0 else None)
    ic = session = None
ax.axvline(1, color='#777', ls='--', lw=1)
ax.set_yticks(range(len(labels)), labels)
ax.set_xlabel('Floor / crater mean counted S (99% bootstrap interval)')
ax.set_title('Drastic prime conductor gaps: independent target samples on fixed curves', loc='left')
ax.legend(frameon=False)
ax.grid(axis='x', alpha=0.2)
fig.text(0.01, 0.01,
        '8 independent one-target workloads per curve; 20,000 resamples. '
        'Count diagnostics with unpriced native work.',
        fontsize=9)
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig(root / 'drastic-gap.svg')
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


make_pdf(root, 'Drastic prime conductor gaps: ECDLP measured on both ends',
         (root / 'README.md').read_text(), 'drastic-gap.svg')

for session in analysis['sessions']:
    cell = Path(session['session']).name
    rd = root / 'reports' / cell
    fig, ax = plt.subplots(figsize=(9, 3.0))
    y = 0
    ylabels = []
    for arm in ['rho-neg', 'rho-neg-aa', 'ic-mitm']:
        arm_data = session['arms'].get(arm)
        if not arm_data or 'curves' not in arm_data:
            continue
        for r in arm_data['curves']:
            if r.get('mean_s') is None:
                continue
            ax.plot([r['min_s'], r['max_s']], [y, y], color='#455a64', lw=2)
            ax.plot(r['mean_s'], y, 'o', color='#00695c')
            ylabels.append(f"{r['slug']} {arm}")
            y += 1
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=7)
    ax.set_xlabel('counted S (mean dot, observed range line)')
    ax.set_title(f'{cell}: recorded counts', loc='left')
    ax.grid(axis='x', alpha=0.2)
    fig.tight_layout()
    fig.savefig(rd / 'counts.svg')
    plt.close(fig)
    make_pdf(rd, f'Run report: {cell}', (rd / 'README.md').read_text(), 'counts.svg')
print('rendered drastic-gap.svg, report.pdf and run reports')
