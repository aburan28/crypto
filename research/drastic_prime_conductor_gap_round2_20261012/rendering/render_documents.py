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

fig, ax = plt.subplots(figsize=(9, 4.4))
labels, rows = [], []
for p in analysis['pairs']:
    labels.append(f"D=-3 ell={p['ell']}")
    rows.append(p)
for i, p in enumerate(rows):
    for key, delta, color, name in [
        ('same_method_rho_neg', 0.18, '#5e35b1', 'rho.negation both (control)'),
        ('best_method_floor_neg_over_crater_cubic', 0.06, '#00695c', 'best method per curve'),
        ('within_crater_neg_over_cubic', -0.06, '#ef6c00', 'within crater'),
    ]:
        q = p[key]
        if q.get('quotient') is None:
            continue
        lo, hi = q['ci99']
        ax.plot([lo, hi], [i + delta, i + delta], color=color, lw=2)
        ax.plot(q['quotient'], i + delta, 'o', color=color,
                label=name if i == 0 else None)
ax.axvline(1, color='#777', ls='--', lw=1)
ax.axvline(3 ** 0.5, color='#b71c1c', ls=':', lw=1.5,
          label='sqrt(3) = 1.7321 (predicted)')
ax.set_yticks(range(len(labels)), labels)
ax.set_xlabel('Floor / crater mean counted S (99% bootstrap interval)')
ax.set_title('Drastic prime conductor gaps, round 2: the j = 0 crater', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.grid(axis='x', alpha=0.2)
fig.text(0.01, 0.01,
        '32 independent one-target workloads per curve; 20,000 resamples. '
        'Count diagnostics with counted rotations; no wall or hardness claim.',
        fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(root / 'drastic-round2-gap.svg')
plt.close(fig)

# the round-1 lead replication panel
r = analysis['replication_d43_ell1016839']['same_method_rho_neg']
fig, ax = plt.subplots(figsize=(9, 2.2))
if r.get('quotient') is not None:
    lo, hi = r['ci99']
    ax.plot([lo, hi], [0, 0], color='#5e35b1', lw=3)
    ax.plot(r['quotient'], 0, 'o', color='#5e35b1', label='32 fresh targets (round 2)')
    ax.plot([0.528596], [0.15], 'x', color='#b71c1c', label='round 1 (8 targets)')
ax.axvline(1, color='#777', ls='--', lw=1)
ax.set_yticks([])
ax.set_xlabel('floor / crater mean counted S (99% interval)')
ax.set_title('Round-1 lead d43-ell1016839, replicated: the 0.529 reading regressed to one', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.grid(axis='x', alpha=0.2)
fig.tight_layout()
fig.savefig(root / 'replication.svg')
plt.close(fig)


def quote(s):
    return json.dumps(s, ensure_ascii=False)


def make_pdf(directory, title, text, graphs):
    images = ''.join(f'#image({quote(g)}, width: 100%)\n\n' for g in graphs)
    source = f'''#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= {title}

#text({quote(text)})

{images}'''
    (directory / 'report.typ').write_text(source)
    typst.compile(str(directory / 'report.typ'), output=str(directory / 'report.pdf'), root=str(root))


make_pdf(root, 'Drastic prime conductor gaps, round 2: the j = 0 crater',
         (root / 'README.md').read_text(), ['drastic-round2-gap.svg', 'replication.svg'])

for cell in ['d3-ell1016839-rho', 'd3-ell1016839-cubic', 'd3-ell2033609-rho',
             'd3-ell2033609-cubic', 'replicate-d43-ell1016839']:
    rd = root / 'reports' / cell
    recs = [json.loads(l) for l in (root / 'sessions' / cell / 'records.jsonl').read_text().splitlines() if l.strip()]
    fig, ax = plt.subplots(figsize=(9, 2.6))
    by = {}
    for r in recs:
        if r['warmup'] or r['outcome']['status'] != 'verified':
            continue
        by.setdefault(f"{r['arm']} {r['workload']['curve']['slug']}", []).append(r['cost']['s'])
    for y, (k, vs) in enumerate(by.items()):
        ax.plot([min(vs), max(vs)], [y, y], color='#455a64', lw=2)
        ax.plot(sum(vs) / len(vs), y, 'o', color='#00695c')
        ax.text(min(vs), y + 0.15, k, fontsize=7)
    ax.set_yticks(range(len(by)), [''] * len(by))
    ax.set_xlabel('counted S (mean dot, observed range line)')
    ax.set_title(f'{cell}: recorded counts', loc='left', fontsize=9)
    ax.grid(axis='x', alpha=0.2)
    fig.tight_layout()
    fig.savefig(rd / 'counts.svg')
    plt.close(fig)
    make_pdf(rd, f'Run report: {cell}', (rd / 'README.md').read_text(), ['counts.svg'])
print('rendered drastic-round2-gap.svg, replication.svg, report.pdf and run reports')
