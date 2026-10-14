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

fig, ax = plt.subplots(figsize=(9, 4.6))
for i, r in enumerate(analysis['results']):
    for key, delta, color, name in [
        ('same_method_paired', 0.15, '#5e35b1', 'same method, paired'),
        ('best_method_paired', 0.0, '#00695c', 'best method, paired'),
        ('within_crater_paired', -0.15, '#ef6c00', 'within crater, paired'),
    ]:
        q = r.get(key)
        if not q or q.get('quotient') is None:
            continue
        lo, hi = q['ci99']
        ax.plot([lo, hi], [i + delta, i + delta], color=color, lw=2)
        ax.plot(q['quotient'], i + delta, 'o', color=color,
                label=name if i == 0 else None)
ax.axvline(1, color='#777', ls='--', lw=1)
ax.axvline(3 ** 0.5, color='#b71c1c', ls=':', lw=1.5, label='sqrt(3)')
ax.axvline(2 ** 0.5, color='#1565c0', ls=':', lw=1.5, label='sqrt(2)')
ax.set_yticks(range(len(analysis['results'])),
             [f"{r['family']} ell={r['ell']}" for r in analysis['results']])
ax.set_xlabel('Floor / crater mean counted S (99% paired bootstrap interval)')
ax.set_title('Drastic prime conductor gap, round 3: paired planted-shared targets', loc='left')
ax.legend(frameon=False, fontsize=8)
ax.grid(axis='x', alpha=0.2)
fig.text(0.01, 0.01,
         'The same planted scalar on both curves; joint bootstrap over indices, 20,000 resamples. '
         'Count diagnostics; no wall or hardness claim.', fontsize=9)
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(root / 'drastic-round3-gap.svg')
plt.close(fig)

def quote(s):
    return json.dumps(s, ensure_ascii=False)

source = f'''#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set heading(numbering: none)
= Drastic prime conductor gap, round 3: paired targets

#text({quote((root / 'README.md').read_text())})

#image({quote('drastic-round3-gap.svg')}, width: 100%)
'''
(root / 'report.typ').write_text(source)
typst.compile(str(root / 'report.typ'), output=str(root / 'report.pdf'), root=str(root))
print('rendered drastic-round3-gap.svg and report.pdf')
