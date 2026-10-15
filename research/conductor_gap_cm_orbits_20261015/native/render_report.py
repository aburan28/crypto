"""Rendering only: consume recorded/native-derived data; no research analysis."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table
from reportlab.lib.styles import getSampleStyleSheet
R=Path(__file__).resolve().parents[1]
s=json.loads((R/'results/orbit-summary.json').read_text())
a=s['aggregate']
labels={0:'Crater CM map',1:'Floor scalar',2:'Crater scalar',3:'Crater scalar A/A'}
fig,ax=plt.subplots(figsize=(8,4))
ax.bar([labels[v['role']] for v in a],[v['mean_cpu_seconds'] for v in a],color=['#258491','#c58236','#486999','#7890af'])
ax.set_ylabel('Mean thread CPU seconds / 8,192-point construction')
ax.set_title('100 conductor gaps; 20 shuffled repetitions per arm')
ax.tick_params(axis='x',labelrotation=12)
ax.text(0,-.32,'L0, concurrent work; descriptive CPU values. Raw source: results/orbits.jsonl',transform=ax.transAxes,fontsize=9)
fig.tight_layout();(R/'figures').mkdir(exist_ok=True)
fig.savefig(R/'figures/orbit-controls.svg',bbox_inches='tight')
fig.savefig(R/'figures/orbit-controls.png',dpi=150,bbox_inches='tight');plt.close(fig)
styles=getSampleStyleSheet();doc=SimpleDocTemplate(str(R/'Conductor_Gap_Native_Progress.pdf'),pagesize=A4,rightMargin=40,leftMargin=40)
story=[Paragraph('Conductor-gap research: native progress',styles['Title']),Paragraph('15 October 2026 — completed observations and pending work',styles['Normal']),Spacer(1,16)]
for text in [
'8,000 construction trials completed across 100 prime gaps. All 2,000 same-crater map/scalar comparisons produced identical complete point/coefficient records. Timing observations are L0 and descriptive.',
'The fast square-root Vélu codomain constructor matched all 100 known direct codomains exactly. Independent PARI checks certified 15 additional floor models, including a 91-bit field with prime conductor gap 17,592,186,219,863, and ten D = −43 controls.',
'The small ecbench adapter control verified all six runs after a public-fixture transport bug was fixed. Both the failed and successful sessions are retained. Independent Rust replay checks factor-base scalars, relations, rank and final target recovery.',
'The 100-pair campaign is running: 300 IC attempts and 200 independent same-point rho attempts. No final solve rate or IC/rho speedup is reported here. Work counters include unpriced components.',
'Fresh-target, query-seed, coordinate-model and additional-floor replication remains pending. The 95/96-bit constructor rows were deferred after shared memory approached its 16 GiB limit. No completed 96-bit paired benchmark or general ECDLP hardness gap is claimed.',
'This PDF renders recorded values only. REPORT.md, results/CONDITIONS.md, the native sources, raw JSONL and certificate receipts define the evidence and accounting limits.'
]: story.extend([Paragraph(text,styles['BodyText']),Spacer(1,10)])
story.insert(5,Image(str(R/'figures/orbit-controls.png'),width=500,height=260))
doc.build(story)
