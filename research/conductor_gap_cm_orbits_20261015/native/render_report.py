"""Render recorded/native-derived findings. No research arithmetic or statistics."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

R=Path(__file__).resolve().parents[1]
load=lambda p: json.loads((R/p).read_text())
orbit=load('results/isolated-summary.json')
primary=load('results/paired-primary.json')
replication=load('results/replication-metrics.json')
labels={0:'Crater CM map',1:'Floor scalar',2:'Crater scalar',3:'Crater scalar A/A',4:'Crater recurrence',5:'Floor recurrence'}
(R/'figures').mkdir(exist_ok=True)
fig,ax=plt.subplots(figsize=(8,4))
ax.bar([labels[v['role']] for v in orbit['aggregate']],[v['mean_cpu_seconds'] for v in orbit['aggregate']],color=['#258491','#c58236','#486999','#7890af','#478b65','#76ad7a'])
ax.set_ylabel('Mean thread CPU seconds / 8,192-point construction')
ax.set_title('100 conductor gaps; 20 shuffled repetitions per arm')
ax.tick_params(axis='x',labelrotation=25)
fig.text(.05,.01,'CPU 0 pinned; local monitor uncontended. Physical-host/frequency isolation unverified.',fontsize=8)
fig.tight_layout(rect=(0,.08,1,1))
fig.savefig(R/'figures/orbit-controls.svg',bbox_inches='tight')
fig.savefig(R/'figures/orbit-controls.png',dpi=150,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(6,5))
ax.scatter([v['crater']['queries'] for v in primary['paired']],[v['floor']['queries'] for v in primary['paired']],s=20,alpha=.75,color='#258491')
ax.plot([0,1e9],[0,1e9],color='#777',linestyle='--',label='Equal query counts')
ax.set(xlim=(0,1e9),ylim=(0,1e9),xlabel='Crater CM queries',ylabel='Floor scalar queries',title='100 completed primary curve pairs')
ax.legend();fig.tight_layout()
fig.savefig(R/'figures/paired-query-counts.svg',bbox_inches='tight')
fig.savefig(R/'figures/paired-query-counts.png',dpi=150,bbox_inches='tight');plt.close(fig)

pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBold','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
styles=getSampleStyleSheet()
for style in styles.byName.values(): style.fontName='DejaVuBold' if style.name in ('Title','Heading1','Heading2') else 'DejaVu'
styles['BodyText'].fontSize=9;styles['BodyText'].leading=13
styles['Title'].fontSize=20;styles['Title'].leading=25
story=[]
def para(text): story.extend([Paragraph(text,styles['BodyText']),Spacer(1,9)])
def heading(text): story.extend([Paragraph(text,styles['Heading1']),Spacer(1,8)])
def table(rows,widths):
    wrapped=[[Paragraph(str(v),styles['BodyText']) for v in row] for row in rows]
    t=Table(wrapped,colWidths=widths,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e5f0ef')),('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,0),.6,colors.grey)]))
    story.extend([t,Spacer(1,12)])
def interval(cell,key):
    v=cell[key];ci=v.get('cluster_percentile_95')
    return f"{v['floor_over_crater_ratio_of_means']:.3f}"+(f" [{ci[0]:.3f}, {ci[1]:.3f}]" if ci else ' (descriptive)')

story.append(Paragraph('Conductor-gap research: completed native campaigns',styles['Title']))
para('15 October 2026 · Exact recoveries, structural certificates and remaining work')
para('No systematic crater/floor mean timing difference is resolved in the 100-pair D = −19 population. A cheap CM endomorphism is concrete, while a generic recurrence gives the floor most of the same orbit-construction benefit. Intrinsic ECDLP hardness remains unestablished.')
table([['Completed experiment','Verified attempts / result'],['100-pair primary campaign','500 / 500: 300 folded IC and 200 same-point rho'],['Replication in 30 sessions','265 / 265: 155 IC and 110 rho'],['Same-crater map/scalar comparisons','145 exact agreements: 100 primary plus 45 replication'],['Independent replay','473 IC base/relation/matrix replays; 797 PARI target recoveries'],['Orbit construction','44,000 trials, including A/A, recurrence, pinned repeats and fresh-process batches']], [235,275])
cpu=primary['floor_over_crater_cpu'];queries=primary['floor_over_crater_queries']
para(f"Primary floor/crater native-backend CPU: <b>{cpu['ratio_of_means']:.3f}</b>, paired 95% interval <b>[{cpu['paired_percentile_95'][0]:.3f}, {cpu['paired_percentile_95'][1]:.3f}]</b>. Query ratio: {queries['ratio_of_means']:.3f} [{queries['paired_percentile_95'][0]:.3f}, {queries['paired_percentile_95'][1]:.3f}]. Ten thousand whole-gap bootstrap resamples; all 100 endpoint pairs included.")
story.append(Image(str(R/'figures/paired-query-counts.png'),width=390,height=325))
para('CPU observations are L0 diagnostics; endpoint strata ran mostly in separate queue waves. Bootstrap intervals do not remove host/wave confounding. Native-backend CPU includes orbit/table setup and native recovery. The measured process also includes parameter derivation and independent replay. Native prices remain incomplete; no complete-price IC/rho ratio or S=0 result is admitted.')

story.append(PageBreak());heading('1. The structural effect and its generic control')
para('The tested pair family has p = 5ℓ² − ℓ + 1, t = ℓ − 2, prime N = 5ℓ² − 2ℓ + 4, D<sub>K</sub> = −19 and endpoint conductors 1 and ℓ. The 100 distinct prime gaps span 33,555,419–35,282,963; completed ECDLP comparisons use 53-bit fields.')
para('On the crater a rational degree-5 CM map acts as [λ], where λ² − λ + 5 = 0 modulo N. The floor order has no degree-5 endomorphism. Nevertheless P[k+2] = P[k+1] − [5]P[k], after one initial [λ] multiplication per seed, generates the same linked orbit with generic arithmetic on either endpoint.')
story.append(Image(str(R/'figures/orbit-controls.png'),width=505,height=258))
para('Pinned means per 8,192-point construction: crater CM 5.464 ms; crater recurrence 6.859 ms; floor recurrence 6.857 ms; direct floor scalar multiplication 42.141 ms. The local monitor reported no contention: 0.34 other CPU seconds during 300.43 wall seconds. Physical-host and frequency isolation are unverified.')
para('Twenty fresh-process batches add 12,000 construction records. Each batch contains all 600 model/arm jobs and uses precomputed CM parameters; it is not a separate process per arm. These batches ran at L0. First observed passes and later shuffled passes are recorded separately. Cold/warm cohorts are not a causal timing comparison.')
para('Four public seed logs and 2,048 levels yield 8,192 factor-base points. Signed pair sums produce equations in four seed logs and the target. Scalar-proportional column parts are rejected; the matrix may pin the target before determining every seed log. Independent replay checks the base, each relation, rank and [d]G = Q.')
para('For 136 cases, independent PARI point-order and Weil-pairing witnesses certify crater-twist invariant factors ℓ and qℓ, versus a cyclic floor twist of order qℓ² (q = 5 or 11). This gives derived rational degree-ℓ kernel-line counts ℓ+1 versus one. These counts are not a complete isogeny enumeration. Both original ECDLP endpoints retain the same prime order N.')

story.append(PageBreak());heading('2. Replication and the current curve')
rows=[['Cohort / gap clusters','Query ratio, floor/crater','Native CPU ratio, floor/crater']]
for cell in replication['comparisons']:
    rows.append([f"{cell['family']} / {cell['curve_gap_clusters']}",interval(cell,'queries'),interval(cell,'native_backend_cpu_seconds')])
table(rows,[170,170,170])
para('Intervals use 10,000 whole-gap cluster resamples, retaining all targets/rounds within each sampled gap. The five-round pilot has only one gap and therefore no curve-population interval. Twenty additional floor vertices completed another 40 IC/rho attempts.')
para('The D = −43 native CPU interval narrowly excludes one in this ten-pair exploratory subset; the excess CPU accompanies extra queries. Intervals are unadjusted for multiple comparisons, and timings are L0. This is a replication lead, not a general conductor-hardness conclusion. The coordinate control transports G and Q but reruns prefix seed selection, so it also changes geometric seed points.')
table([['Pair 46 parameter','Exact value'],['p','5947698030651343'],['N = r; cofactor','5947697996161643; 1'],['t; prime gap ℓ','34489701; 34489703'],['Crater y² = x³ + ax + b','a = 5757116383569871; b = 3560173328417442'],['Floor y² = x³ + ax + b','a = 2087370704375809; b = 1402887937802154']], [175,335])
para('The generic-recurrence pilot recovered both pair-46 targets and reproduced the direct-scalar factor-base hashes and query counts exactly. Across all completed sessions, 797 target recoveries passed independent PARI replay, including prime/Hasse witnesses and rejected scalar mutations. All 473 IC transcripts passed native public-base, relation and rank replay, including corrupted-sign rejection.')
para('Six initial small-control preflight failures exposed loss of public fixtures in measured-child normalization. CurveSpec::for_child fixes that transport. Failed receipts remain preserved. The full Rust library run passed 4,116 tests (117 ignored); current fixture and CM-map checks and pinned Rust 1.98 strict Clippy also pass. Conditional Lean recurrence/orbit/norm lemmas are axiom-free; CM classification and point-map bridges are not claimed as formalized.')

story.append(PageBreak());heading('3. Larger gaps, scope and next steps')
para('The fast codomain constructor matched all 100 direct references; the batched variant also passed forced small-batch comparisons. Independent PARI checks certify 15 new floor models through 91 bits and 20 extra floor vertices. The largest independently certified gap in this snapshot is 17,592,186,219,863 over p = 1547425079971670554692473983.')
para('Under a 6 GiB address-space cap, the 95-bit case reached two of four giant batches before its 1,800-second wall cap; two 96-bit cases aborted on FLINT allocation failures. A separate retry of the same rational-family 96-bit input uses an 11 GiB cap and 3,600 seconds. A completed floor is recorded only after independent certification. Full isogeny point maps and dual-point-map certificates remain open.')
para('The older near-limit candidate, ℓ = 129149576694539 and p = 79228162511798403263235114731, is separate. Its nonzero kernel point needs extension degree 2,226,716,839,561 (twist point degree 4,453,433,679,122). A dense 96-bit basis would need about 26.7 TB per field element. This is a prerequisite barrier for the current rational-point method; other kernel representations are not excluded. Its floor is unconstructed.')
para('The earlier depth-20 degree-2 volcano remains a separate completed export: 21 levels, 1,048,576 explicit chosen-trace models and all intermediate strata. It is a smooth conductor ratio traversed by small steps. The large-prime benchmark samples endpoints and 20 additional floors, not every vertex of every prime-degree volcano. Earlier prefix-IC censored results are preserved and are not converted into solve times.')
for text in [
'Repeat D = −43 on more fresh targets and gaps, with the generic degree-11 recurrence control; preserve failures and use curve-level clustering.',
'Transport the four geometric seed points in coordinate controls, as well as G and Q, to isolate arithmetic from prefix-seed selection.',
'Price all native components before a complete IC/rho work claim. Retain cold setup and all five exclusive online phases, including independent replay; earn isolation for primary wall-time comparisons.',
'Complete and independently certify the expanded 96-bit construction. The original non-rational-kernel candidate needs a compatible kernel representation.',
'A prime-field Semaev/F4/SAT complexity claim requires a separately implemented and independently checked backend. The completed folded MITM experiment makes no such claim.'
]: para('• '+text)
para('Sources: REPORT.md, native source, frozen input/binary manifests, results/paired-primary.json, results/replication-metrics.json, results/isolated-summary.json, replay and PARI certificates, and per-session reports. The original six-page summary and evidence ledger remain in inputs/historical-evidence-summary.zip. Rebuild: python native/render_report.py (rendering only). PR: github.com/aburan28/crypto/pull/1888.')
def footer(canvas,doc):
    canvas.setFont('DejaVu',8);canvas.drawString(40,25,'15 October 2026 · Completed, censored and pending work distinguished');canvas.drawRightString(A4[0]-40,25,str(doc.page))
SimpleDocTemplate(str(R/'Conductor_Gap_Native_Progress.pdf'),pagesize=A4,leftMargin=40,rightMargin=40,topMargin=35,bottomMargin=40).build(story,onFirstPage=footer,onLaterPages=footer)
