"""Presentation only: plot recorded native values; never select or calculate statistics."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
import fitz
R=Path(__file__).resolve().parent
S=json.loads((R/'../../selection/selection_freeze.json').read_text())
P=json.loads((R/'native_profiles.json').read_text())
plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(9,4))
for i,c in enumerate(S['candidates']):
    ax.bar(i,c['all_attempt_mean_wall_s'],color='#19747e' if c['eligible'] else '#b75b48',hatch=None if c['eligible'] else '//')
    ax.text(i,c['all_attempt_mean_wall_s']+20,f"{c['all_attempt_mean_wall_s']:.3f}",ha='center',fontsize=9)
ax.set_xticks(range(5),['C0 K600 wide','C1 K600 dense','C2 window600','C3 prefix1271','C4 prefix1798\nFAILED, ineligible'])
ax.set_ylabel('Native mean actual process wall (seconds)')
ax.set_ylim(0,1450);ax.set_title('Training only: two fresh processes per candidate; all attempts retained')
fig.tight_layout();fig.savefig(R/'training_costs.svg');fig.savefig(R/'training_costs.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4.5))
complete=[v for v in P['runs'] if 'producer_phase_ms' in v]
for i,v in enumerate(complete):
    t=v['producer_phase_ms']
    ax.bar(i,t['index_build']/1000,color='#d5a54b',label='Index build' if i==0 else None)
    ax.bar(i,t['rank_stage']/1000,bottom=t['index_build']/1000,color='#19747e',label='Rank stage, including trace overhead' if i==0 else None)
ax.set_xticks(range(len(complete)),[v['id'].replace('train_','t') for v in complete],rotation=30,ha='right')
ax.set_ylabel('Recorded producer time (seconds)');ax.set_title('Preparation components, each completed training process')
ax.legend(fontsize=9);fig.tight_layout();fig.savefig(R/'training_phases.svg');fig.savefig(R/'training_phases.png',dpi=150);plt.close(fig)
styles=getSampleStyleSheet();story=[]
def text(s):story.append(Paragraph(s,styles['BodyText']));story.append(Spacer(1,8))
def title(s):story.append(Paragraph(s,styles['Title']));story.append(Spacer(1,12))
def table(rows,widths=None):
 t=Table(rows,colWidths=widths,repeatRows=1);t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),7),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dcebef')),('GRID',(0,0),(-1,-1),.25,colors.grey)]));story.append(t);story.append(Spacer(1,12))
title('Degree-83 IC basis/setup: training stage')
text('Run started 2026-10-15 15:16 UTC. Published worker/controller/protocol freeze: e07da4efdc092d6d7d1d7b660b29225a665d0775. Selection committed at 19664f02fdcb5f3e5674f5a26e989a117be2f0bc before fresh holdouts were generated.')
text('Observed: all 18 declared processes attempted; 16 completed and two C4 processes failed allocation. Eight completed IC transcripts passed independent C++ rank/log and scalar verification. Three rho A/A work pairs matched. C3 is selected by the predeclared training rule. No fresh-holdout end-to-end winner is established.')
story.append(Image(str(R/'training_costs.png'),width=510,height=227))
text('Means come directly from the native selector. C4 bars show failed-attempt cost and are ineligible; they do not indicate a fast solution. No uncertainty interval was declared for two-point training selection.')
text('Exact legacy curve: a=b=1, GF(2^83) modulus z^83+z^7+z^4+z^2+1, subgroup r=8569786107849059; this is a 53-bit prime subgroup. Same CPU3, one OS thread, 12 GiB RLIMIT_AS and 2100-second complete-process cap for both methods.')
story.append(PageBreak());title('What preparation costs')
story.append(Image(str(R/'training_phases.png'),width=510,height=255))
text('The native profile separates anchor construction, relation search excluding witness checks, witness checks, modular row construction, elimination and backsolve. Producer setup/rank clocks include instrumentation. They are diagnostics inside the actual complete-process wall, not replacements for it.')
rows=[['Run','Probes','Search (s)','Elimination (s)','Backsolve (s)']]
for v in complete:
 t=v['phase_sum_ns'];rows.append([v['id'],str(v['probes']),f"{t['relation_search_excluding_checks']/1e9:.6f}",f"{t['elimination']/1e9:.6f}",f"{t['backsolve']/1e9:.6f}"])
table(rows)
text('Every completed training run obtained one useful rank-gaining relation per attempt: 600 for C0/C1/C2 and 1271 for C3. These are producer counters whose exported modular rows/rank gains and four-point witnesses were independently replayed; search probes and durations were not independently rerun.')
story.append(PageBreak());title('Complete process costs and failures')
rows=[['Process','Status','Wall (s)','CPU (s)','Peak RSS (bytes)']]
for v in S['process_costs']:
 rows.append([v['id'], 'PASS' if v['status']=='COMPLETED_VERIFIED' else 'FAILED',f"{v['wall_s']:.6f}",f"{v['cpu_total_s']:.6f}",str(v['peak_rss_bytes'])])
table(rows)
text('Both C4 stderr records say: memory allocation of 12884901888 bytes failed. Native process receipts retain the allocation failure, return/signal state and raw evidence. The memory-limit-hit indicator remains null; peak RSS is the wait4 whole-process measurement, not requested allocation size.')
text('External checker wall/CPU/RSS is recorded per successful IC run in records/training-audit/audit_meter.json and excluded from solver costs. Both rho training singleton costs and all six A/A process costs remain visible above; no historic IC time is mixed with new rho time.')
story.append(PageBreak());title('Scope, evidence and next stage')
text('C0/C1 use the identical retained K600 mathematical base, isolating wide versus dense layout. C2 uses the next 600 public-scan orbits, verified disjoint from K600. Its entire 1200-orbit construction then discard/reindex is charged. C3/C4 enlarge the public prefix to 1271/1798. Dense packing is a representation change, not a new mathematical base.')
text('An initial draft first-600 candidate equalled K600 and was corrected before timing. Two training points excluded 1911 old/base signed orbits. After the committed selection freeze, eight fresh points were generated outside 1913 orbits. Cold holdouts and three independently seeded eight-target preparation cohorts remain pending.')
text('Public input validation, solve, final generic scalar replay and answer serialization/write/flush are included in matching aligned query clocks; full process wall/CPU/RSS remains primary. Original retained K600 construction and OS-cache flushing are outside scope. Endpoint VmRSS/HWM is not exclusive phase peak; trace timing does not calibrate all instrumentation overhead.')
text('No 1685-target crossover is executed. No calibrated S, L2/independent-host admission, global optimum, asymptotic or deployed-security conclusion is claimed. Training selection is not holdout confirmation. All failed/preflight records remain preserved.')
for label,path in [('Protocol','../../PROTOCOL.md'),('Selection','../../selection/selection_freeze.json'),('Native profiles','native_profiles.json'),('Native profiler source','../../postprocessing/profile_training.cpp'),('Audit receipts','../../records/training-audit/audit_meter.json'),('Raw seals','../../selection/TRAINING_RAW_SHA256SUMS')]:text(f'{label}: {path}')
text('Study PR: https://github.com/aburan28/crypto/pull/1899. Required CI and human review remain merge gates. Report source plots recorded native JSON values; Python computes no selection, statistics or experimental result.')
SimpleDocTemplate(str(R/'n83_basis_training_report.pdf'),pagesize=(612,792),rightMargin=48,leftMargin=48,topMargin=45,bottomMargin=45).build(story)
doc=fitz.open(R/'n83_basis_training_report.pdf')
qa={'pages':len(doc),'all_pages_rendered':True,'empty_pages':[],'out_of_bounds_text_blocks':[]}
for i,p in enumerate(doc):
 if not p.get_text().strip():qa['empty_pages'].append(i+1)
 for b in p.get_text('blocks'):
  if b[0]<0 or b[1]<0 or b[2]>p.rect.width+.1 or b[3]>p.rect.height+.1:qa['out_of_bounds_text_blocks'].append([i+1,list(b[:4])])
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(R/f'qa-page-{i+1}.png')
(R/'PDF_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(qa)
