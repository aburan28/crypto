"""Presentation of recorded native statistics only; no new scientific calculations."""
import json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/n83-report-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
import fitz
R=Path(__file__).resolve().parent
S=json.loads((R/'native_statistics.json').read_text());A=json.loads((R/'independent_receipt.json').read_text())
roles=['control','selected','rho'];labels=['C0 retained K600 wide','C3 prefix1271 dense','Optimized rho'];palette=['#586a84','#19747e','#cc8740']
plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(9,4))
for i,role in enumerate(roles):
 y=S['process_summary'][role]['wall_s']['mean'];ax.bar(i,y,color=palette[i]);ax.text(i,y*1.12,f'{y:.6f} s',ha='center')
ax.set_xticks(range(3),labels);ax.set_yscale('log');ax.set_ylim(1,2500);ax.set_ylabel('Actual complete-process mean wall (seconds, log scale)');ax.set_title('Fresh cold holdouts: eight matched target/preparation blocks');fig.tight_layout();fig.savefig(R/'cold_mean_costs.svg');fig.savefig(R/'cold_mean_costs.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4))
for i,role in enumerate(roles):
 rows=[v for v in S['all_process_rows'] if v['role']==role];ax.plot([v['block']+1 for v in rows],[v['wall_s'] for v in rows],marker='o',label=labels[i],color=palette[i])
ax.set_yscale('log');ax.set_xlabel('Fresh target block');ax.set_ylabel('Actual complete-process wall (seconds, log scale)');ax.set_title('Every fresh cold process; matched rank seed per target');ax.legend();fig.tight_layout();fig.savefig(R/'cold_target_costs.svg');fig.savefig(R/'cold_target_costs.png',dpi=150);plt.close(fig)
styles=getSampleStyleSheet();story=[]
def text(s):story.append(Paragraph(s,styles['BodyText']));story.append(Spacer(1,8))
def title(s):story.append(Paragraph(s,styles['Title']));story.append(Spacer(1,12))
def table(rows):
 t=Table(rows,repeatRows=1);t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTSIZE',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),5),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dcebef')),('GRID',(0,0),(-1,-1),.25,colors.grey)]));story.append(t);story.append(Spacer(1,10))
title('Degree-83 basis/setup: fresh cold holdouts')
text('Stage started 2026-10-15 about18:42UTC. All26 declared processes completed and passed independent verification:26 canonical scalars,16 full IC rank/log audits and one identical rho A/A work pair; zero invalid/missing. This report is the cold stage, not completion of the full study.')
text('Observed complete-process means: selected C3 611.139215195s, K600 control1223.782148706s, optimized rho2.107724219s. All ratios below compare means of actual complete-process costs. Preparation is included.')
story.append(Image(str(R/'cold_mean_costs.png'),width=510,height=227))
rows=[['Native comparison','Mean cost ratio','Conditional95% interval']]
for key,label in [('selected_over_rho','C3/rho'),('control_over_rho','C0/rho'),('selected_over_control','C3/C0')]:
 v=S['comparisons'][key];rows.append([label,f"{v['mean_cost_ratio']:.6f}",f"[{v['conditional_95_percent_interval'][0]:.6f}, {v['conditional_95_percent_interval'][1]:.6f}]"])
table(rows)
text('C3 has lower measured cold cost than K600, but optimized rho has lower cold cost than both IC arms on this panel. For IC/rho, below1 favors IC; C3/C0 is explicitly an IC-basis comparison. No end-to-end IC-over-rho win is established.')
story.append(PageBreak());title('Distribution and all target blocks')
story.append(Image(str(R/'cold_target_costs.png'),width=510,height=227))
rows=[['Role','Mean(s)','Median(s)','p90(s)','p95(s)','Max(s)']]
for role in roles:
 v=S['process_summary'][role]['wall_s'];rows.append([role]+[f"{v[k]:.6f}" for k in ['mean','median','p90','p95','max']])
table(rows)
text('Native analysis uses midpoint conventional medians and type7 linear-interpolated quantiles. The20,000-draw paired percentile bootstrap uses seed202610157801 and resamples all eight matched target/preparation blocks. Intervals are conditional on this panel on one guest; they are not independent-host or population-wide guarantees.')
rows=[['Role','MeanCPU(s)','Mean aligned query(s)','MaxRSS(bytes)']]
for role in roles:
 v=S['process_summary'][role];rows.append([role,f"{v['cpu_s']['mean']:.6f}",f"{v['aligned_online_s']['mean']:.6f}",str(int(v['peak_rss_bytes']['max']))])
table(rows)
text('Aligned queries include generic public-Q validation, solve/recovery, generic scalar replay and answer serialization/write/flush in both methods. They exclude input-file loading and their own timing-receipt emission. Online comparisons remain secondary because actual preparation dominates this boundary. Memory-limit-hit remains null.')
story.append(PageBreak());title('Every measured process, including control')
rows=[['Process','Wall(s)','CPU(s)','PeakRSS(bytes)']]
for v in A['process_costs']:rows.append([v['id'],f"{v['wall_s']:.6f}",f"{v['cpu_total_s']:.6f}",str(v['peak_rss_bytes'])])
table(rows)
text(f"All26 completed; zero stage failures. Main rho retains{S['main_rho_capped_walks_retained']} internal capped walks in all costs. The A/A pair had identical deterministic work; larger/smaller process wall ratio{S['AA_wall_larger_over_smaller']:.6f}. A/A time differences are not forced equal. Checker wall/CPU/RSS is separately metered per IC process and excluded from solver costs.")
story.append(PageBreak());title('Scope, evidence and remaining study')
text('Exact legacy curve: a=b=1, GF(2^83) polynomial z^83+z^7+z^4+z^2+1, prime subgroup r=8569786107849059 (53bits). Common native CPU3 reservation, oneOSthread,12GiB RLIMIT_AS and2100-second process cap. EPYC9V74 KVM guest; physical-host exclusivity and independent-host replication are not established.')
text('C3 prefix1271dense was selected only on training and committed at19664f02fdcb5f3e5674f5a26e989a117be2f0bc before eight fresh holdouts were generated. Frozen remaining37process/100answer matrix and separate controller were published1e93e4ca1 before this stage. Rank seeds101/103/107/109/113/127/131/137 match across the two IC arms within each target block. Method order cycles prospectively.')
text('Training retains two failed C4 allocations. Planned full scope is55processes/118answers; these are not actual verified totals. Eleven already-frozen cohort/A/A processes still follow:three independently seeded8-target shared-setup triplets plus the final A/A pair. Full-study analysis, final addendum, updated comprehensive PDF and evidence ZIP remain pending.')
text('Independent native verification replays each canonical scalar and validates exact target order, seeds/policy/envelope, raw/matrix/binary identities, aligned-clock receipt boundaries and deterministic A/A work. Frozen independent checker verifies every IC anchor/four-point witness/intermediate, full modular row/rank gain, backsolve, representative[log]G and target log reconstruction. Exported probe/time counters are observations, not independently rerun search timings.')
text('Original retained K600 construction and OS-cache flushing are unpriced. Endpoint VmRSS/HWM is not exclusive phase peak. No calibratedS, L2, independent-host, global optimum, asymptotic or deployed-security conclusion. The projected1685-target crossover remains unobserved and is not executed. Historical studies and the prior87-page PDF remain unchanged.')
for path in ['../../PROTOCOL.md','native_statistics.json','independent_receipt.json','../../stage_verification/check_stage.cpp','../../postprocessing/analyze_stage.cpp','../../postprocessing/statistics.hpp'] :text('Evidence/source: '+path)
text('PR https://github.com/aburan28/crypto/pull/1899 remains draft. No ledger promotion or merge; current CI and required human review remain gates. Python presents recorded native values and performs unit/format conversions only.')
SimpleDocTemplate(str(R/'n83_basis_cold_holdout_report.pdf'),pagesize=(612,792),rightMargin=48,leftMargin=48,topMargin=45,bottomMargin=45).build(story)
d=fitz.open(R/'n83_basis_cold_holdout_report.pdf');qa={'pages':len(d),'all_pages_rendered':True,'empty_pages':[],'out_of_bounds_text_blocks':[]}
for i,p in enumerate(d):
 if not p.get_text().strip():qa['empty_pages'].append(i+1)
 for b in p.get_text('blocks'):
  if b[0]<0 or b[1]<0 or b[2]>p.rect.width+.1 or b[3]>p.rect.height+.1:qa['out_of_bounds_text_blocks'].append([i+1,list(b[:4])])
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(R/f'qa-page-{i+1}.png')
(R/'PDF_QA.json').write_text(json.dumps(qa,indent=2)+'\n');print(qa)
