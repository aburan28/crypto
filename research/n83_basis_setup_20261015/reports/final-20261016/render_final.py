"""Presentation only: recorded native values, plots, layout and PDF QA.
No experimental selection, aggregation, scientific ratios or bootstrap in Python.
"""
import os,json,html,shutil
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/n83-final-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitz
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors
R=Path(__file__).resolve().parent; study=R.parent.parent
F=json.loads((R/'native_final.json').read_text())
C=json.loads((study/'reports/run-20261015T1842Z/native_statistics.json').read_text())
T=json.loads((study/'selection/selection_freeze.json').read_text()) if (study/'selection/selection_freeze.json').exists() else json.loads(Path('/workspace/n83_basis_setup/training-fixtures/selection_freeze.json').read_text())
roles=['control','selected','rho'];labels=['C0 K600 wide','C3 prefix1271 dense','Optimized rho'];palette=['#586a84','#147d83','#cc8740']
plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
def save(fig,name):
 fig.tight_layout();fig.savefig(R/(name+'.svg'));svg=R/(name+'.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n');fig.savefig(R/(name+'.png'),dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(9,3.8))
for i,role in enumerate(roles):
 v=C['process_summary'][role]['wall_s']['mean'];ax.bar(i,v,color=palette[i]);ax.text(i,v*1.13,f'{v:.6f} s',ha='center')
ax.set_xticks(range(3),labels);ax.set_yscale('log');ax.set_ylim(1,2500);ax.set_ylabel('Actual complete-process mean (s, log scale)');ax.set_title('Eight fresh cold target/preparation blocks; all costs charged');save(fig,'cold_end_to_end')
fig,ax=plt.subplots(figsize=(9,3.8))
for i,v in enumerate(T['candidates']):
 ax.bar(i,v['all_attempt_mean_wall_s'],color='#147d83' if v['eligible'] else '#aa6670',hatch='' if v['eligible'] else '//');ax.text(i,v['all_attempt_mean_wall_s']+35,f"{v['all_attempt_mean_wall_s']:.3f}",ha='center')
ax.set_xticks(range(5),['C0','C1','C2','C3','C4 failed']);ax.set_ylabel('Mean attempted-process wall (s)');ax.set_ylim(0,1500);ax.set_title('Two training processes per candidate; C4 is ineligible');save(fig,'training_candidates')
fig,ax=plt.subplots(figsize=(9,3.8))
for i,role in enumerate(roles[:2]):
 rows=[x for x in F['all_process_rows'] if x['stage']=='cohort' and x['role']==role];ax.plot(range(1,4),[v['wall_s'] for v in rows],marker='o',label=labels[i],color=palette[i])
ax.set_xticks([1,2,3]);ax.set_xlabel('Independent preparation cluster (rank seed149 /151 /157)');ax.set_ylabel('Actual eight-target process wall (s)');ax.set_title('Cohort IC only: all three rho processes rejected8-point input');ax.legend();save(fig,'cohort_ic_only')
profiles=F['profiles'];picks=[p for p in profiles if p['id'] in ['cohort_0_control','cohort_0_selected']]
fig,axes=plt.subplots(1,2,figsize=(10,3.8))
for p,col in zip(picks,palette):
 rows=p['attempt_rows'];axes[0].plot([v['attempt'] for v in rows],[v['rank_after'] for v in rows],label=p['role'],color=col);axes[1].plot([v['attempt'] for v in rows],[v['memory']['rss_bytes']/1e9 for v in rows],label=p['role'],color=col)
axes[0].set(xlabel='Recorded attempt',ylabel='Reconstructed rank',title='Rank gain: cohort0');axes[1].set(xlabel='Recorded attempt',ylabel='Endpoint RSS (GB)',title='Endpoint samples; not phase peaks');axes[0].legend();axes[1].legend();save(fig,'rank_memory')
fig,ax=plt.subplots(figsize=(10,4.3))
for i,p in enumerate(profiles):
 if p['stage']!='cohort':continue
 x=p['id'];ax.bar(x,p['summary']['timing_ms']['index_build']/1000,color='#6c91a4',label='Index build' if x=='cohort_0_control' else None);ax.bar(x,p['phase_sum_ns']['relation_search_excluding_checks']/1e9,bottom=p['summary']['timing_ms']['index_build']/1000,color='#147d83',label='Relation search' if x=='cohort_0_control' else None)
ax.tick_params(axis='x',labelrotation=25);ax.set_ylabel('Recorded component wall (s)');ax.set_title('Dominant components, not a complete exclusive process partition');ax.legend();save(fig,'preparation_components')
styles=getSampleStyleSheet();story=[];pages=[]
def text(s):story.append(Paragraph(s,styles['BodyText']));story.append(Spacer(1,8))
def page(title):
 if story:story.append(PageBreak())
 pages.append(title);story.append(Paragraph(title,styles['Title']));story.append(Spacer(1,12))
def table(rows,widths=None):
 t=Table(rows,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,1),(-1,-1),'Helvetica'),('FONTSIZE',(0,0),(-1,-1),7.3),('LEADING',(0,0),(-1,-1),9),('BOTTOMPADDING',(0,0),(-1,-1),5),('TOPPADDING',(0,0),(-1,-1),5),('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dcebef')),('GRID',(0,0),(-1,-1),.2,colors.grey)]));story.append(t);story.append(Spacer(1,10))
def graph(name):story.append(Image(str(R/(name+'.png')),width=510,height=215))
def fmt(v):return f'{v:.6f}'
page('Degree-83 IC basis/setup study')
text('Final evidence publication,16 October2026. Exact legacy a=1/tail149 degree83 curve; 53-bit prime subgroup r=8569786107849059. Question: can a better basis and relation strategy reduce preparation enough for a measured end-to-end IC win?')
graph('cold_end_to_end')
text('Observed cold means: C3 611.139215s; retained K600 C0 1223.782149s; optimized rho2.107724s. C3/rho289.952172, conditional95% interval[245.967290,375.396780]. C3/C0 0.499386[0.484206,0.514186]. Optimized rho has lower complete-process cost on all eight matched cold targets. All IC/rho ratios below1 favor IC.')
text('55 attempted processes:50 completed,5 failed;92 independently verified canonical answers (72IC,20rho),30 full rank/log transcripts,5 identical rho A/A work pairs. Planned118 answers were not achieved. Three8-target rho processes rejected the input contract; no cohort IC/rho ratio or winner is available. The two C4 training allocation failures also remain in the evidence.')
text('This addendum closes publication of the attempted bounded campaign and explicitly retains the incomplete cohort comparison. Historical87 pages are preserved in the appended comprehensive report. PR1899; required CI and human review remain merge gates.')
page('Prospective design and basis selection')
table([['Candidate','Mathematical basis','Representation'],['C0','Retained K600','Wide control'],['C1','Identical retained K600','Dense table'],['C2','Next600 orbits after first600','Dense table'],['C3','Public-scan prefix1271','Dense table'],['C4','Public-scan prefix1798','Dense table']])
graph('training_candidates')
text('C0/C1 isolates representation; C1/C2 isolates point selection. C3/C4 enlarge the prefix. Dense packing is engineering, not a new mathematical basis. The first600 draft was exactly K600 and was corrected to a disjoint window before any timing. C2 charges construction of1200 then removal/reindexing of600. No pre-solved logs are supplied.')
text('Training uses two public targets outside1911 excluded signed orbits, rank seeds29/43 matched across bases. Selection is arithmetic mean whole-process wall over both full-rank-verified runs; failed candidates are ineligible and ties lexical. C3 selected at622.6626231075s versus C0 1250.669494426s. Selection committed19664f02f before generating eight holdouts outside1913 excluded signed orbits. Worker/protocol freeze e07da4efd; holdout matrix/controller freeze1e93e4ca1.')
page('Matched cold results and uncertainty')
text('Exact polynomial: z^83+z^7+z^4+z^2+1, a=b=1; G=(6614244651398184860292015,8868260410029249391539529), lambda=2998004823950457, signed Frobenius orbit166. Slug icv1-f2m83-t6151469093347-9eba6cd5. No tail45 or a=0 substitution.')
rows=[['Ratio of native mean costs','Estimate','Conditional95% interval']]
for k,label in [('selected_over_rho','C3/rho'),('control_over_rho','C0/rho'),('selected_over_control','C3/C0')]:
 v=C['comparisons'][k];rows.append([label,fmt(v['mean_cost_ratio']),'['+', '.join(fmt(x) for x in v['conditional_95_percent_interval'])+']'])
table(rows)
text('All26 cold processes/answers passed;16 full IC rank/log audits and one A/A pair passed with zero invalid/missing. Primary paired bootstrap resamples eight matched target/preparation blocks,20,000 draws, seed202610157801. Quantiles use conventional midpoint median and type7 linear interpolation. These are conditional intervals for this panel on one guest, not independent-host guarantees.')
text('Both methods charge actual process launch, loading, preparation, solve, generic scalar verification and output. IC original materialized-K600 selection/construction and OS-cache flushing remain unpriced. C3 constructs its public-scan base in every process. The native CPU3 reservation does not make the physical host exclusive; PSI unavailable and unmovable daemons are logged. No calibrated S or L2 admission is claimed.')
text('Aligned query clocks match: before generic public-Q curve/subgroup validation through solve/recovery, generic scalar replay and answer serialization/write/flush. Both exclude input-file parsing and their own receipt emission. The cold aligned-query mean C3/rho ratio is '+fmt(F['online_cold_selected_over_rho'])+'; this is secondary and does not remove preparation from the end-to-end finding.')
text('No old IC times are paired with new rho times, no gains are multiplied across studies, and no projected1685-target experiment is performed. The old32-target matched cohort and this failed8-target cohort remain separate historical regimes.')
page('Final cohorts: results and failed reference')
graph('cohort_ic_only')
v=F['cohort_selected_over_control'];text('Three genuinely reused8-target IC preparation clusters completed. C0 mean1150.628970s; C3 mean591.192669s. Native C3/C0 '+fmt(v['mean_cost_ratio'])+' ['+', '.join(fmt(x) for x in v['conditional_95_percent_interval'])+'],20,000 paired cluster draws, seed202610157802. Only three preparation clusters limit uncertainty. This IC-basis comparison is not an IC/rho comparison.')
text('All three rho processes returned2 in about0.006s before any answer: inherited worker lines107-108 accept only a singleton or exactly32 points. The predeclared8-point matrix and frozen worker were incompatible. This implementation mismatch should have been caught before timing. No solves are retried; raw records, failed audit and analyzer refusal remain unchanged. Failed-launch wall time is never treated as rho solve time.')
text('Six cohort IC rank/log audits and the final rho A/A pair pass:8/11 processes and50/74 canonical answers, including48IC and2rho. The full-stage analyzer correctly refuses to generate complete cohort statistics. Supplemental native analysis reports the six verified IC processes only and leaves every rho field null. See problems/n83-eight-target-rho-contract.md for the unresolved acceptance conditions.')
page('All cold target aligned-clock rows')
rows=[['Target','C0 query(s)','C3 query(s)','Rho query(s)','C3/rho']]
for v in C['online_target_rows']:rows.append([str(v['target_index']),fmt(v['control']['mean_s']),fmt(v['selected']['mean_s']),fmt(v['rho']['mean_s']),fmt(v['selected_over_rho_mean_ratio'])])
table(rows)
text('Each row is one fresh public target, one complete preparation per IC method, and one independent rho solve. Clock boundaries match in both methods. These are native recorded intervals; native ratios of mean query costs are secondary. Full-process target rows and exact public-Q identities are in the matrix, sidecars and raw answers. Index0..7 follows the frozen holdout target order.')
text('Main cold rho retains253 internal capped walks, all charged. Capped walks that later recover are distinct from failed processes; no slow cases or internal restarts are dropped. All scalar answers are canonical and independently replayed against the exact public Q.')
page('All cohort targets: three recorded repeats')
rows=[['Target / role','Repeat1(s)','Repeat2(s)','Repeat3(s)','Native mean(s)']]
for v in F['cohort_online_target_rows']:
 for role in ['control','selected']:
  x=v[role];rows.append([str(v['target_index'])+' / '+role]+[fmt(n) for n in x['recorded_repeats_s']]+[fmt(x['mean_s'])])
table(rows)
text('Every target has three verified IC query intervals under independent preparations. Rho repeats, online ratios and intervals are unavailable: all three cohort rho workers rejected the8-point input before answering. There is no substitution from cold or historical32-target rho measurements.')
text('Rank seeds149/151/157; target rotations0/3/5, reversal in cohort2 as numbered by protocol. Exact order and per-Q rho seeds are in the frozen matrix. No target-dependent rho DP/walk state sharing was introduced. Repeated-target rows remain secondary to actual complete-process costs.')
for metric,title,unit in [('wall_s','Complete-process wall distributions','seconds'),('cpu_s','Complete-process CPU distributions','seconds'),('peak_rss_bytes','Whole-process peak RSS distributions','bytes'),('aligned_online_s','Aligned-query distributions (secondary)','seconds')]:
 page(title);rows=[['Stage / role','n','Mean','Median','p90','p95','Max']]
 for key,v in F['summaries'].items():
  d=v[metric];formatvalue=(lambda x:str(int(x))) if metric=='peak_rss_bytes' else (lambda x:f'{x:.5f}')
  rows.append([key,str(d['n'])]+[formatvalue(d[k]) for k in ['mean','median','p90','p95','max']])
 table(rows)
 text('Units: '+unit+'. Native arithmetic means, midpoint medians and type7 linear-interpolated p90/p95. Failed C4 and cohort-rho attempts are excluded from successful-cost distributions and retained individually in the attempt tables. A/A controls have their own spread table. No complete cohort-rho distribution exists.')
 text('Whole-process CPU is user+system from wait4; peak RSS is whole-process maximum. Endpoint memory samples below are diagnostic and cannot establish exclusive phase peaks. Memory-limit-hit flag is null. Online distributions pool recorded query intervals as descriptive values; they do not create additional independent preparation clusters.')
for stage in ['training','cold','cohort']:
 rows=[v for v in F['all_process_rows'] if v['stage']==stage]
 for offset in range(0,len(rows),15):
  page(stage.title()+' attempted processes '+str(offset+1)+'-'+str(min(offset+15,len(rows))))
  data=[['Process','Status','Wall(s)','CPU(s)','PeakRSS(bytes)']]
  for v in rows[offset:offset+15]:data.append([v['id'],'complete' if v['completed'] else 'FAILED',fmt(v['wall_s']),fmt(v['cpu_s']),str(v['peak_rss_bytes'])])
  table(data);text('Every declared process appears, including failures, internal-cap costs and controls. Values copied from native final aggregation; execution and raw hashes bind the records. Failed process times are attempted costs, never successful solve observations.')
  if stage=='training':text('Both C4 processes failed allocation of12884901888 bytes under the common12GiB address-space envelope; partial traces and outputs are preserved. C4 is ineligible. This wide-state allocation policy failure does not prove larger bases or alternative packed-state layouts impossible.')
  if stage=='cohort':text('All three rho entries returned2 without answers. Their tiny launch costs cannot be compared with IC solve costs. This stage has no complete matched IC/rho result. The final two singleton rho A/A controls did complete and have identical work.')
page('Where preparation time and memory go')
graph('preparation_components');graph('rank_memory')
text('For all30 completed IC processes, native profiles export every attempt: probes, states, found/useful relation, rank before/after and gain, exclusive anchor/search/check/row/elimination times and endpoint RSS/HWM. Backsolve and prior trace formatting/write/flush overhead are recorded separately. Producer rank/process residuals include uncalibrated instrumentation; these panels are not a complete exclusive cost partition.')
for offset in range(0,len(profiles),10):
 page('Native rank/setup profiles '+str(offset+1)+'-'+str(min(offset+10,len(profiles))))
 rows=[['Process','Attempts / gain','Probes','Index(s)','Search(s)','Elim(s)','Backsolve(s)']]
 for p in profiles[offset:offset+10]:
  ph=p['phase_sum_ns'];rows.append([p['id'],str(p['attempts'])+' / '+str(p['rank_gain']),str(p['probes']),f"{p['summary']['timing_ms']['index_build']/1000:.3f}",f"{ph['relation_search_excluding_checks']/1e9:.3f}",f"{ph['elimination']/1e9:.4f}",f"{ph['backsolve']/1e9:.6f}"])
 table(rows)
 text('Index and relation search dominate these recorded preparation profiles; elimination and backsolve are much smaller. Every useful relation/rank gain is independently reconstructed. Probe/state/time counters are producer observations, not independently repeated search timings. Per-attempt details, all other exclusive phases and memory curves are in native_final.json/profiles and immutable rank.jsonl.')
 rows=[['Process','State layout / slot','State bytes','Root-table bytes']]
 for p in profiles[offset:offset+10]:
  s=p['summary'];rows.append([p['id'],s['root_state_layout']+' / '+str(s['root_state_slot_bytes']),str(s['root_state_allocated_bytes']),str(s['root_table_allocated_bytes'])])
 table(rows)
 text('C3 uses packed_u128 root table slots16 bytes with4GiB table allocation and wide_u128_roots states with6GiB allocation. Retained memory policy is recorded, not a claim of optimal packing. Original K600 construction, OS-cache flushing and instrumentation calibration remain outside the priced boundary.')
page('Independent audits and rho A/A controls')
rows=[['Pair / stage','First wall(s)','Second wall(s)','Larger / smaller','Work']]
for v in F['AA_spreads']:rows.append([v['pair']+' / '+v['stage'],fmt(v['first_wall_s']),fmt(v['second_wall_s']),fmt(v['larger_over_smaller']),'identical'])
table(rows)
text('Native recursion compares every deterministic A/A outcome field; clocks need not match. Five pairs span training, cold and cohort phases. Independent rank checker costs are excluded from solver process costs and separately wait4-metered; per-check wall/user/system CPU/peakRSS and logs are in training-audit, cold-audit and cohort-audit. No audit was repeated to replace a failed result.')
text('Full audits independently reconstruct every modular relation row and rank gain, backsolve all base logs, check every representative[log]G, anchored target/four-point witness/intermediate and final answer/log reconstruction. Training selection independently checked all16 canonical answers; cold26 and cohort50 were independently checked. Thirty transcripts cover72 IC answers, and all20 available rho answers pass canonical scalar replay.')
text('The cohort checker initially reports “Cannot read input file” for missing answer files. Raw process receipts and stderr establish why those files are absent. That diagnostic is retained unchanged; no checker correction is needed. The full-stage statistical gate refuses incomplete evidence as intended. The new final aggregator binds all attempts, raw hashes and previous audit/selection seals, recursively checks all five A/A pairs, and emits null for the missing batch comparison.')
page('Scope, acceptance and unresolved obligations')
table([['Requirement','Actual evidence / status'],['Basis diversity and training-only selection','Verified C0..C3; C4 failed allocations retained'],['Fresh matched cold complete costs','Verified8 triplets; rho lower cost'],['Rank/log export and independent audit','Verified30 full transcripts /72IC answers'],['Three8-target matched cohorts','PARTIAL: sixIC and AA passed; three rho failed'],['Full actual accounting','55 attempts /50complete /5failed /92 answers'],['Unqualified IC-over-rho winner','Not established; cold ratio is above1'],['New solve / projected1685-target run','Not performed'],['Independent host, calibrated S, L2','Unset'],['Current CI and required human review','Merge gates, not scientific admission']],[195,315])
text('One guest, three shared preparation clusters, physical host not exclusive. Isolation records and unavailable PSI/unmovable daemons are preserved. No population-wide, global-optimum, exponent, asymptotic or deployed-security conclusion follows from this one field and subgroup. Larger-base feasibility under other layouts remains open.')
text('The missing batch reference is an implementation/workload-contract gap. A separately versioned prospective correction would need positive8-target native acceptance and duplicate/malformed tests before freezing and timing. No replacement campaign is run in this task. This publication preserves the attempted study; it does not mark the failed matched-cohort acceptance requirement complete.')
page('Source and evidence map')
for s in ['PROTOCOL.md — immutable prospective parameters and acceptance','provenance/SOURCE_SHA256SUMS and BINARY_SHA256SUMS — training identities','provenance/holdout/matrix.json and FREEZE_SHA256SUMS — remaining37 inputs/controller','selection/selection_freeze.json — pre-holdout C3 selection','records/training-campaign, cold-campaign, cohort-campaign — every raw process','records/training-audit, cold-audit, cohort-audit — independent receipts/meters','reports/run-20261015T1516Z — inspected4-page training stage PDF','reports/run-20261015T1842Z — inspected4-page cold stage PDF/native statistics','reports/run-20261015T2311Z — cohort stage report and failed receipt','postprocessing/final_aggregate.cpp — new native binding/profile/statistics export','postprocessing/test_final_aggregate.cpp — coverage/A-A/SHA tamper guards','checks/final — frozen-source seals, native checks, final PDF QA/preservation','reports/final-20261016/native_final.json — source-linked values shown here','delivery — archive manifest, reconstruction paths and CRC/full-hash receipt','problems/n83-eight-target-rho-contract.md — unresolved reference obligation']:
 text(html.escape(s))
text('Study paths relative to research/n83_basis_setup_20261015 except problems/. Reports, editable figures, SVG and150dpiPNG copies accompany benchmark index20261015-n83-basis-setup. Canonical diagnostic changes preserve historical entries and do not promote this result into admitted S/floor/ledger rankings. PR https://github.com/aburan28/crypto/pull/1899 remains dependent on parent1828 until integration. No reviewer is requested without the owner naming one; required human review including Cursor remains a gate.')
text('Python only formats recorded native values, draws figures, renders/inspects PDF pages and packages existing files. No scientific result or experimental logic is recomputed in this presentation. The comprehensive file appends this standalone addendum to the unchanged87-page report (SHAaf881641988a98bec3177e286694233f7b03731e82b6517a51ea78fc5e19793a).')
def footer(canvas,doc):
 canvas.setFont('Helvetica',8);canvas.drawString(48,26,'N83 basis/setup | final evidence | 2026-10-16 | PR1899');canvas.drawRightString(564,26,str(doc.page))
pdf=R/'n83_basis_setup_addendum.pdf';SimpleDocTemplate(str(pdf),pagesize=(612,792),rightMargin=48,leftMargin=48,topMargin=40,bottomMargin=42).build(story,onFirstPage=footer,onLaterPages=footer)
d=fitz.open(pdf);qa={'pages':len(d),'intended_sections':pages,'empty_pages':[],'out_of_bounds':[],'pages_rendered':[]}
for i,p in enumerate(d):
 if len(p.get_text().strip())<50:qa['empty_pages'].append(i+1)
 for b in p.get_text('blocks'):
  if b[0]<0 or b[1]<0 or b[2]>612.1 or b[3]>792.1:qa['out_of_bounds'].append([i+1,list(b[:4])])
 p.get_pixmap(matrix=fitz.Matrix(1,1)).save(R/f'qa-page-{i+1:02}.png');qa['pages_rendered'].append(i+1)
(R/'PDF_QA.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))
# Append old pages without altering the historical input.
old=Path('/workspace/report_output/unified_curve_basis_report_matched_ic_rho.pdf')
merged=fitz.open(old);assert len(merged)==87;merged.insert_pdf(d);merged.save(R/'unified_curve_basis_report_basis_setup.pdf')
for name in ['n83_basis_setup_addendum.pdf','unified_curve_basis_report_basis_setup.pdf']:shutil.copy2(R/name,Path('/workspace/report_output')/name)
