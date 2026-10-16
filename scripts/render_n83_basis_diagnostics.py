#!/usr/bin/env python3
"""Presentation only: copy frozen native N83 values into canonical diagnostics.
No new scientific arithmetic, estimators, bootstrap, or ledger admission.
"""
import json,re,html
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STUDY='research/n83_basis_setup_20261015'
C=json.loads((ROOT/STUDY/'reports/run-20261015T1842Z/native_statistics.json').read_text())
F=json.loads((ROOT/STUDY/'reports/final-20261016/native_final.json').read_text())
assert F['cohort_IC_over_rho'] is None and F['canonical_answers']==92
source='https://github.com/aburan28/crypto/blob/main/'+STUDY+'/reports/final-20261016/native_final.json'
coldsource='https://github.com/aburan28/crypto/blob/main/'+STUDY+'/reports/run-20261015T1842Z/native_statistics.json'
p=ROOT/'docs/ic/progress-timeline.json';d=json.loads(p.read_text());ids=['n83-basis-cold-c3','n83-basis-cold-c0','n83-basis-online-c3']
d['series']=[v for v in d['series'] if v['id'] not in ids]
entries=[(ids[0],'Training-selected C3, eight actual cold singleton blocks',C['comparisons']['selected_over_rho']['mean_cost_ratio'],C['comparisons']['selected_over_rho']['conditional_95_percent_interval'],coldsource,'actual complete-process wall including preparation'),(ids[1],'Retained K600 control, eight fresh cold singleton blocks',C['comparisons']['control_over_rho']['mean_cost_ratio'],C['comparisons']['control_over_rho']['conditional_95_percent_interval'],coldsource,'actual complete-process wall including preparation'),(ids[2],'C3 aligned online query mean, secondary diagnostic',F['online_cold_selected_over_rho'],None,source,'aligned validation/solve/replay/output clocks; secondary')]
for identity,label,value,interval,src,unit in entries:
 point={'date':'2026-10-16','measurement_date':'2026-10-15','publication_status':'PR1899 unmerged; date provisional and must refresh on actual merge.','value':value,'reference':'matched','unit':unit,'label':label,'detail':'Native IC/rho; below1 favors IC. CPU3/oneOSthread/12GiB/2100s,one guest. C3 selected before holdouts; old K600 construction and OS-cache flushing unpriced. Three cohort rho input-contract failures are not replaced. No calibrated S/L2/independent-host or ledger admission; provisional publication date pending PR1899 merge.','class':'engineering','anchor':'#n83-basis-setup-20261015','source':src}
 if interval is not None:point['interval_95']=interval
 d['series'].append({'id':identity,'label':label,'short_label':identity,'role':'context','barrier_test':False,'note':'Native IC/rho values in local measured units, outside admitted best/floor tables. Historical prior-study values remain separate.','points':[point]})
d['updated']='2026-10-16';p.write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
p=ROOT/'docs/index-calculus-scoreboard.html';s=p.read_text();s=re.sub(r'<script type="application/json" id="progress-data">.*?</script>','<script type="application/json" id="progress-data">'+json.dumps(d,indent=2,ensure_ascii=False)+'</script>',s,flags=re.S)
s=re.sub(r'<!-- BEGIN n83-basis-setup-20261015 -->.*?<!-- END n83-basis-setup-20261015 -->\n?','',s,flags=re.S)
rows=''
for k,label in [('selected_over_rho','C3/rho'),('control_over_rho','K600/rho'),('selected_over_control','C3/K600')]:
 v=C['comparisons'][k];rows+=f'<tr><td>{label}</td><td>{v["mean_cost_ratio"]:.9f}</td><td>[{v["conditional_95_percent_interval"][0]:.9f}, {v["conditional_95_percent_interval"][1]:.9f}]</td></tr>'
panel='<!-- BEGIN n83-basis-setup-20261015 -->\n<section class="panel" id="n83-basis-setup-20261015"><h2>Degree-83 basis/setup: C3 halves K600 cost; rho cheaper cold <span class="chip">engineering diagnostic</span></h2><p class="panel-summary">Exact legacy a=1/tail149,53-bit subgroup. Eight fresh matched cold blocks: C3 mean611.139215s, K6001223.782149s, optimized rho2.107724s. Below1 favors IC for IC/rho. C3/K600 explicitly compares bases.</p><div class="table-scroll"><table><thead><tr><th>Cold boundary</th><th>Native mean cost ratio</th><th>Conditional95% interval</th></tr></thead><tbody>'+rows+'</tbody></table></div><p>Training-selected C3 prefix1271dense, selection committed before holdouts. Actual55 attempts/50complete/5failed,92canonical answers,30full rank/log audits,5identical A/A workpairs. Two C4 allocation failures retained. Three8-target rho cohorts rejected input because frozen worker accepts1 or32points: cohort IC/rho unavailable, no replacement times or retry. Six IC cohorts independently verified; C3/K6000.513799569 [0.505576379,0.522482173], only3prep clusters.</p><p>One guest, CPU3,one thread,12GiB,2100s; physical host not exclusive. Online aligned query ratio C3/rho'+format(F['online_cold_selected_over_rho'],'.9f')+' is secondary. Original retained K600 construction/OS-cache flushing unpriced, memory-limit flagnull. S/floor/admitted-best rankings unchanged; no L2/independent-host/global optimum/asymptotic/security admission. Prior studies stay separate, with no timing mixing or multiplication of gains. PR1899 unmerged; publication date provisional. <a href="'+source+'">Native full accounting/profiles</a> · <a href="'+coldsource+'">Native cold statistics</a> · <a href="https://github.com/aburan28/crypto/blob/main/'+STUDY+'/REPORT.md">Final report and unresolved reference gap</a>.</p></section>\n<!-- END n83-basis-setup-20261015 -->\n'
s=s.replace('<!-- BEGIN n83-basis-setup-20261015 -->','') if False else s
# Place a top-level diagnostic immediately before the progress chart panel.
anchor=s.index('<details class="dash-details"><summary>Every plotted point as a table</summary>')
# Existing enclosing progress section starts earlier; insert before its section element.
progress=s.rfind('<section',0,anchor)
s=s[:progress]+panel+s[progress:]
start=s.index('<summary>Every plotted point as a table</summary>');end=s.index('</tbody>',start)
newrows='<!-- BEGIN n83-basis-table -->\n'
for identity,label,value,interval,src,unit in entries:
 extra='' if interval is None else f' <small>95%: {interval[0]:.9f}&ndash;{interval[1]:.9f}</small>'
 newrows+=f'<tr><td>2026-10-16</td><td>{html.escape(label)}</td><td>Measured2026-10-15, PR1899 unmerged; landing date provisional</td><td class="n">{value:.9f}&times;{extra}</td><td>matched, local diagnostic</td><td><span class="chip">engineering</span></td><td><a href="#n83-basis-setup-20261015">panel</a> · <a href="{src}">native source</a></td></tr>\n'
newrows+='<!-- END n83-basis-table -->\n'
s=re.sub(r'<!-- BEGIN n83-basis-table -->.*?<!-- END n83-basis-table -->\n?', '', s, flags=re.S)
s=re.sub(r'<tr data-n83-basis=".*?</tr>\n?','',s,flags=re.S);start=s.index('<summary>Every plotted point as a table</summary>');end=s.index('</tbody>',start);s=s[:end]+newrows+s[end:];p.write_text(s)
print('Copied native N83 cold/online diagnostics; missing cohort rho kept unavailable')
