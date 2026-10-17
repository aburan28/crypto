"""Render recorded evidence only; never runs an algorithm or timing harness."""
import json,pathlib,statistics,math,datetime,hashlib,zipfile
D=pathlib.Path(__file__).parent;root=D.parents[1]
ids=['native-quintic-rational-repeated-support-20261017-v1','native-fixed-u-affine-cell-pruning-20261017-v1','native-quintic-quadratic-exceptional-support-20261017-v2']
records=[]
for eid in ids:
 p=root/'lab/ledger/attempts'/eid/'0001';r=json.loads((p/'scientific_result.json').read_text());proc=json.loads((p/'process.json').read_text());record={'id':eid,'status':r['status'],'result_path':str(p/'scientific_result.json'),'cold_seconds_diagnostic':proc['wall_seconds'],'evidence':r['evidence']}
 if 'seconds_internal_diagnostic' in r['evidence']['rows'][0]:
  rows=[]
  for row in r['evidence']['rows']:
   a,b=map(statistics.median,row['seconds_internal_diagnostic']);aa=row['aa_internal_seconds'];rows.append({'field':row['field'],'internal_gain':a/b,'aa_spread':(max(aa)-min(aa))/min(aa)})
  record['analysis']={'rows':rows,'geometric_gain':math.exp(sum(math.log(x['internal_gain'])for x in rows)/len(rows)),'noise_limit_exceeded':any(x['aa_spread']>.1 for x in rows),'qualified_speed_claim':False}
 records.append(record)
now=datetime.datetime.now(datetime.timezone.utc).isoformat();nextq=[{'id':'actual-affine-family-U-rank-native-v1','budget_seconds':90,'memory_mib':4096,'cores':1,'acceptance':'Independently reproduce actual rank93/nullity73; compute U ranks on the four fixed h,k cells. If U rank<=8 admit <=256 conditional cells; otherwise refuse exponential partitioning and profile variable-U bilinear incidence.'},{'id':'conditional-support-v3-cold-cost-native','budget_seconds':90,'memory_mib':4096,'cores':1,'acceptance':'New unopened source-bound hard-family candidates and independent arithmetic reference; include field/domain/input setup, cold output, all per-field regressions; expand A/A to distinguish any noisy fields. No inherited stage gain or original all-roots qualification.'},{'id':'odd-prime-domain-native-adapter','budget_seconds':90,'memory_mib':4096,'cores':1,'acceptance':'Exact interval/lift root extraction on odd prime fields, repeated/nonrational and wrapped-coefficient controls; extension subspace recurrence alone does not cover prime intervals.'}]
report={'date':now,'records':records,'next_queue':nextq,'source_pr':'https://github.com/aburan28/crypto/pull/2049','goal_complete':False,'requested_gain':8,'qualified_gain':None,'original_system_complete':False,'novelty_established':False,'publication':'Reports and raw evidence local pending explicit approval; source-only PR updated.','accounting':'Internal support timings exclude domain/field/input construction and are amortized batches; cold worker time includes all stages but has no paired cold reference. No qualified isolation or external-fastest comparison.','cell_transfer':'Six synthetic fixed-U affine cells leave dimension8. No claim that the actual norm family is reduced to2^8. See TRANSFER_REVIEW.md.'}
(D/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['Rootfinding bounded cycle: three native experiments','',f'Date: {now}','','All registered workers finished in their declared component scope; check individual statuses below. No original-system completion, novelty or qualified8x result.','']
for r in records:
 lines += [r['id']+': '+r['status']]
 if 'analysis'in r:lines += [f"Internal mixed-predicate geometric gain: {r['analysis']['geometric_gain']:.4f}; A/A above10%: {r['analysis']['noise_limit_exceeded']}. Unqualified descriptive comparison."]
 else:lines += ['Fixed-U cell dimensions: '+', '.join(x['field']+': '+str(x['parameter_dimension'])+' -> '+str(x['surviving_dimension'])for x in r['evidence']['rows'])+'. Conditional synthetic affine space only.']
lines += ['',report['accounting'],'',report['cell_transfer'],'','Next: actual-family U ranks, cold complete component costs and the separate odd prime interval adapter.','Reports/raw evidence remain local because automatic approval review rejected publication without explicit report-payload authorization. Source code is committed in PR2049.']
(D/'REPORT.md').write_text('\n'.join(lines)+'\n')
# Presentation only: plot recorded aggregates, with limitations on the figure.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(8,4.5));rows1=records[0]['analysis']['rows'];rows2=records[2]['analysis']['rows'];x=list(range(6));ax.bar([i-.18 for i in x],[r['internal_gain']for r in rows1],width=.35,label='v1');ax.bar([i+.18 for i in x],[r['internal_gain']for r in rows2],width=.35,label='v2');ax.set_xticks(x,[r['field']for r in rows1]);ax.axhline(1,color='gray',linestyle=':');ax.set_ylabel('Internal reference / candidate median');ax.set_title('Supplied-quintic support predicate: diagnostic only\n5 internal rounds, domain setup excluded; no external/Groebner or8x claim',fontsize=10);ax.legend();fig.tight_layout();fig.savefig(D/'support.svg');fig.savefig(D/'support.png',dpi=130);plt.close(fig)
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image
from reportlab.lib.styles import getSampleStyleSheet
styles=getSampleStyleSheet();story=[Paragraph('Exact rootfinding component cycle',styles['Title'])]
for line in lines[2:]:
 if line:story += [Paragraph(line.replace('&','&amp;').replace('<','&lt;'),styles['BodyText']),Spacer(1,7)]
story.append(Image(str(D/'support.png'),width=490,height=276));SimpleDocTemplate(str(D/'RESULTS.pdf')).build(story)
with (root/'lab/decisions.jsonl').open('a')as f:f.write(json.dumps({'date':now,'experiments':ids,'decision':'Keep exact branches; v1 timing is unqualified due GF409 A/A. Fixed-U whole-cell reduction is conditional, not actual norm-family pruning. v2 is separate evidence; no multiplied stage gains.','report':str(D/'REPORT.json'),'next_queue':nextq,'goal_complete':False})+'\n')
hp=root/'lab/hypotheses.json';h=json.loads(hp.read_text());key='hypotheses'if 'hypotheses'in h else next(k for k,v in h.items()if isinstance(v,list)and v and isinstance(v[0],dict));h[key].extend({'id':r['id'],'status':r['status'],'report':str(D/'REPORT.json'),'scope':'component only; not original-family completion'}for r in records);h[key].extend(dict(q,status='next_bounded_design',priority=1)for q in nextq);h['updated_utc']=now;hp.write_text(json.dumps(h,indent=2)+'\n')
cp=root/'lab/CYCLE_REPORT.json';c=json.loads(cp.read_text());c['registered_experiment_count']=122;c['updated_utc']=now;c['latest_bounded_cycle_report']=str(D/'REPORT.json');c['next_priority']=nextq;c['research_program_state']='active_research';c['actual_attempts'].extend({'experiment_id':r['id'],'status':r['status'],'result_path':r['result_path'],'summary':'Native component experiment; no original-system completion or qualified gain.'}for r in records);c['journal_records']=json.loads((root/'research/quintic_rational_support_20261017_v2/RUN_READY.json').read_text())['journal_records'];cp.write_text(json.dumps(c,indent=2)+'\n')
files=[]
for directory in [D,root/'research/quintic_rational_support_20261017_v1',root/'research/quintic_fixed_u_cells_20261017_v1',root/'research/quintic_rational_support_20261017_v2']:
 files += [p for p in directory.rglob('*')if p.is_file()and '__pycache__'not in p.parts]
files += [root/'research/quintic_resultant_20261016_v1/src/main.rs']
for eid in ids:files += [p for p in (root/'lab/ledger/attempts'/eid/'0001').iterdir()if p.is_file()]
(D/'ARTIFACT_SHA256.json').write_text(json.dumps({str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files},indent=2)+'\n');files.append(D/'ARTIFACT_SHA256.json')
z=pathlib.Path('/workspace/rootfinding_conditional_cycle_20261017_v2.zip')
with zipfile.ZipFile(z,'x',zipfile.ZIP_DEFLATED)as bundle:
 for p in files:bundle.write(p,str(p.relative_to(root)))
print(json.dumps({'report':str(D/'REPORT.json'),'pdf':str(D/'RESULTS.pdf'),'archive':str(z),'sha256':hashlib.sha256(z.read_bytes()).hexdigest(),'comparisons':[r.get('analysis')for r in records]}))
