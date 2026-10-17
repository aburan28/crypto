# Presentation of frozen native outputs only; no solver or research oracle.
import json, datetime, textwrap, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/rootfinding_cycle9_plotcache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitz
from reportlab.pdfgen import canvas

r=Path(__file__).resolve().parents[2]
d=Path(__file__).parent
v=r/'research/hybrid_norm_xor_20261017_v1'
eid='native-cut-aware-hybrid-norm-profile-n83-20261017-v1'
a=r/'lab/ledger/attempts'/eid/'0001'
result=json.loads((a/'scientific_result.json').read_text())
process=json.loads((a/'process.json').read_text())
e=result.get('evidence',{})
arms=e.get('arms',[])
control=json.loads((v/'CONTROLS.json').read_text())
report={'schema':'rootfinding_cycle_report_v1','experiment_id':eid,'native_result':result,'controls':control,'process_record':str((a/'process.json').relative_to(r)),'source_freeze':str((v/'FREEZE.json').relative_to(r)),'experiment_count':1,'registered_experiment_count':134,'goal_complete':False,'release_gate_claim':False,'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'next_experiment':'native-hybrid-initialization-and-new-hard-search-v2' if e.get('prediction_at_least10percent_internal_reduction') else 'native-shared-parity-basis-optimization-v2'}
(d/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['Native cut-aware norm representation profile','',
 'Hypothesis: at least one bounded affine support limit reduces actual native internal coordinates by10% versus the exact limit1 reference. This is a representation-count prediction, not a speedup claim.',
 '', 'The same174-primary-bit actual GF(2^83),d18 norm predicate was used under every arm. The source, fixture hash, seven limits and150000 preclean-coordinate admission estimate cap were frozen before this target profile. FourGiB, CPU0 and90-second outer budget were retained. No SAT search or Gaussian initialization ran.',
 '', 'Independent controls passed73216 literal polynomial-norm assignments,116480 hybrid complete-primary assignment checks, and seven native model sets of28 models each. Unique functional extension gives a bijection on primary models. Constants, zero and repeated sources are retained.',
 '', 'Recorded arms:']
for x in arms:
 lines.append('limit '+str(x['support_limit'])+': external '+str(x['external_variables'])+', XOR rows '+str(x['xor_rows'])+', actual native internal '+str(x.get('native_internal_variables_after_flush','not loaded: preregistered cap'))+'.')
lines += ['', 'Worker status: '+result['status']+'. Profile complete: '+str(e.get('profile_complete'))+'. Preregistered10% prediction met: '+str(e.get('prediction_at_least10percent_internal_reduction'))+'.',
 '', 'The native count includes cut coordinates created after explicit pending-addition flush. Raw coefficient counts and preclean estimates remain in the result. A skipped dense arm is capped, not an empty root set or an UNSAT result. No original-root count, timing gain, new field-family performance coverage or novelty is claimed.',
 '', 'Decision and next experiment: '+report['next_experiment']+'. Retain every arm, use full native matrix statistics and literal witness replay before hard search, and keep original S3 ghost/exception/auxiliary coverage and the matching genuine F4/F5/F6 eightfold end-to-end goal unresolved.',
 '', 'Novelty review: the author CryptoMiniSat documentation and the NIST Boyar/Matthews/Peralta2008 record establish earlier XOR/circuit-minimization work. PRIMARY_SOURCES.md limits this review to opened documentation and abstract; no new priority is claimed.',
 '', 'Evidence: research/hybrid_norm_xor_20261017_v1/ and lab/ledger/attempts/'+eid+'/0001/. Canonical performance graphs remain unchanged because this is not a qualified end-to-end result. New report publication remains pending the earlier separate authorization; source code was committed to PR2049.']
text='\n\n'.join(lines)
(d/'REPORT.md').write_text(text+'\n')
fig,axes=plt.subplots(1,2,figsize=(11,3.7))
labels=[str(x['support_limit']) for x in arms]
axes[0].bar(range(len(arms)),[x['external_variables'] for x in arms],color='#237a94')
idx=[i for i,x in enumerate(arms) if 'native_internal_variables_after_flush' in x]
axes[1].bar(idx,[arms[i]['native_internal_variables_after_flush'] for i in idx],color='#c7752e')
for i,x in enumerate(arms):
 if i not in idx:axes[1].text(i,0,'capped',rotation=90,ha='center',va='bottom')
for ax,title in zip(axes,['External variables','Actual native internal variables after flush']):
 ax.set_xticks(range(len(arms)),labels);ax.set_xlabel('Affine support limit (1 = reference)');ax.set_ylabel('Variables (count)');ax.set_title(title);ax.grid(axis='y',alpha=.2)
fig.suptitle('Actual GF(2^83) norm circuit; one fixed instance, count profile only')
fig.tight_layout();fig.savefig(d/'variables.svg');fig.savefig(d/'variables.png',dpi=150);plt.close(fig)
c=canvas.Canvas(str(d/'RESULTS.pdf'),pagesize=(612,792));y=760
for para in text.split('\n'):
 for line in textwrap.wrap(para,100)or['']:
  if y<45:c.showPage();y=755
  c.setFont('Helvetica',8);c.drawString(32,y,line);y-=11
if y<240:c.showPage();y=755
c.drawImage(str(d/'variables.png'),25,y-200,width=562,height=189,preserveAspectRatio=True)
c.save()
pdf=fitz.open(d/'RESULTS.pdf')
for page in pdf:page.get_pixmap(matrix=fitz.Matrix(.7,.7))
print(json.dumps({'status':result['status'],'report':str(d),'pdf_pages':len(pdf)}))
