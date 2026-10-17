from pathlib import Path
import json,gzip,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Image,PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.utils import ImageReader
D=Path(__file__).resolve().parents[2];R=D.parents[1];B=R/'benchmarks/20261016-generalized-rootfinding';B.mkdir(exist_ok=True);run=Path(__file__).resolve().parent;snapshot=D/'historical'
def load(p):
 p=Path(p)
 return json.loads(p.read_text()) if p.exists() else json.loads(gzip.decompress(Path(str(p)+'.gz').read_bytes()))
ids=['batch-separator-manyfield-exact-coldpairs-20261016-v1','support-squarefree-cost-manyfields-20261016-v1','support-squarefree-cost-actual83-20261016-v1']
records=[load(snapshot/'lab/ledger/attempts'/eid/'0001/scientific_result.json') for eid in ids]
data={'source_files':['historical/lab/ledger/attempts/'+eid+'/0001/scientific_result.json' for eid in ids],'batch':{k:v for k,v in records[0]['evidence'].items() if k not in ['records']},'support_many':records[1]['evidence'],'support_actual':records[2]['evidence'],'scope':'Historical Python, descriptive shared-host calibration; no qualified wall or IC/rho claim.'}
(run/'plot-data.json').write_text(json.dumps(data,indent=2)+'\n')
# Rendering recorded data only: no arithmetic/solver experiments are run.
fig,ax=plt.subplots(figsize=(8,10));rows=data['batch']['cases'];ax.barh([x['case'] for x in rows],[x['median_old_over_new'] for x in rows],color='#236d91');ax.invert_yaxis();ax.axvline(1,color='gray',linestyle=':');ax.axvline(5,color='#bc4749',linestyle='--');ax.set_xlabel('Recorded median prior/candidate time ratio (descriptive only)');ax.set_title('Historical complete separator: 20 fields / 38 curve contexts\n3 cold pairs; domain size at most 4; unqualified isolation',fontsize=10);ax.tick_params(axis='y',labelsize=7);fig.tight_layout();fig.savefig(run/'batch.svg');fig.savefig(run/'batch.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,5));rows=data['support_many']['rows'];ax.barh([x['case'] for x in rows],[x['median_old_over_shortcut'] for x in rows],color='#76579b');ax.invert_yaxis();ax.axvline(1,color='gray',linestyle=':');ax.set_xlabel('Recorded power/GCD-shortcut time ratio (descriptive only)');ax.set_title('Historical support verification: 16 fields, 12 supplied quintics each\n3 alternating internal rounds; not cold solver comparisons',fontsize=10);fig.tight_layout();fig.savefig(run/'support-many.svg');fig.savefig(run/'support-many.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4));row=data['support_actual']['rows'][0];ax.plot([1,2,3],row['old_power_seconds'],'o-',label='Quotient-ring power');ax.plot([1,2,3],row['gcd_shortcut_seconds'],'s-',label='Squarefree shortcut');ax.set_xticks([1,2,3]);ax.set_ylabel('Recorded seconds per 16 supplied quintics');ax.set_xlabel('Internal round');ax.set_title('Actual GF(2^83) fixed-parameter support checks\nBoth criteria rejected all 16; not search or root counts',fontsize=10);ax.legend();fig.tight_layout();fig.savefig(run/'support-actual.svg');fig.savefig(run/'support-actual.png',dpi=150);plt.close(fig)
styles=getSampleStyleSheet()
def pdf(path,title,paras,graphs):
 story=[Paragraph(title,styles['Title'])]
 for t in paras:story.extend([Paragraph(t,styles['BodyText']),Spacer(1,9)])
 for name in graphs:
  iw,ih=ImageReader(str(run/(name+'.png'))).getSize();h=460*ih/iw
  if h>540:story.append(PageBreak())
  story.extend([Image(str(run/(name+'.png')),width=460,height=h),Spacer(1,9)])
 SimpleDocTemplate(str(path),pagesize=(595,842),rightMargin=45,leftMargin=45,topMargin=40,bottomMargin=40).build(story)
common=['Historical Python provenance; source and input identities are pinned by TRANSPORT.json and the two frozen SOURCE_MANIFEST.json files. This PR introduces no new solver measurements.','Host: managed Linux x86_64 cloud worker, one CPU and 4 GiB worker budget. Qualified uncontended isolation, A/A control and independent external review were not established. Values are descriptive; no confidence interval was recorded.','All three frozen workers passed. Source snapshots and raw stdout/stderr/process/preregistration/result records are preserved byte-for-byte after gzip decoding. The full older journal is identified by digest but not bundled; this is not a complete resumable lab history.']
for eid,graph in zip(ids,['batch','support-many','support-actual']):
 pdf(run/(graph+'-run-report.pdf'),eid,common+['Command and accounting: see historical/lab/ledger/attempts/'+eid+'/0001/preregistration.json and process.json. Seeds, field parameters, exact command arrays and deadline boundaries remain frozen there.'],[graph])
pdf(D/'RESULTS.pdf','Generalized rootfinding: historical checkpoint',common+['Requested target remains 5x against genuine matched all-roots F4/F5/F6, with the original immutable release gate. That goal is incomplete.','Batch calibration: whole-cold median ratio 1.3830038110595322 against the prior cached separator, 20 fields / 38 curve contexts, three cold pairs. Outputs agree across arms. Tiny domains do not establish hardness.','Support checks: actual GF83 fixed-parameter 16-quintic ratio 1.1408942779538267; both criteria rejected all candidates. Several NIST prime and odd-extension contexts regress under the blanket squarefree shortcut. No original-root count or search completion is claimed.','Historical independent oracle algorithms checked 15 complete small systems and 4,448 support cases. This native archival audit checks payload hashes, pinned closure and recorded paired output equality; it does not independently re-solve equations.','Next work: native solver/harness port before new research comparisons; separately test inversion-free squarefree detection and bounded sparse contraction width for the hard support search. No novel-priority, asymptotic, security or IC/rho claim.'],['batch','support-many','support-actual'])
for svg in run.glob('*.svg'):
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
for p in run.iterdir():
 if p.suffix in ['.pdf','.svg','.png']:shutil.copyfile(p,B/p.name)
shutil.copyfile(D/'RESULTS.pdf',B/'RESULTS.pdf')
(run/'REPORT.md').write_text('# Historical run reports\n\nThree frozen experiments, all passed, with no new timing runs in this PR. See the three run-report PDFs, recorded plot-data.json, editable rendering source, and immutable raw records under ../../historical/lab/ledger/attempts/.\n\nScope: historical Python development calibration; isolation unqualified, no A/A or matched genuine F4/F5/F6 panel. The graph values are descriptive only. Full-system counts remain unset for support experiments.\n')
print('reports and benchmark copies rendered')
