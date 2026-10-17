# Presentation only: uses recorded native evidence; no research computation.
import json,datetime,textwrap
from pathlib import Path
import fitz
from reportlab.pdfgen import canvas
r=Path(__file__).resolve().parents[2];d=Path(__file__).parent
eid='native-long-XOR-interface-complete-controls-20261017-v1'
a=r/'lab/ledger/attempts'/eid/'0001'
result=json.loads((a/'scientific_result.json').read_text());assert result['status']=='passed',result
text='''Native long-XOR interface validation

Hypothesis: a fresh source-pinned wrapper can accept long XOR rows while retaining exact functional-primary models despite native internal cut variables. The unsupported implied_by entry point is disabled. No native core or Gaussian settings are changed.

The frozen worker passed. It enumerated all128 expected even-parity8-bit models, checked against all256 primary vectors by independent integer-popcount oracle. All256 declared AND outputs and the parity output were replayed. The265-variable XOR row was accepted; an empty true-parity row was rejected. Native statistics recorded265 declared and396 internal coordinates, with0 auxiliary decisions. This exercises internal coordinates on this fixture, not arbitrary long-row correctness.

Native integration test1/1 passed with its pinned library. It is ignored by default in ordinary CI. All26 native configuration readback values were exact. Wrapper/core hashes,159 core-source hashes, input contract and30-second/4GiB/CPU0 budget were frozen. The streaming compatibility coordinator verified the journal; the outer supervisor charged child startup and complete lifetime.

This is interface correctness evidence, not Gaussian admission, actual compressed-norm validation, a hard-root solve, novelty or a speed comparison. The single-row fixture initialized no Gaussian matrix. Original all-roots coverage and the matching genuine F4/F5/F6 end-to-end gate remain open. No release claim is made.

Source review found that add_every_combination_xor retains the original XOR but also introduces CNF cut coordinates. The265-variable fixture created131 such coordinates. A dense full-cone interface may therefore increase native work despite fewer external variables. Exact actual counts after cleaning remain unmeasured.

Next: first profile a cut-aware hybrid cone representation; then encode the canonical affine definitions and AND clauses of the exact norm compiler; compare complete tiny primary model sets, including zero/repeated sources, against the independent literal norm oracle. Only then freeze a new hard search. Do not retry the unchanged old timeout or infer speed from representation counts.

Publication review and explicit authorization extended this interaction beyond the nominal15-minute cycle. No scientific worker budget was enlarged. One new bounded interface experiment ran. All earlier attempts are preserved.

Evidence: research/native_long_xor_interface_20261017_v1/ and lab/ledger/attempts/native-long-XOR-interface-complete-controls-20261017-v1/0001/. Canonical performance graphs remain unchanged because no new qualified performance result exists.
'''
(d/'REPORT.md').write_text(text)
report={'schema':'rootfinding_cycle_report_v1','experiment_id':eid,'result':result,'original_goal_complete':False,'release_gate_claim':False,'next_experiment':'native-affine-cone-norm-integration-v2','timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'interface validation only','core_unchanged':True}
(d/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n')
svg='<svg xmlns="http://www.w3.org/2000/svg" width="760" height="210"><rect width="760" height="210" fill="white"/><g font-family="sans-serif" fill="#153c55" font-size="16"><text x="20" y="30">8 primary bits; 265-variable parity row</text><text x="20" y="70">256 declared AND outputs, all replayed</text><text x="20" y="110">396 internal coordinates; 0 auxiliary decisions</text><text x="20" y="150">128 complete models match independent exhaustive oracle</text><text x="20" y="190">Interface validation only: no speed measurement</text></g></svg>'
(d/'interface.svg').write_text(svg)
sv=fitz.open(stream=svg.encode(),filetype='svg');sv[0].get_pixmap(matrix=fitz.Matrix(150/72,150/72)).save(d/'interface.png')
c=canvas.Canvas(str(d/'RESULTS.pdf'),pagesize=(612,792));y=760
for para in text.split('\n'):
 for line in textwrap.wrap(para,100)or['']:
  if y<45:c.showPage();y=755
  c.setFont('Helvetica',8);c.drawString(32,y,line);y-=11
if y<210:c.showPage();y=755
c.drawImage(str(d/'interface.png'),30,y-155,width=550,height=152,preserveAspectRatio=True)
c.save();pdf=fitz.open(d/'RESULTS.pdf')
for page in pdf:page.get_pixmap(matrix=fitz.Matrix(.8,.8))
print(json.dumps({'report':str(d),'pages':len(pdf),'worker_status':result['status']}))
