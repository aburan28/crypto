import argparse,json,time,statistics,hashlib,random
from pathlib import Path
from .core import product_h,prefix_parameters,support,squarefree
from research.fast_curve_fibers_20261015_v1.panel import cases
from rootfinder.binaryfield import BinaryField
from research.advance.five_divisor import prepare_five_model,conditional_five_family,five_coefficients_formula
from rootfinder import polynomial as p

def run(mode,started):
 rows=[]
 if mode=='many':
  selected=set();jobs=[]
  for name,f,curve in cases():
   key=(f.order,f.degree)
   if key in selected or f.order<8:continue
   selected.add(key)
   if f.characteristic==2 and f.degree>=18:D=list(range(8));domain={'kind':'prefix','betas':prefix_parameters(f,[1<<i for i in range(18)])};outside=1<<18
   else:D=list(range(min(8,f.order-1)));domain={'kind':'explicit','polynomial':product_h(f,D)};outside=len(D)
   patterns=[[0,1,2,3,4],[0]*5,[1]*4+[2],[0,0,1,1,2],[2,2,3,4,4],[3]*5];Hs=[product_h(f,rs) for rs in patterns]+[product_h(f,rs[:-1]+[outside]) for rs in patterns];expected=[True]*6+[False]*6;jobs.append((name,f,domain,Hs,expected))
 else:
  path=Path('research/advance/instances/n83-d18-m5-T001-offset0.json');blob=path.read_bytes();assert hashlib.sha256(blob).hexdigest()=='5c9be468c207cc2b013123b4265ee0963e96490c16aad1b1650a3528c0ce0c25';d=json.loads(blob);f=BinaryField(83,0x800000000200000000007);basis=tuple(int(x,16) for x in d['source_basis_hex']);assert basis==tuple(1<<i for i in range(18));model=prepare_five_model(f,int(d['curve']['a_hex'],16),int(d['curve']['b_hex'],16),(int(d['target']['x_hex'],16),int(d['target']['y_hex'],16)),0,basis);h=0x4208336d0d441df79c400;k=0;fam=conditional_five_family(model,h,k);assert fam['nullity']==73
  rng=random.Random(20261016);Hs=[]
  for i in range(16):
   v=fam['particular']
   for b in fam['basis']:
    if rng.randrange(2):v^=b
   Hs.append(five_coefficients_formula(model,h,k,v&((1<<83)-1),v>>83))
  jobs=[('actual83-fixedh-k-relaxed16',f,{'kind':'prefix','betas':prefix_parameters(f,basis)},Hs,None)]
 for name,f,domain,Hs,expected in jobs:
  seconds={'power':[],'gcdshortcut':[]};answers={}
  for rep in range(3):
   for method in (('power','gcdshortcut') if rep%2==0 else ('gcdshortcut','power')):
    t=time.monotonic();values=[support(f,H,domain,method=='gcdshortcut') for H in Hs];seconds[method].append(time.monotonic()-t)
    if method in answers:assert answers[method]==values
    answers[method]=values
  assert answers['power']==answers['gcdshortcut']
  if expected is not None:assert answers['power']==expected
  rows.append(dict(case=name,order=str(f.order),degree=f.degree,candidates=len(Hs),squarefree_candidates=sum(squarefree(f,H) for H in Hs),accepted=answers['power'],old_power_seconds=seconds['power'],gcd_shortcut_seconds=seconds['gcdshortcut'],median_old_over_shortcut=statistics.median(seconds['power'])/statistics.median(seconds['gcdshortcut']),original_root_count_exact=None,auxiliary_root_count_exact=None))
  if time.monotonic()-started>=70:raise TimeoutError('charged support panel')
 return dict(rows=rows,mode=mode,release_gate_claim=False,scope='exact support verification of suppliedquintics, not search or enumeration',original_root_count_exact=None,auxiliary_root_count_exact=None)
def main():
 t=time.monotonic();q=argparse.ArgumentParser();q.add_argument('--experiment-id',required=True);q.add_argument('--mode',choices=['many','actual'],required=True);a=q.parse_args()
 try:e=run(a.mode,t);status='passed';summary='Exact nilpotence and squarefree-shortcut support criteria compared with all repeated roots retained.'
 except Exception as x:e={'release_gate_claim':False,'original_root_count_exact':None};status='inconclusive';summary=type(x).__name__+': '+str(x)
 r=dict(schema='scientific_result_v1',experiment_id=a.experiment_id,status=status,summary=summary,evidence=e,next_actions=['Use measured inversion/GCDoverhead to selectshortcut; ifslower refuseblanketpromotion. Compactverificationstillneedsconditional search acrossrelaxation.']);s=json.dumps(r,sort_keys=True)
 if time.monotonic()-t>=75:r.update(status='inconclusive',summary='Latesupportserialization.',evidence={'release_gate_claim':False,'original_root_count_exact':None});s=json.dumps(r,sort_keys=True)
 print(s)
if __name__=='__main__':main()
