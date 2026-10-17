import argparse,json,time
from unittest.mock import patch
from .core import BatchSeparator
from research.fast_curve_fibers_20261015_v1.panel import cases
from research.fast_curve_fibers_20261015_v1.core import FastFiber,Context
from research.advance import general_chain as gc
from research.advance.test_general_chain import independent_s3

def run(variant):
 rows=[]
 for name,f,curve in cases():
  t=time.monotonic();domain=tuple(range(min(4,f.order)));targets=range(min(8,f.order));answers=[];checked=set();ctx=Context(f,curve)
  def verifier(result,target):
   assert result['search_complete'] and result['output_complete']
   assert result['solution_count_exact']==len(result['solutions']) and len(set(map(tuple,result['solutions'])))==len(result['solutions'])
   for r in result['solutions']:
    assert len(r)==8 and all(x in domain for x in r[:5])
    for v in r:f.validate(v)
    x1,x2,x3,x4,x5,u1,u2,u3=r
    for link in [(x1,x2,u1),(u1,x3,u2),(u2,x4,u3),(u3,x5,target)]:
     key=tuple(sorted(link))
     if key not in checked:
      assert independent_s3(f,*link,curve)==0;checked.add(key)
  def factory(field,coefficients=None,**kwargs):return FastFiber(field,coefficients,context=ctx,**kwargs)
  if variant=='new':
   batch=BatchSeparator(f,domain,curve,max_seconds=30)
   for target in targets:
    r=batch.solve(target);verifier(r,target);answers.append(dict(target=target,count=r['solution_count_exact'],roots=r['solutions']))
  else:
   with patch.object(gc,'GeneralS3Fiber',factory):
    for target in targets:
     r=gc.solve_general_chain(f,domain,target,curve=curve,max_output=100000,max_seconds=30);verifier(r,target);answers.append(dict(target=target,count=r['solution_count_exact'],roots=r['solutions']))
  rows.append(dict(case=name,order=str(f.order),degree=f.degree,seconds=time.monotonic()-t,answers=answers,independently_checked_unique_edges=len(checked)))
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--variant',choices=['old','new'],required=True);a=p.parse_args();t=time.monotonic();print(json.dumps(dict(variant=a.variant,rows=run(a.variant),seconds=time.monotonic()-t),sort_keys=True))
