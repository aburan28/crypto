import argparse,json,time
from unittest.mock import patch
from research.advance import general_chain as gc
from research.fast_curve_fibers_20261015_v1.core import FastFiber,Context
from research.curve_portfolio_20261015_v1.fields import BudgetedBinaryField
from lab.research.nist_fields import CertifiedPrimeField,preset

def cases():
 for name in ['P192','P224','P521']:
  r=preset(name);yield name,CertifiedPrimeField(r['p']),r['curve']
 for n,m,a,b in [(83,0x800000000200000000007,0,1),(163,(1<<163)|(1<<7)|(1<<6)|(1<<3)|1,1,1),(283,(1<<283)|(1<<12)|(1<<7)|(1<<5)|1,3,7),(571,(1<<571)|(1<<10)|(1<<5)|(1<<2)|1,0,1)]:yield 'binary-'+str(n),BudgetedBinaryField(n,m),{'a2':a,'b':b}

def job(variant):
 rows=[]
 for name,f,curve in cases():
  t=time.monotonic();ctx=Context(f,curve);answers=[]
  def factory(field,coefficients=None,**kwargs):return FastFiber(field,coefficients,context=ctx,**kwargs)
  def batch():
   for target in range(8):
    r=gc.solve_general_chain(f,(0,1,2,3),target,curve=curve,max_output=100000,max_seconds=30)
    if not r['search_complete'] or not r['output_complete']:raise RuntimeError('incomplete chain '+name)
    answers.append(dict(target=target,root_count=r['solution_count_exact'],roots=r['solutions']))
  if variant=='new':
   with patch.object(gc,'GeneralS3Fiber',factory):batch()
  else:batch()
  rows.append(dict(case=name,seconds=time.monotonic()-t,answers=answers,cache_hits=ctx.hits,cache_misses=ctx.misses))
 return rows
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--variant',choices=['old','new'],required=True);a=p.parse_args();t=time.monotonic();print(json.dumps(dict(variant=a.variant,rows=job(a.variant),seconds=time.monotonic()-t),sort_keys=True))
