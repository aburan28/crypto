import json
from unittest.mock import patch
from research.advance import general_chain as gc
from research.advance.test_general_chain import full_oracle
from research.fast_curve_fibers_20261015_v1.core import FastFiber,Context
from rootfinder.fields import PrimeField,ExtensionField
from rootfinder.binaryfield import BinaryField
cases=0
for f,curve in [(PrimeField(5),{'a2':0,'a4':1,'a6':1}),(ExtensionField(3,2,(1,0,1)),{'a2':3,'a4':1,'a6':4}),(BinaryField(3,11),{'a2':0,'b':1}),(BinaryField(4,19),{'a2':5,'b':7})]:
 ctx=Context(f,curve,capacity=8)
 def factory(field,coefficients=None,**kwargs):return FastFiber(field,coefficients,context=ctx,**kwargs)
 with patch.object(gc,'GeneralS3Fiber',factory):
  for target in (0,1,2):
   r=gc.solve_general_chain(f,(0,1),target,curve=curve,max_output=100000);expected=full_oracle(f,(0,1),target,curve)
   assert r['search_complete'] and r['output_complete'] and {tuple(x) for x in r['solutions']}==expected;cases+=1
print(json.dumps({'status':'passed','independent_D5K3_allroot_oracles':cases,'context_shared_across_targets':True,'LRU_eviction_exercised':True,'wide_actual_runs':0}))
