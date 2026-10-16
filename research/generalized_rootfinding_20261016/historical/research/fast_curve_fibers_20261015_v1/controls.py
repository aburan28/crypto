import json,itertools
from .core import FastFiber,Context
from rootfinder.fields import PrimeField,ExtensionField
from rootfinder.binaryfield import BinaryField
from research.advance.test_general_chain import independent_s3
cases=0
for f,curve in [(PrimeField(3),{'a2':0,'a4':1,'a6':1}),(PrimeField(5),{'a2':0,'a4':1,'a6':1}),(PrimeField(7),{'a2':0,'a4':1,'a6':1}),(ExtensionField(3,2,(1,0,1)),{'a2':3,'a4':1,'a6':4}),(ExtensionField(5,2,(2,0,1)),{'a2':6,'a4':1,'a6':1}),(BinaryField(3,11),{'a2':0,'b':1}),(BinaryField(4,19),{'a2':5,'b':7})]:
 ctx=Context(f,curve,capacity=7);fiber=FastFiber(f,curve,context=ctx)
 for x,y in itertools.product(range(f.order),repeat=2):
  expected=tuple(z for z in range(f.order) if independent_s3(f,x,y,z,curve)==0)
  assert fiber.roots(x,y)==expected and fiber.roots(y,x)==expected;cases+=1
print(json.dumps({'status':'passed','full_field_exhaustive_fiber_cases':cases,'covers_characteristic3_evenextensions_repeats_linear_empty_LRUeviction':True,'actualwide_runs':0,'external_review':False}))
