import json
from .fields import BudgetedBinaryField
from .adapter import curve_capability,solve_chain,solve_polynomials
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField
n=0
for degree,mod in [(1,3),(2,7),(3,11),(4,19),(8,283)]:
 a=BudgetedBinaryField(degree,mod);b=BinaryField(degree,mod)
 for x in range(min(a.order,16)):
  assert a.square(x)==b.square(x)
  for y in range(min(a.order,16)):assert a.mul(x,y)==b.mul(x,y);n+=1
for degree,mod in [(3,9),(4,21)]:
 try:BudgetedBinaryField(degree,mod)
 except ValueError:pass
 else:raise AssertionError('reducible accepted')
f=BinaryField(3,11)
try:curve_capability(f,'koblitz',{'a2':2,'b':1})
except ValueError:pass
else:raise AssertionError('invalid Koblitz')
r=solve_chain(PrimeField(5),(0,1),2,family='prime',curve={'a2':0,'a4':1,'a6':1},max_output=0)
assert not r['complete']
r=solve_polynomials(PrimeField(5),[[(1,(2,)),(4,(0,))]],1,(0,1,2,3,4));assert r['complete'] and {tuple(x) for x in r['roots']}=={(1,),(4,)}
print(json.dumps({'status':'passed','wide_backend_matches_existing_small_multiplications':n,'invalid_moduli_and_Koblitz_refused':True,'output_cap_does_not_claim_allroots':True}))
