from itertools import product,combinations_with_replacement
from rootfinder.binaryfield import BinaryField
from rootfinder.fields import PrimeField
from .core import product_h,support,prefix_parameters
import json
cases=0
for f,D in [(BinaryField(2,7),(0,1)),(PrimeField(3),(0,1))]:
 valid={product_h(f,rs) for rs in combinations_with_replacement(D,5)};domain={'kind':'explicit','polynomial':product_h(f,D)}
 for cs in product(range(f.order),repeat=5):
  H=cs+(1,);expected=H in valid
  assert support(f,H,domain,False)==expected and support(f,H,domain,True)==expected;cases+=1
f=BinaryField(4,19);betas=prefix_parameters(f,[1,2]);domain={'kind':'prefix','betas':betas}
for rs in combinations_with_replacement(range(4),5):assert support(f,product_h(f,rs),domain);cases+=1
print(json.dumps({'status':'passed','independent_enumerated_support_cases':cases,'repeatedroots_inseparable_char3_and2':True,'actualwide_runs':0,'external_independent_review':False}))
