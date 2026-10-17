from itertools import product,combinations_with_replacement
from rootfinder.fields import PrimeField
from research.support_cost_20261016_v1.core import product_h,support
import json
f=PrimeField(5);valid={product_h(f,rs) for rs in combinations_with_replacement((0,1),5)};domain={'kind':'explicit','polynomial':product_h(f,(0,1))};count=0
for cs in product(range(5),repeat=5):
 H=cs+(1,);expected=H in valid
 assert support(f,H,domain,False)==expected and support(f,H,domain,True)==expected;count+=1
print(json.dumps({'status':'passed','all_monic_quintics_GF5':count,'zero_derivative_quintics_included':True,'actualwide_runs':0}))
