from math import lcm,prod
ELL=680564733841876926932320129493409985129
S=196511074115861092422032515080945363956
FACTORS=[2,3,11,109,131,263,32326729,21234899465981031419669]
assert prod([2**3,*FACTORS[1:]])==ELL-1
orders=[]
for j in (3,4):
 a=(1+pow(S,j,ELL))%ELL
 order=ELL-1
 for q in FACTORS:
  while order%q==0 and pow(a,order//q,ELL)==1:order//=q
 orders.append(order)
 print({'jump':j,'multiplier':a,'order':order,'index':(ELL-1)//order})
assert lcm(*orders)==ELL-1
print({'joint_order':lcm(*orders),'full_multiplicative_group':True})
