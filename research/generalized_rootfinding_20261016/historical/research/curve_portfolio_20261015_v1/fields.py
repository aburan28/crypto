import time
from rootfinder.binaryfield import BinaryField,_bit_remainder,_bit_square,_bit_gcd

class BudgetedBinaryField(BinaryField):
 """Explicit-modulus binary backend without a fixed256 degree ceiling.
 Exact Rabin proof; construction is bounded by caller's validation deadline.
 No automatic modulus search or arbitrary-scale performance claim.
 """
 __slots__=()
 def __init__(self,degree,modulus,*,validation_seconds=30):
  if type(degree) is not int or degree<1 or type(modulus) is not int or modulus<2 or modulus.bit_length()!=degree+1:raise ValueError('explicit positive degree and monic polynomial required')
  if type(validation_seconds) not in (int,float) or not 0<validation_seconds<float('inf'):raise ValueError('positive finite validation budget')
  started=time.monotonic()
  def check():
   if time.monotonic()-started>=validation_seconds:raise TimeoutError('exact field validation unresolved')
  if degree>1:
   if not modulus&1 or modulus.bit_count()%2==0:raise ValueError('reducible modulus')
   left=degree;factors=set();d=2
   while d*d<=left:
    check()
    if left%d==0:
     factors.add(d)
     while left%d==0:left//=d
    d+=1
   if left>1:factors.add(left)
   checkpoints={degree//q for q in factors};power=2
   for j in range(1,degree+1):
    check();power=_bit_remainder(_bit_square(power),modulus)
    if j in checkpoints and _bit_gcd(power^2,modulus)!=1:raise ValueError('reducible modulus')
   if power!=2:raise ValueError('reducible modulus')
  self.characteristic=2;self.degree=degree;self.order=1<<degree;self.modulus=modulus
