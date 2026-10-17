from collections import OrderedDict
from research.advance.general_chain import GeneralS3Fiber,_Unresolved
from research.advance.finite_fiber import BinaryQuadraticFiber
from rootfinder.binaryfield import BinaryField
class Context:
 def __init__(self,field,curve,capacity=4096):
  if type(capacity) is not int or capacity<0:raise ValueError('capacity')
  self.field=field;self.curve=dict(curve);self.capacity=capacity;self.cache=OrderedDict();self.binary=None;self.nonresidue=None;self.hits=0;self.misses=0
 def power(self,a,n):return pow(a,n,self.field.characteristic) if self.field.degree==1 else self.field.pow(a,n)
 def sqrt(self,a):
  f=self.field;q=f.order
  if a==0:return 0
  if q%4==3:
   x=self.power(a,(q+1)//4);return x if f.square(x)==a else None
  if self.power(a,(q-1)//2)!=1:return None
  odd=q-1;s=0
  while odd%2==0:odd//=2;s+=1
  if self.nonresidue is None:
   for z in range(2,min(q,130)):
    if self.power(z,(q-1)//2)==f.neg(1):self.nonresidue=z;break
   if self.nonresidue is None:raise _Unresolved('bounded quadratic nonresidue search unresolved')
  c=self.power(self.nonresidue,odd);x=self.power(a,(odd+1)//2);t=self.power(a,odd);m=s
  while t!=1:
   tt=t;i=0
   while tt!=1 and i<m:tt=f.square(tt);i+=1
   if i==m:raise ArithmeticError('Tonelli invariant')
   b=self.power(c,1<<(m-i-1));x=f.mul(x,b);c=f.square(b);t=f.mul(t,c);m=i
  if f.square(x)!=a:raise ArithmeticError('sqrt replay')
  return x
class FastFiber(GeneralS3Fiber):
 def __init__(self,field,curve=None,*,context=None,**kw):
  super().__init__(field,curve,**kw);self.context=context or Context(field,self.curve)
  if self.context.field is not field or self.context.curve!=self.curve:raise ValueError('context binding')
 def roots(self,x,y):
  f=self.field;x=f.validate(x);y=f.validate(y);key=(min(x,y),max(x,y));ctx=self.context
  if key in ctx.cache:ctx.hits+=1;ctx.cache.move_to_end(key);return ctx.cache[key]
  ctx.misses+=1;c,b,a=self.coefficients(x,y)
  if f.characteristic==2 and isinstance(f,BinaryField):
   if ctx.binary is None:ctx.binary=BinaryQuadraticFiber(f)
   roots=ctx.binary.roots(a,b,c)
  elif f.characteristic==2:roots=super().roots(x,y)
  elif a==0:
   if b==0:
    if c==0:raise ArithmeticError('unrestricted fiber')
    roots=()
   else:roots=(f.mul(f.neg(c),f.inv(b)),)
  else:
   disc=f.sub(f.square(b),f.mul(4%f.characteristic,f.mul(a,c)));d=ctx.sqrt(disc)
   if d is None:roots=()
   else:
    inv=f.inv(f.mul(2%f.characteristic,a));roots=tuple(sorted({f.mul(f.add(f.neg(b),d),inv),f.mul(f.sub(f.neg(b),d),inv)}))
  if any(self.evaluate(x,y,z)!=0 for z in roots):raise ArithmeticError('originalS3 replay')
  roots=tuple(sorted(set(roots)))
  if ctx.capacity:
   ctx.cache[key]=roots
   if len(ctx.cache)>ctx.capacity:ctx.cache.popitem(last=False)
  return roots
