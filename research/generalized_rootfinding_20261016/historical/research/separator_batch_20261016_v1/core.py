import time
from research.fast_curve_fibers_20261015_v1.core import FastFiber
from research.advance.general_chain import _Unresolved

class BatchSeparator:
 """Complete target-independent leftseparator, exactfinite auxiliaryfield.
 Bounded materialization; does not enumerate hard2^18source domains.
 """
 def __init__(self,field,domain,curve,*,max_domain=4096,max_left_records=100000,max_root_calls=1000000,max_seconds=30):
  if any(type(v) is not int or v<0 for v in [max_domain,max_left_records,max_root_calls]):raise ValueError('integer budgets')
  if not isinstance(max_seconds,(int,float)) or not 0<max_seconds<float('inf'):raise ValueError('time budget')
  self.started=time.monotonic();self.max_seconds=max_seconds;self.fiber=FastFiber(field,curve);self.field=field;self.domain=[];self.max_root_calls=max_root_calls;self.calls=0;self.left={};self.left_records=0;self.checked=set()
  for x in domain:
   self.check()
   if len(self.domain)>=max_domain:raise _Unresolved('source domain preflightcap')
   self.domain.append(field.validate(x))
  if len(self.domain)!=len(set(self.domain)):raise ValueError('duplicate domain')
  self.domain=tuple(sorted(self.domain))
  k=len(self.domain)
  if k>max_domain or 4*k**3>max_left_records or k*k+2*k**3>max_root_calls:raise _Unresolved('completeleftseparator preflight cap')
  for x1 in self.domain:
   for x2 in self.domain:
    for u1 in self.roots(x1,x2):
     self.certify(x1,x2,u1)
     for x3 in self.domain:
      for u2 in self.roots(u1,x3):
       self.certify(u1,x3,u2);self.left.setdefault(u2,[]).append((x1,x2,x3,u1));self.left_records+=1
       if self.left_records>max_left_records:raise ArithmeticError('left record bound')
 def check(self):
  if time.monotonic()-self.started>=self.max_seconds:raise _Unresolved('chargedbatch deadline')
 def roots(self,x,y):
  self.check()
  if self.calls>=self.max_root_calls:raise _Unresolved('rootcallcap')
  self.calls+=1;return self.fiber.roots(x,y)
 def certify(self,x,y,z):
  key=tuple(sorted((x,y,z)))
  if key not in self.checked:
   if self.fiber.evaluate(x,y,z):raise ArithmeticError('originaledgeverification')
   self.checked.add(key)
 def solve(self,target,*,max_output=100000,max_right_records=100000):
  if type(max_output) is not int or max_output<0 or type(max_right_records)is not int or max_right_records<0:raise ValueError('output/recordcap')
  self.check();self.field.validate(target);k=len(self.domain)
  if 4*k*k>max_right_records:raise _Unresolved('rightpreflightcap')
  outputs=[];count=0;right_records=0
  try:
   for x5 in self.domain:
    for u3 in self.roots(x5,target):
     self.certify(u3,x5,target)
     for x4 in self.domain:
      for u2 in self.roots(u3,x4):
       self.certify(u2,x4,u3);right_records+=1
       if right_records>max_right_records:raise ArithmeticError('rightrecordbound')
       matches=self.left.get(u2,());count+=len(matches)
       for x1,x2,x3,u1 in matches[:max(0,max_output-len(outputs))]:
        self.check();outputs.append([x1,x2,x3,x4,x5,u1,u2,u3])
   outputs.sort();self.check()
   return dict(search_complete=True,output_complete=len(outputs)==count,solution_count_exact=count,solution_count_lower_bound=count,solutions=outputs,verified=True,verification='fourcertifiedoriginalS3edgescomposedunderexactseparator')
  except _Unresolved as x:return dict(search_complete=False,output_complete=False,solution_count_exact=None,solution_count_lower_bound=count,solutions=outputs,verified=True,reason=str(x))
