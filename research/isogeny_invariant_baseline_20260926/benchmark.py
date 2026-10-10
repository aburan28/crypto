"""Exact tiny-field construction and two-summand incidence baseline; no DLP solver."""
import argparse,json,time,itertools,collections,platform
from functools import lru_cache
class Field:
 def __init__(self,q,mod):
  self.q=q;self.mod=[x%q for x in mod];self.n=len(mod)-1;self.N=q**self.n
 def digits(self,a):
  out=[]
  for _ in range(self.n):out.append(a%self.q);a//=self.q
  return out
 def enc(self,a):return sum((x%self.q)*self.q**i for i,x in enumerate(a))
 def add(self,a,b):return self.enc([x+y for x,y in zip(self.digits(a),self.digits(b))])
 def neg(self,a):return self.enc([-x for x in self.digits(a)])
 def sub(self,a,b):return self.add(a,self.neg(b))
 def mul(self,a,b):
  aa=self.digits(a);bb=self.digits(b);c=[0]*(2*self.n-1)
  for i,x in enumerate(aa):
   for j,y in enumerate(bb):c[i+j]+=x*y
  for i in range(len(c)-1,self.n-1,-1):
   v=c[i]%self.q
   for j in range(self.n):c[i-self.n+j]-=v*self.mod[j]
  return self.enc(c[:self.n])
 def pow(self,a,e):
  r=1
  while e:
   if e&1:r=self.mul(r,a)
   a=self.mul(a,a);e//=2
  return r
 def inv(self,a):
  assert a;return self.pow(a,self.N-2)
 def div(self,a,b):return self.mul(a,self.inv(b))
def orbit_hist(S,f):
 unseen=set(S);h=collections.Counter()
 while unseen:
  a=min(unseen);o={a};b=f(a)
  while b!=a:
   assert b in S and b not in o;o.add(b);b=f(b)
  unseen-=o;h[len(o)]+=1
 return dict(sorted(h.items()))
def construction(q,mod,D,a):
 F=Field(q,mod);w=q;t=time.perf_counter()
 # Rabin irreducibility criterion, with polynomial Euclidean gcd.
 def trim(v):
  while v and v[-1]==0:v.pop()
  return v
 def gcd(a,b):
  a=trim(a[:]);b=trim(b[:])
  while b:
   r=a[:]
   while len(r)>=len(b):
    shift=len(r)-len(b);c=r[-1]*pow(b[-1],-1,q)%q
    for j in range(len(b)):r[j+shift]=(r[j+shift]-c*b[j])%q
    trim(r)
   a,b=b,r
  return a
 assert F.pow(w,q**F.n)==w
 primes=[p for p in range(2,F.n+1) if F.n%p==0 and all(p%d for d in range(2,p))]
 for p in primes:assert len(gcd(F.mod,F.digits(F.sub(F.pow(w,q**(F.n//p)),w))))==1
 # Verify fiber [n](w)=a using (w+sqrt(D))**n.
 A,B=1,0
 for _ in range(F.n):A,B=F.add(F.mul(A,w),F.mul(D,B)),F.add(A,F.mul(B,w))
 assert F.div(A,B)==a
 # Enumerate projectively normalized nonzero denominators; retain redundancy counts.
 den=[1]+[F.add(c,w) for c in range(q)]
 nums=[F.add(c,F.mul(d,w)) for c in range(q) for d in range(q)]
 multiplicities=collections.Counter(F.div(f,g) for f in nums for g in den)
 S=set(multiplicities);assert len(S)==q**3
 frob={x:F.pow(x,q) for x in S};assert set(frob.values())==S
 translations=[]
 for u in range(q):
  if F.div(F.add(F.mul(u,w),D),F.add(w,u))==F.pow(w,q):translations.append(u)
 assert len(translations)==1
 out={'q':q,'n':F.n,'modulus_ascending':mod,'D':D,'fiber_a':a,'irreducible':True,'frobenius_translation_u':translations[0],'set_size':len(S),'raw_normalized_parameter_count':sum(multiplicities.values()),'representation_multiplicity_histogram':dict(collections.Counter(multiplicities.values())),'frobenius_orbits':orbit_hist(S,lambda x:frob[x]),'seconds':time.perf_counter()-t}
 return F,S,out

def tiny(ca=1,cb=1):
 F,S,out=construction(3,[1,1,0,-1,1],2,1);N=F.N
 # Cache complete tiny-field tables. These are explicitly charged as setup, not a solver speed claim.
 add=[[F.add(a,b) for b in range(N)] for a in range(N)];mul=[[F.mul(a,b) for b in range(N)] for a in range(N)];neg=[F.neg(a) for a in range(N)];inv=[0]+[F.inv(a) for a in range(1,N)]
 A=lambda a,b:add[a][b];M=lambda a,b:mul[a][b];sub=lambda a,b:A(a,neg[b]);sq=lambda a:M(a,a)
 roots=collections.defaultdict(list)
 for y in range(N):roots[sq(y)].append(y)
 # Smooth ordinary-form curve in characteristic 3: derivative of x^3+x+1 is 1.
 ys={x:roots[A(A(M(sq(x),x),M(ca,x)),cb)] for x in range(N)}
 pts=[None]+[(x,y) for x in range(N) for y in ys[x]];idx={p:i for i,p in enumerate(pts)}
 def padd(P,Q):
  if P is None:return Q
  if Q is None:return P
  x,y=P;u,v=Q
  if x==u and A(y,v)==0:return None
  slope=M(sub(v,y),inv[sub(u,x)]) if x!=u else M(ca,inv[M(2,y)])
  z=sub(sub(sq(slope),x),u)
  return z,sub(M(slope,sub(x,z)),y)
 table=[[idx[padd(P,Q)] for Q in pts] for P in pts]
 assert all(table[0][i]==i and table[i][0]==i for i in range(len(pts)))
 # Exhaustive associativity certifies the tiny group table used for counts.
 assert all(table[table[i][j]][k]==table[i][table[j][k]] for i in range(len(pts)) for j in range(len(pts)) for k in range(len(pts)))
 def f3(x,y,z):
  xy=M(x,y);ss=A(x,y)
  return sub(A(M(sq(sub(x,y)),sq(z)),sub(sq(sub(xy,ca)),M(cb,ss))),M(2,M(A(M(ss,A(xy,ca)),M(2,cb)),z))) # 4=1 mod 3
 # Independent summation polynomial verification for all liftable x and affine targets.
 valid=[x for x in range(N) if ys[x]];f3checks=0
 for x in valid:
  for y in valid:
   possible={padd((x,u),(y,v)) for u in ys[x] for v in ys[y]}
   possible_x={p[0] for p in possible if p is not None}
   for z in valid:
    assert (f3(x,y,z)==0)==(z in possible_x);f3checks+=1
 frob=lambda x:F.pow(x,3)
 def measure(name,V):
  B=[i for i,p in enumerate(pts) if p is not None and p[0] in V];counts=[0]*len(pts);distinct=[0]*len(pts)
  for i in B:
   for j in B:
    counts[table[i][j]]+=1
    if pts[i][0]!=pts[j][0]:distinct[table[i][j]]+=1
  invariant={frob(x) for x in V}==V
  result={'name':name,'field_set_size':len(V),'frobenius_invariant':invariant,'liftable_x':sum(bool(ys[x]) for x in V),'affine_base_points':len(B),'covered_affine_targets':sum(c>0 for c in counts[1:]),'covered_affine_targets_distinct_x':sum(c>0 for c in distinct[1:]),'ordered_pair_count_histogram_affine':dict(sorted(collections.Counter(counts[1:]).items())),'counts_in_point_order':counts,'distinct_x_counts_in_point_order':distinct}
  if invariant:
   result['field_orbits']=orbit_hist(V,frob)
   result['point_orbits']=orbit_hist(set(B),lambda i:idx[(frob(pts[i][0]),frob(pts[i][1]))])
  return result
 comparisons=[]
 # All hyperplanes, indexed by normalized nonzero linear functionals.
 for c in itertools.product(range(3),repeat=4):
  if not any(c) or next(x for x in c if x)!=1:continue
  V={x for x in range(N) if sum(a*b for a,b in zip(c,F.digits(x)))%3==0}
  r=measure('kernel_'+''.join(map(str,c)),V);r['functional']=c;comparisons.append(r)
 assert len(comparisons)==40
 fraction=measure('torus_fractions_k1',S)
 witness=next((a,b,A(a,b)) for a in sorted(S) for b in sorted(S) if A(a,b) not in S)
 out.update({'curve':f'y^2=x^3+{ca}x+{cb}','curve_order':len(pts),'affine_target_count':len(pts)-1,'nonadditivity_witness':witness,'f3_equivalence_checks':f3checks,'fraction':fraction,'all_40_hyperplanes':comparisons,'points':pts})
 return out
if __name__=='__main__':
 parser=argparse.ArgumentParser(description='Exact eligible torus and exhaustive tiny-curve incidence benchmark')
 parser.add_argument('--output',required=True,help='New raw JSON result path')
 args=parser.parse_args()
 start=time.perf_counter();small=tiny()
 # Published Couveignes--Lercier q=13,d=7 example.
 mod=[-64,4,-48,10,-40,3,-56,1]
 F,S,large=construction(13,mod,2,8)
 assert F.pow(13,13)==F.div(F.add(F.mul(4,13),2),F.add(13,4))
 sweep=[tiny(a,b) for a in [1,2] for b in range(3) if (a,b)!=(1,1)]
 result={'other_five_curves':sweep,'scope':'Exact torus construction and tiny two-summand incidence; no polynomial solver timing or DLP recovery','tiny':small,'published_13_7':large,'python':platform.python_version(),'total_seconds':time.perf_counter()-start}
 with open(args.output,'x') as f:json.dump(result,f,indent=2)
 print(f"Wrote {args.output}: 6 tiny curves, {sum(c['f3_equivalence_checks'] for c in [small]+sweep)} lift checks; construction size {large['set_size']} over F_13^7")
