import json,re,itertools,hashlib,importlib.util
from pathlib import Path
from sage.all import GF,PolynomialRing,EllipticCurve
base=Path(__file__).resolve().parents[3]
root=base/'research/sat_factor_base_review_20260908/continuation-05-sota-gates'
v2=json.loads((root/'stage-1-seed-20260909-v2/result.json').read_text());v3=json.loads((root/'stage-1-wdsat-rebuild-20260909-v3/result.json').read_text())
def parse_anf(text):
 lines=text.splitlines();n=int(lines[0].split()[2]);rows=[]
 assert len(lines)-1==int(lines[0].split()[3])
 for line in lines[1:]:
  t=line.split();assert t[0]=='x' and t[-1]=='0';eq={0};i=1
  while t[i]!='0':
   v=t[i];i+=1
   if v=='T':mask=0
   elif v.startswith('.'):
    count=int(v[1:]);mask=sum(1<<(int(w)-1) for w in t[i:i+count]);i+=count
   else:mask=1<<(int(v)-1)
   if mask in eq:eq.remove(mask)
   else:eq.add(mask)
  rows.append(eq)
 return n,rows

def eval_rows(rows,bits):
 assignment=sum(int(b)<<j for j,b in enumerate(bits))
 return all(sum((m&assignment)==m for m in row)%2==0 for row in rows)
def reduced(expr,n):
 values={}
 for exps,c in expr.dict().items():
  mask=sum(1<<i for i,v in enumerate(exps) if v)
  values[mask]=values.get(mask,0)+c
 out=[set() for _ in range(n)]
 for mask,c in values.items():
  for j,b in enumerate(c._vector_()):
   if b:out[j].add(mask)
 return out

def s3(a,b,c):return (a*b+a*c+b*c)**2+a*b*c+1
result=[]
for old,new in zip(v2['instances'],v3['instances']):
 assert old['cell']==new['cell'];cell=old['cell'];n=cell['n'];ell=cell['ell'];kind=cell['basis'];name=f'n{n}-l{ell}-m3-{kind}'
 if 'manifest' not in old:assert n==67 and old['generator']['returncode']!=0;continue
 m=old['manifest'];assert m['exports']==new['manifest']['exports'];a=root/'stage-1-seed-20260909-v2'/name;b=root/'stage-1-wdsat-rebuild-20260909-v3'/name
 for f in ['instance.anf','instance.xor.cnf','instance.magma']:assert (a/f).read_bytes()==(b/f).read_bytes()
 nv,rows=parse_anf((a/'instance.anf').read_text());assert nv==m['source_variables']
 fp=PolynomialRing(GF(2),'z');z=fp.gen();F=GF(2**n,'a',modulus=z**n+sum(z**i for i in m['irreducible_low_terms']));E=EllipticCurve(F,[1,1,0,0,1]);O=E(0)
 basis=[F.from_integer(int(x)) for x in m['factor_base_basis_bitmasks']]
 if kind=='standard':assert [int(v.to_integer()) for v in basis]==[1<<i for i in range(ell)]
 else:
  exps=m['factor_base_predicate']['linearised_exponents'];factor=sum(z**i for i in exps);assert (z**n+1)%factor==0 and factor.is_irreducible()
  assert all(sum((v**(2**i) for i in exps),F(0))==0 for v in basis)
 def point(d):return E(F.from_integer(int(d['x'])),F.from_integer(int(d['y'])))
 target=point(m['target']);planted=[point(p) for p in m['planted_points']];assert sum(planted,O)==target
 values=[F(0)]
 for v in basis:values += [x+v for x in values]
 assert len(set(values))==2**ell and all(p[0] in values for p in planted)
 pts=[]
 def lift(x):
  if not x:return [E(0,1)]
  c=x+1+x**-2
  if c.trace():return []
  h=sum((c**(4**i) for i in range((n+1)//2)),F(0));assert h*h+h==c
  return [E(x,x*h),E(x,x*(h+1))]
 for x in values:pts.extend(lift(x))
 assert len(pts)==m['direct_meet_in_the_middle']['factor_points']
 # Independently reconstruct the source equations, not just their serialization.
 ring=PolynomialRing(F,nv,'b');vs=ring.gens();xs=[sum((vs[i*ell+j]*v for j,v in enumerate(basis)),ring(0)) for i in range(3)]
 if kind=='ggmp':
  u=sum((vs[3*ell+j]*F.gen()**j for j in range(n)),ring(0));expected=[r for expr in [s3(xs[0],xs[1],u),s3(u,xs[2],target[0])] for r in reduced(expr,n) if r]
 else:
  sig=[sum(xs),xs[0]*xs[1]+xs[0]*xs[2]+xs[1]*xs[2],xs[0]*xs[1]*xs[2]];expected=[];offset=3*ell;es=[]
  for i in range(3):
   length=(i+1)*(ell-1)+1;coords=reduced(sig[i],n)
   for j in range(length):expected.append(coords[j]^{1<<(offset+j)})
   es.append(sum((vs[offset+j]*F.gen()**j for j in range(length)),ring(0)));offset+=length
  e1,e2,e3=es;t=target[0]
  expr=t**4+e1**4+e3**4+e2**4*t**4+e3**3*t+e3*e2**2*t**3+e3*e1**2*t+e3*t**3+e1**2*e3**2*t**2+e3**2*t**4+e3**2+e2**2*t**2
  expected += [r for r in reduced(expr,n) if r]
 assert rows==expected,(name,'source rows')
 # Independently reconstruct all AND gates and XOR constraints from source rows.
 text=(a/'instance.xor.cnf').read_text();lines=text.splitlines();maxvar=int(lines[0].split()[2]);cnfs=[];xor=[]
 for line in lines[1:]:
  t=line.split();assert t[-1]=='0'
  (xor if t[0]=='x' else cnfs).append(list(map(int,t[1:-1] if t[0]=='x' else t[:-1])))
 products=sorted({tuple(j+1 for j in range(nv) if mask>>j&1) for row in rows for mask in row if mask.bit_count()>1})
 amap={p:nv+1+i for i,p in enumerate(products)};expectedcnf=[]
 for p,g in amap.items():expectedcnf.append([-v for v in p]+[g]);expectedcnf.extend([[v,-g] for v in p])
 assert cnfs==expectedcnf and maxvar==nv+len(products) and len(lines)-1==int(lines[0].split()[3])
 xoreqs=[]
 for literals in xor:
  eq={0}
  for lit in literals:
   mask=1<<(abs(lit)-1) if abs(lit)<=nv else sum(1<<(v-1) for v in products[abs(lit)-nv-1])
   if mask in eq:eq.remove(mask)
   else:eq.add(mask)
   if lit<0:
    if 0 in eq:eq.remove(0)
    else:eq.add(0)
  xoreqs.append(eq)
 assert xoreqs==rows
 checks=[]
 for label,folder,record in [('wdsat',b,next(x for x in new['solvers'] if x['solver']=='wdsat')),('cryptominisat',a,next(x for x in old['solvers'] if x['solver']=='cryptominisat'))]:
  if record['status']!='sat':checks.append(dict(solver=label,status=record['status']));continue
  out=(folder/f'{label}.stdout').read_text()
  if label=='wdsat':bits=list(map(int,next(v.strip() for v in out.splitlines() if len(v.strip())>=nv and set(v.strip())<={'0','1'})))[:nv]
  else:
   d={abs(int(t)):int(t)>0 for line in out.splitlines() if line.startswith('v ') for t in line.split()[1:] if t!='0'};assert len(d)==maxvar;bits=[int(d[i]) for i in range(1,maxvar+1)]
   assert all(any((bits[abs(lit)-1]==1)==(lit>0) for lit in clause) for clause in cnfs)
   assert all(sum((bits[abs(lit)-1]==1)==(lit>0) for lit in literals)%2==1 for literals in xor)
  assert eval_rows(rows,bits[:nv])
  xx=[sum((basis[j]*bits[i*ell+j] for j in range(ell)),F(0)) for i in range(3)]
  choices=[lift(x) for x in xx];witness=next((points for points in itertools.product(*choices) if sum(points,O)==target),None)
  checks.append(dict(solver=label,status='sat',source_model_valid=True,rational_point_decomposition=witness is not None,xs=[int(v.to_integer()) for v in xx]))
 rec=dict(cell=name,n=n,ell=ell,predicate_basis=[int(v.to_integer()) for v in basis],factor_points=len(pts),target=m['target'],planted_points=m['planted_points'],source_rows_reconstructed=True,exports_semantically_equivalent=True,checks=checks)
 print(json.dumps(rec),flush=True);result.append(rec)
print('Stage-1 independent reconstruction complete; internal replay evidence is retained separately in INDEPENDENT_INTERNAL_CHECK_RESULT.json.',flush=True)
