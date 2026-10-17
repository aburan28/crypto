from rootfinder import polynomial as p

def product_h(f,roots):
 H=(1,)
 for a in roots:H=p.mul(H,(f.neg(a),1),f)
 return H

def prefix_parameters(f,basis):
 values=list(basis);betas=[]
 for i in range(len(values)):
  b=values[i]
  if not b:raise ValueError('dependentbinarybasis')
  betas.append(b)
  for j in range(i+1,len(values)):values[j]=f.add(f.square(values[j]),f.mul(b,values[j]))
 return betas

def remainder(f,H,domain):
 if domain['kind']=='prefix':
  r=p.divmod_poly((0,1),H,f)[1]
  for b in domain['betas']:r=p.divmod_poly(p.add(p.mul(r,r,f),tuple(f.mul(b,a) for a in r),f),H,f)[1]
  return r
 return p.divmod_poly(domain['polynomial'],H,f)[1]
def squarefree(f,H):
 derivative=tuple(f.mul(i%f.characteristic,c) for i,c in enumerate(H) if i)[::1]
 return p.gcd_poly(H,derivative,f)==(1,)
def support(f,H,domain,shortcut=False):
 r=remainder(f,H,domain)
 if shortcut and squarefree(f,H):return not r
 power=1
 if f.characteristic==2:
  while power<len(H)-1:r=p.divmod_poly(p.mul(r,r,f),H,f)[1];power*=2
 else:r=p.powmod_poly(r,len(H)-1,H,f)
 return not r
