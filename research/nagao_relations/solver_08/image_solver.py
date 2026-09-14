"""Hybrid support through charged image-space membership and linear root recovery."""
import time
import curves
import nagaoannihilator as scalar


def batchInverse(f,values):
    if not values:return []
    products=[];acc=f.one()
    for value in values:
        if not value:raise ZeroDivisionError('batch inverse of zero')
        products.append(acc);acc=f.mul(acc,value)
    acc=f.inv(acc);out=[None]*len(values)
    for i in range(len(values)-1,-1,-1):
        out[i]=f.mul(acc,products[i]);acc=f.mul(acc,values[i])
    return out


def asRoots(f,curve,c):
    if f.trace(c):return []
    w=curve.halfTrace(c)
    return [w,f.add(w,f.one())]


class Search:
    def __init__(self,f,curve,d,target,deadline):
        self.f=f;self.curve=curve;self.d=d;self.target=target;self.deadline=deadline
        self.images={}
        for uBits in range(1,1<<d):
            self.checkDeadline();u=f.fromCoords(uBits);pivots={}
            for i in range(d):
                pre=f.fromCoords(1<<i);value=f.add(f.sqr(pre),f.mul(u,pre))
                for k,(basis,lift) in sorted(pivots.items(),reverse=True):
                    if f.toCoords(value)>>k&1:
                        value=f.add(value,basis);pre=f.add(pre,lift)
                if value:pivots[f.toCoords(value).bit_length()-1]=(value,pre)
            if len(pivots)!=d-1:raise ArithmeticError('image rank mismatch')
            self.images[uBits]=sorted(pivots.items(),reverse=True)
        r,s=target
        self.zs=[f.fromCoords(i) for i in range(1,1<<d) if f.fromCoords(i)!=r]
        values=self.zs+[f.add(r,z) for z in self.zs]
        invs=batchInverse(f,values);count=len(self.zs)
        # alpha=z/(r+z), beta=(r+z)/z².  The AS input is K*beta/b².
        self.geometry=[(z,f.mul(z,invs[count+i]),f.mul(f.add(r,z),f.sqr(invs[i]))) for i,z in enumerate(self.zs)]
        self.checkDeadline()

    def checkDeadline(self):
        if time.perf_counter()>=self.deadline:raise TimeoutError('all-phase budget exhausted')

    def candidates(self):
        f=self.f;r,s=self.target
        for hBits in range(1<<self.d):
            self.checkDeadline();h2=f.fromCoords(hBits)
            for b in asRoots(f,self.curve,f.add(h2,r)):
                if not b:continue
                h,_=scalar.residualNorm(f,self.target,0,b)
                invB=f.inv(b);invB2=f.sqr(invB)
                for z,alpha,beta in self.geometry:
                    self.checkDeadline()
                    if z==h2:continue
                    k=scalar.polyEval(f,h,z)
                    t=f.mul(b,alpha)
                    for w in asRoots(f,self.curve,f.mul(f.mul(k,beta),invB2)):
                        yield f.mul(t,w),b,invB,z

    def recover(self,a,b,invB,z):
        f=self.f;r,s=self.target
        h,c=scalar.residualNorm(f,self.target,a,b)
        if scalar.polyEval(f,h,z):raise ArithmeticError('conditioned root missing')
        if not h[0] or not scalar.polyEval(f,h,r):return None
        u=f.add(h[2],z);uBits=f.toCoords(u)
        if uBits not in self.images:raise ArithmeticError('invalid root sum')
        v=f.add(h[1],f.mul(z,u))
        remainder=v;w=0
        for k,(basis,lift) in self.images[uBits]:
            if f.toCoords(remainder)>>k&1:
                remainder=f.add(remainder,basis);w=f.add(w,lift)
        if remainder:return None
        # Q splits inside V; the linear solve also supplies one preimage w.
        xs=[z,w,f.add(w,u)]
        coords=tuple(sorted(f.toCoords(x) for x in xs))
        if len(set(coords))!=3:return None
        if any(x==0 or x>=1<<self.d or x==f.toCoords(r) for x in coords):raise ArithmeticError('support exclusion mismatch')
        total=None
        for x in xs:
            y=f.mul(f.add(f.add(f.sqr(x),f.mul(a,x)),c),invB)
            if not self.curve.onCurve((x,y)):raise ArithmeticError('off-curve recovery')
            total=self.curve.add(total,(x,y))
        if total!=self.target:raise ArithmeticError('signed sum mismatch')
        return coords


def cell(n,d,targetCoords,mode,budget):
    start=time.perf_counter();deadline=start+budget
    f=scalar.CountedField(n);curve=curves.Curve(f);target=tuple(f.fromCoords(x) for x in targetCoords)
    solutions=set();seen=set();duplicates=0;candidates=0;first=None;complete=False
    try:
        f.phase='setup';search=Search(f,curve,d,target,deadline)
        f.phase='search'
        for a,b,invB,z in search.candidates():
            candidates+=1
            if (a,b) in seen:duplicates+=1;continue
            seen.add((a,b));f.phase='support_extract_verify';xs=search.recover(a,b,invB,z);f.phase='search'
            if xs is not None:
                solutions.add(xs)
                if first is None:first=time.perf_counter()-start
                if mode=='first':break
        else:complete=True
    except TimeoutError:pass
    elapsed=time.perf_counter()-start
    status='first' if mode=='first' and solutions else ('complete' if complete else 'timeout')
    return {'n':n,'d':d,'target':targetCoords,'variant':'quadratic-image','mode':mode,'status':status,
        'solutions':[list(x) for x in sorted(solutions)],'verified_unique_relations':len(solutions),
        'candidate_functions':candidates,'duplicate_functions':duplicates,'first_verified_seconds':first,
        'all_phase_seconds':elapsed,'within_budget':elapsed<=budget,'field_api_counts':f.report(),
        'full_dlp_S':None,'rho_ratio':None}
