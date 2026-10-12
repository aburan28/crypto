"""General-coefficient successor to solver_08; frozen predecessor unchanged."""
import time
import algebra
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


class Search:
    def __init__(self,f,curve,d,target,deadline):
        self.f=f;self.curve=curve;self.d=d;self.target=target;self.deadline=deadline
        self.images={}
        for uBits in range(1,1<<d):
            self.checkDeadline();u=f.fromCoords(uBits)
            cols=[]
            for i in range(d):
                pre=f.fromCoords(1<<i)
                cols.append((f.add(f.sqr(pre),f.mul(u,pre)),pre))
            image=algebra.LinearImage(f,cols)
            if len(image.rows)!=d-1:raise ArithmeticError('image rank mismatch')
            self.images[uBits]=image
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
            for b in self.curve.asRoots(f.add(h2,r)):
                if not b:continue
                h,_=algebra.residualNorm(f,self.curve,self.target,0,b)
                invB=f.inv(b);invB2=f.sqr(invB)
                for z,alpha,beta in self.geometry:
                    self.checkDeadline()
                    if z==h2:continue
                    k=scalar.polyEval(f,h,z)
                    t=f.mul(b,alpha)
                    for w in self.curve.asRoots(f.mul(f.mul(k,beta),invB2)):
                        yield f.mul(t,w),b,invB,z

    def recover(self,a,b,invB,z):
        f=self.f;r,s=self.target
        h,c=algebra.residualNorm(f,self.curve,self.target,a,b)
        if scalar.polyEval(f,h,z):raise ArithmeticError('conditioned root missing')
        if not h[0] or not scalar.polyEval(f,h,r):return None
        u=f.add(h[2],z);uBits=f.toCoords(u)
        if uBits not in self.images:raise ArithmeticError('invalid root sum')
        v=f.add(h[1],f.mul(z,u))
        w=self.images[uBits].preimage(v)
        if w is None:return None
        # Q splits inside V; the linear solve also supplies one preimage w.
        xs=[z,w,f.add(w,u)]
        coords=tuple(sorted(f.toCoords(x) for x in xs))
        if len(set(coords))!=3:return None
        if any(x==0 or x>=1<<self.d or x==f.toCoords(r) for x in coords):raise ArithmeticError('support exclusion mismatch')
        total=None;points=[]
        for x in xs:
            y=f.mul(f.add(f.add(f.sqr(x),f.mul(a,x)),c),invB)
            if not self.curve.onCurve((x,y)):raise ArithmeticError('off-curve recovery')
            points.append([f.toCoords(x),f.toCoords(y)])
            total=self.curve.add(total,(x,y))
        if total!=self.target:raise ArithmeticError('signed sum mismatch')
        return coords,points

