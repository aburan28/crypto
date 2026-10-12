"""Exact tiny-group controls for the independent dimension/length sweep."""
import importlib.util
import math
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'coefficient_pullback'))
import experiment as prior
import curves
import nagaoannihilator as scalar

spec=importlib.util.spec_from_file_location('sweep_cold_reference',HERE.parent/'coefficient_pullback/e2e.py')
cold=importlib.util.module_from_spec(spec);spec.loader.exec_module(cold)


class Curve(curves.Curve):
    def __init__(self,f,b):
        super().__init__(f);self.b=b

    def onCurve(self,p):
        if p is None:return True
        x,y=p;f=self.f
        return f.add(f.mul(y,y),f.mul(x,y))==f.add(f.mul(f.mul(x,x),x),self.b)

    def pointFromX(self,x):
        if not x:return None
        f=self.f;inv=f.inv(f.mul(x,x))
        c=f.add(x,inv if self.b==f.one() else f.mul(self.b,inv))
        if f.trace(c):return None
        p=(x,f.mul(x,self.halfTrace(c)))
        if not self.onCurve(p):raise ArithmeticError('invalid point construction')
        return p


def fieldCurve(profile):
    f=scalar.CountedField(profile['bits'])
    return f,Curve(f,f.fromCoords(profile['curve_b_coords']))


def support(f,profile,ell):
    if profile['q_bits']==1:
        return [f.fromCoords(x) for x in range(1<<ell)]
    qb=profile['q_bits']
    subfield=[f.fromCoords(x) for x in range(1<<f.m)
              if f.frob(f.fromCoords(x),qb)==f.fromCoords(x)]
    if len(subfield)!=1<<qb:raise ArithmeticError('wrong embedded subfield')
    span={0};basis=[]
    candidates=[f.one()]+[f.fromCoords(x) for x in range(1,1<<f.m) if f.fromCoords(x)!=f.one()]
    for candidate in candidates:
        if candidate in span:continue
        multiples=[f.mul(candidate,c) for c in subfield]
        span={f.add(a,b) for a in span for b in multiples}
        basis.append(candidate)
        if len(basis)==ell:break
    if len(span)!=(1<<qb)**ell:raise ArithmeticError('support dimension mismatch')
    return sorted(span,key=f.toCoords)


def factors(curve,values):
    return [p for x in values if x and (p:=curve.pointFromX(x)) is not None]


def primeFactors(value):
    out=[];p=2
    while p*p<=value:
        if value%p==0:
            out.append(p)
            while value%p==0:value//=p
        p+=1
    if value>1:out.append(value)
    return out


def groupSetup(profile):
    f,curve=fieldCurve(profile)
    zero=(0,f.frob(curve.b,f.m-1))
    points=[zero]
    for x in range(1,1<<f.m):
        p=curve.pointFromX(f.fromCoords(x))
        if p is not None:points.extend((p,curve.neg(p)))
    full=len(points)+1;primes=primeFactors(full);order=max(primes);h=full//order
    if math.gcd(h,order)!=1:raise ValueError('cofactor projection is not invertible modulo selected prime')
    A=next((p for p in points if curve.mul(p,full) is None
            and all(curve.mul(p,full//r) is not None for r in primes)),None)
    if A is None:raise ValueError('cyclic group generator unavailable')
    G=curve.mul(A,h)
    if G is None or curve.mul(G,order) is not None:raise ArithmeticError('invalid prime-subgroup generator')
    return f,curve,points,full,order,h,A,G


class Truth:
    """Validation-only cyclic coordinates. Never passed to a measured solver."""
    def __init__(self,profile):
        self.profile=profile
        self.f,self.curve,self.points,self.full,self.order,self.h,self.A,self.G=groupSetup(profile)
        self.logs={None:0};p=None
        for i in range(1,self.full):
            p=self.curve.add(p,self.A)
            if p in self.logs:raise ArithmeticError('premature group cycle')
            self.logs[p]=i
        if self.curve.add(p,self.A) is not None or set(self.logs)!={None,*self.points}:
            raise ArithmeticError('group coordinate oracle incomplete')
        self.cache={}

    def census(self,ell,k):
        key=(ell,k)
        if key in self.cache:return self.cache[key]
        f=self.f;curve=self.curve;values=support(f,self.profile,ell);base=factors(curve,values)
        exponents=[self.logs[p] for p in base];which={p[0]:i for i,p in enumerate(base)}
        def dp(excluded):
            rows=[[0]*self.full for _ in range(k+1)];rows[0][0]=1
            for i,exponent in enumerate(exponents):
                if i==excluded:continue
                for length in range(k,0,-1):
                    for value,count in enumerate(rows[length-1]):
                        if count:
                            rows[length][(value+exponent)%self.full]+=count
                            rows[length][(value-exponent)%self.full]+=count
            return rows[k]
        countVectors={None:dp(None)}
        counts={}
        for target in self.points:
            excluded=which.get(target[0])
            if excluded not in countVectors:countVectors[excluded]=dp(excluded)
            counts[target]=countVectors[excluded][self.logs[target]]
        projected=set()
        for p in base:
            p=curve.mul(p,self.h)
            if p is not None:projected.add(min(p,curve.neg(p)))
        M=len(base);q=1<<self.profile['q_bits'];n=self.profile['bits']//self.profile['q_bits']
        lam=q**(k*ell-n)/math.factorial(k)
        mean=sum(counts.values())/(self.full-1)
        output={'profile':self.profile['name'],'bits':self.profile['bits'],'q':q,'n':n,'ell':ell,'k':k,
                'curve_order':self.full,'subgroup_order':self.order,'cofactor':self.h,
                'support_size':len(values),'raw_abscissae':M,'raw_signed_points':2*M,'projected_columns':len(projected),
                'affine_targets':self.full-1,'covered_targets':sum(v>0 for v in counts.values()),
                'exact_uniform_coverage':sum(v>0 for v in counts.values())/(self.full-1),
                'exact_mean_representations':mean,'paper_intensity':lam,'paper_sparse_probability_capped':min(1,lam),
                'support_count_ceiling':min(1,2**k*math.comb(M,k)/(self.full-1)) if M>=k else 0,
                'poisson_from_actual_mean':-math.expm1(-mean)}
        self.cache[key]=(output,counts);return output,counts


class Search:
    def __init__(self,f,curve,profile,ell,k,variant,deadline):
        self.f=f;self.curve=curve;self.k=k;self.variant=variant;self.deadline=deadline
        self.stats={'prefix_visits':0,'table_entries':0,'lookups':0,'bucket_checks':0}
        f.phase='factor_base';self.values=support(f,profile,ell);self.base=factors(curve,self.values)
        self.signed=[[(i,1,p),(i,-1,curve.neg(p))] for i,p in enumerate(self.base)]
        self.singles=[(term[2],(term,)) for terms in self.signed for term in terms]
        self.lookup={p:terms for p,terms in self.singles}
        f.phase='table_setup';self.table={};self.left=[]
        if variant=='pair-mitm':
            if k==2:entries=self.singles
            else:
                entries=[]
                for i in range(len(self.base)):
                    for j in range(i+1,len(self.base)):
                        self.check()
                        for a in self.signed[i]:
                            for b in self.signed[j]:entries.append((curve.add(a[2],b[2]),(a,b)))
            for point,terms in entries:self.table.setdefault(point,[]).append(terms)
            self.stats['table_entries']=len(entries)
            self.left=entries if k==4 else self.singles
        elif variant!='enumerate-last':raise ValueError('unknown solver')
        self.check()

    def check(self):
        if time.perf_counter()>=self.deadline:raise TimeoutError('cooperative deadline')

    def verify(self,target,terms):
        f=self.f;f.phase='verification';total=None
        xs=[term[2][0] for term in terms]
        if len(terms)!=self.k or len(set(xs))!=self.k or any(not x or x==target[0] or x not in self.values for x in xs):
            raise ArithmeticError('invalid decomposition domain')
        for _,_,p in terms:
            if not self.curve.onCurve(p):raise ArithmeticError('off-curve factor')
            total=self.curve.add(total,p)
        if total!=target:raise ArithmeticError('wrong point sum')

    def find(self,target):
        self.check();self.f.phase='decomposition'
        if len(self.base)<self.k:return None
        curve=self.curve
        if self.variant=='pair-mitm':
            for point,left in self.left:
                self.check();self.stats['prefix_visits']+=1
                if any(term[2][0]==target[0] for term in left):continue
                needed=curve.add(target,curve.neg(point));self.stats['lookups']+=1
                for right in self.table.get(needed,[]):
                    self.stats['bucket_checks']+=1
                    if left[-1][0]>=right[0][0] or any(term[2][0]==target[0] for term in right):continue
                    terms=left+right;self.verify(target,terms);return terms
            return None
        def walk(start,partial,terms):
            self.check();self.stats['prefix_visits']+=1
            if len(terms)==self.k-1:
                needed=curve.add(target,curve.neg(partial));self.stats['lookups']+=1
                last=self.lookup.get(needed)
                if last is not None and last[0][0]>=start and last[0][2][0]!=target[0]:return terms+last
                return None
            for i in range(start,len(self.base)):
                if self.base[i][0]==target[0]:continue
                for term in self.signed[i]:
                    found=walk(i+1,curve.add(partial,term[2]),terms+(term,))
                    if found is not None:return found
            return None
        answer=walk(0,None,())
        if answer is not None:self.verify(target,answer)
        return answer


def certificate(f,terms):
    return None if terms is None else {'x':[f.toCoords(t[2][0]) for t in terms],
        'signs':[t[1] for t in terms],'points':[[f.toCoords(v) for v in t[2]] for t in terms]}


def stageCell(profile,ell,k,targetCoords,variant,budget):
    started=time.perf_counter();f,curve=fieldCurve(profile);search=None;answer=None
    try:
        search=Search(f,curve,profile,ell,k,variant,started+budget)
        answer=search.find(tuple(f.fromCoords(v) for v in targetCoords));status='found' if answer else 'empty'
    except TimeoutError:status='timeout'
    elapsed=time.perf_counter()-started
    return {'profile':profile['name'],'ell':ell,'k':k,'variant':variant,'target':targetCoords,'status':status,
            'witness':certificate(f,answer),'cold_seconds':elapsed,'within_budget':elapsed<=budget,
            'field_api_counts':f.report(),'search_counters':None if search is None else search.stats,
            'full_dlp_S':None,'rho_ratio':None}


def dlp(profile,ell,k,variant,seed,contract):
    started=time.perf_counter();deadline=started+contract['run_seconds']
    f,curve,_,full,N,h,A,G=groupSetup(profile)
    secret=random.Random(seed+999).randrange(1,N);Q=curve.mul(G,secret);ops=cold.ScalarOps(N)
    phases={'group_setup':time.perf_counter()-started,'factor_table_setup':0.,'collection':0.,'matrix':0.,'individual':0.,'rho':0.,'final_verification':0.}
    attempts=[];search=None
    output={'profile':profile['name'],'ell':ell,'k':k,'variant':variant,'seed':seed,'status':'incomplete',
            'target':[f.toCoords(v) for v in Q],'generator':[f.toCoords(v) for v in G],
            'curve_order':full,'subgroup_order':N,'cofactor':h,'attempts':attempts,'full_dlp_S':None,'rho_ratio':None}
    try:
        if variant=='rho':
            t=time.perf_counter();f.phase='rho';recovered,iters=cold.rho(curve,f,G,Q,ops,seed,deadline)
            phases['rho']+=time.perf_counter()-t;output['rho_iterations']=iters
        else:
            t=time.perf_counter();search=Search(f,curve,profile,ell,k,variant,deadline)
            columns={};mapping=[];f.phase='projection'
            for p in search.base:
                projected=curve.mul(p,h)
                if projected is None:mapping.append(None);continue
                key=min(projected,curve.neg(projected))
                if key not in columns:columns[key]=len(columns)
                mapping.append((columns[key],1 if projected==key else -1))
            phases['factor_table_setup']+=time.perf_counter()-t;output['columns']=len(columns)
            if not columns:raise ValueError('zero projected factor columns')
            if len(search.base)<k:raise ValueError('fewer admissible abscissae than summands')
            matrix=cold.Rows(len(columns),ops);rng=random.Random(seed)
            def attempt(target,phase,s):
                at=time.perf_counter()
                try:
                    answer=search.find(target);status='found' if answer else 'empty'
                except TimeoutError:
                    attempts.append({'phase':phase,'coefficient':s,'target':[f.toCoords(v) for v in target],
                                     'status':'timeout','witness':None,'seconds':time.perf_counter()-at})
                    raise
                attempts.append({'phase':phase,'coefficient':s,'target':[f.toCoords(v) for v in target],
                                 'status':status,'witness':certificate(f,answer),'seconds':time.perf_counter()-at})
                if answer is None:return None
                row=[0]*len(columns)
                for i,sign,_ in answer:
                    if mapping[i] is not None:
                        col,orientation=mapping[i];row[col]=ops.add(row[col],sign*orientation)
                return row
            for _ in range(contract['target_attempt_cap']):
                t=time.perf_counter();f.phase='relation_generation';s=rng.randrange(1,full)
                row=attempt(curve.mul(A,s),'collection',s);phases['collection']+=time.perf_counter()-t
                if row is None:continue
                t=time.perf_counter();matrix.add(row+[s%N]);phases['matrix']+=time.perf_counter()-t
                if len(matrix.pivots)==len(columns):break
            output['rank']=len(matrix.pivots)
            t=time.perf_counter();logs=matrix.solve();phases['matrix']+=time.perf_counter()-t
            recovered=None;rng=random.Random(seed+100000);t=time.perf_counter()
            for _ in range(contract['target_attempt_cap']):
                f.phase='individual_generation';s=rng.randrange(full);target=curve.add(Q,curve.mul(A,s))
                if target is None:continue
                row=attempt(target,'individual',s)
                if row is None:continue
                value=0
                for coefficient,log in zip(row,logs):value=ops.add(value,ops.mul(coefficient,log))
                recovered=ops.mul(ops.add(value,-s),ops.inv(h%N));break
            phases['individual']+=time.perf_counter()-t
            if recovered is None:raise TimeoutError('individual attempt cap')
        t=time.perf_counter();f.phase='final_verification'
        if recovered!=secret or curve.mul(G,recovered)!=Q:raise ArithmeticError('incorrect recovered scalar')
        phases['final_verification']+=time.perf_counter()-t
        output.update(status='verified',recovered_scalar=recovered,planted_scalar=secret)
    except (TimeoutError,ValueError) as exc:output['error']=str(exc)
    elapsed=time.perf_counter()-started
    phases['unassigned_overhead_or_interrupted_phase']=max(0.,elapsed-sum(phases.values()))
    output.update(cold_seconds=elapsed,phase_seconds=phases,field_api_counts=f.report(),scalar_modular_counts=ops.counts,
                  search_counters=None if search is None else search.stats)
    return output
