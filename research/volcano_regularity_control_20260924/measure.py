"""Fixed F_103 educational control. No large-field or discrete-log solver."""
from itertools import product
from math import comb
from pathlib import Path
import json
import random

P = 103
ROOT = (0, 1)
CHILD = (88, 22)
O = None


def add(A, B, curve):
    if A is None: return B
    if B is None: return A
    x, y = A; u, v = B
    if x == u and (y + v) % P == 0: return O
    s = ((3*x*x + curve[0])*pow(2*y, -1, P) if A == B else (v-y)*pow(u-x, -1, P)) % P
    z = (s*s-x-u) % P
    return z, (s*(x-z)-y) % P


def points(curve):
    a,b = curve
    return [O]+[(x,y) for x,y in product(range(P), repeat=2) if (y*y-x*x*x-a*x-b)%P == 0]


def subgroup(G):
    out = [O]; cur = G
    while cur is not None:
        assert cur not in out
        out.append(cur);cur=add(cur,G,ROOT)
    return out


def phi(Q):
    if Q is None or Q == (102,0): return O
    x,y=Q;inv=pow((x+1)%P,-1,P)
    return (x+3*inv)%P, y*(1-3*inv*inv)%P


def univariate_roots(xs):
    c=[1]
    for x in xs:
        out=[0]*(len(c)+1)
        for i,v in enumerate(c):out[i]=(out[i]-x*v)%P;out[i+1]=(out[i+1]+v)%P
        c=out
    return c


def clean(f):return {k:v%P for k,v in f.items() if v%P}


def s3(a,b,t):
    # S3(X,Y,t), in the same two coordinate variables for every curve.
    return clean({(2,2):1,(2,1):-2*t,(1,2):-2*t,(2,0):t*t,(0,2):t*t,
                  (1,1):-2*t*t-2*a,(1,0):-2*t*a-4*b,(0,1):-2*t*a-4*b,
                  (0,0):a*a-4*b*t})


def eval_poly(f,x,y):return sum(c*pow(x,i,P)*pow(y,j,P) for (i,j),c in f.items())%P


def reduce_support(f,c):
    n=len(c)-1; f=dict(f)
    for axis in (0,1):
        while any(k[axis]>=n for k in f):
            key=max((k for k in f if k[axis]>=n),key=lambda k:k[axis])
            val=f.pop(key)
            for i in range(n):
                dst=list(key);dst[axis]=key[axis]-n+i;dst=tuple(dst)
                f[dst]=(f.get(dst,0)-val*c[i])%P
            f=clean(f)
    return f


def rank(rows,width):
    piv={}
    for row in rows:
        row=row[:]
        for i in range(width):
            if not row[i]:continue
            if i in piv:
                v=row[i];row=[(a-v*b)%P for a,b in zip(row,piv[i])]
            else:
                inv=pow(row[i],-1,P);piv[i]=[a*inv%P for a in row];break
    return len(piv)


def regularity(fs):
    top=[]
    for f in fs:
        if f:
            deg=max(sum(k) for k in f);top.append((deg,{k:v for k,v in f.items() if sum(k)==deg}))
    ranks=[]
    for d in range(10):
        rows=[]
        for deg,f in top:
            for shift in range(d-deg+1):
                row=[0]*(d+1)
                for (i,j),v in f.items():row[i+shift]=v
                rows.append(row)
        r=rank(rows,d+1);ranks.append({'degree':d,'rank':r,'dimension':d+1})
        if r==d+1:return d,ranks
    raise AssertionError('Pure support powers guarantee saturation below 10')


def run():
    root_points=points(ROOT);child_points=set(points(CHILD))
    assert len(root_points)==len(child_points)==84
    group=next(g for Q in root_points[1:] if len(g:=subgroup(Q))==21)
    image=[phi(Q) for Q in group]
    assert len(set(image))==21 and set(image)<=child_points
    for A,B in product(group,repeat=2):assert phi(add(A,B,ROOT))==add(phi(A),phi(B),CHILD)
    xs=sorted({Q[0] for Q in group if Q is not None}); assert len(xs)==10
    targets=[next(Q for Q in group if Q is not None and Q[0]==x) for x in xs]
    records=[];censuses=[]
    for size,seed in product((2,3,4),range(8)):
        support=sorted(random.Random(seed).sample(xs,size))
        base=[Q for Q in group if Q is not None and Q[0] in support]
        childbase=[phi(Q) for Q in base]
        childxs=sorted({Q[0] for Q in childbase});assert len(childxs)==size
        counts=[]
        for T in group:
            n1=sum(add(A,B,ROOT)==T for A,B in product(base,repeat=2))
            n2=sum(add(A,B,CHILD)==phi(T) for A,B in product(childbase,repeat=2))
            assert n1==n2;counts.append(n1)
        assert sum(counts)==len(base)**2
        censuses.append({'size_x':size,'seed':seed,'source_x':support,'target_counts':counts})
        for T in targets:
            record={'size_x':size,'points':2*size,'seed':seed,'target':T,'image_target':phi(T)}
            for name,curve,sx,R,B in [('source',ROOT,support,T,base),('neighbor',CHILD,childxs,phi(T),childbase)]:
                c=univariate_roots(sx)
                field=[{(i,0):v for i,v in enumerate(c) if v},{(0,i):v for i,v in enumerate(c) if v}]
                f=s3(*curve,R[0]); reduced=reduce_support(f,c)
                roots={(x,y) for x,y in product(sx,repeat=2) if eval_poly(f,x,y)==0}
                lifted={(A[0],BB[0]) for A,BB in product(B,repeat=2) if add(A,BB,curve)==R}
                assert roots==lifted
                for x,y in product(sx,repeat=2):assert eval_poly(f,x,y)==eval_poly(reduced,x,y)
                raw,raw_ranks=regularity(field+[f]);red,red_ranks=regularity(field+[reduced])
                assert raw==size+1
                record[name]={'support':sx,'raw_d_reg_top':raw,'reduced_d_reg_top':red,
                    'raw_rank_certificate':raw_ranks,'reduced_rank_certificate':red_ranks,
                    's3':[[*k,v] for k,v in sorted(f.items())],
                    'reduced_s3':[[*k,v] for k,v in sorted(reduced.items())],
                    'ordered_x_roots':len(roots),'has_relation':bool(roots)}
            records.append(record)
    report={}
    for encoding in ['raw','reduced']:
        counts={'lower':0,'equal':0,'higher':0}
        for r in records:
            delta=r['neighbor'][encoding+'_d_reg_top']-r['source'][encoding+'_d_reg_top']
            counts['lower' if delta<0 else 'higher' if delta>0 else 'equal']+=1
        report[encoding]=counts
    changed=sum(r[s]['raw_d_reg_top']!=r[s]['reduced_d_reg_top'] for r in records for s in ['source','neighbor'])
    result={'field':103,'curves':{'source':ROOT,'neighbor':CHILD},'subgroup_order':21,
       'paired_records':len(records),'unique_supports':len({tuple(c['source_x']) for c in censuses}),
       'transported_count_comparisons':21*len(censuses),'comparisons':report,
       'within_curve_preprocessing_changes':changed,'records':records,'censuses':censuses,
       'f4_f5_solving_degree':None,'n53_n83_neighbor_measurements':None,'ecdLP_speedup':None}
    Path(__file__).with_name('results.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','censuses')},indent=2))

if __name__=='__main__':run()
